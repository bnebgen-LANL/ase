import numpy as np
from ase.optimize.optimize import Dynamics


class ContourExploration(Dynamics):
    def __init__(self, atoms,
                 maxstep=0.5,
                 parallel_drift=0.1,
                 energy_target=None,
                 angle_limit=None,
                 force_parallel_step_scale=None,
                 remove_translation=False,
                 use_fs=True,
                 initialize_old=True, initialization_step_scale=1e-2,
                 use_target_shift=True, target_shift_previous_steps=10,
                 use_tangent_curvature=False,
                 rng=np.random,
                 force_consistent=None,
                 trajectory=None, logfile=None,
                 append_trajectory=False, loginterval=1):
        """Contour Exploration evolves the system along constant potentials
        energy contours on the potential energy surface. The method uses
        curvature based extrapolation and a potentiostat to correct for
        potential energy errors. It is similar to molecular dynamics but with a
        potentiostat rather than a thermostat and it has no timestep.
        This module was developed in conjuction with the following work:

        M. J. Waters and J. M. Rondinelli, `Contour Exploration with
        Potentiostatic Kinematics` ArXiv 2021

        Parameters:

        atoms: Atoms object
            The Atoms object to operate on. Atomic velocities are required for
            the method. If the atoms object does not contain velocities,
            random ones will be applied.

        maxstep: float
            Used to set the maximum distance an atom can move per
            iteration (default value is 0.5 Å).

        parallel_drift: float
            The fraction of the update step that is parallel to the contour but
            in a random direction. Used to break symmetries.

        energy_target: float
            The total system potential energy for that the potentiostat attepts
            to maintain. (defaults the initial potential energy)

        angle_limit: float or None
            Limits the stepsize to a maximum change of direction angle using the
            curvature. Gives a scale-free means of tuning the stepsize on the
            fly. Typically less than 30 degrees gives reasonable results but
            lower angle limits result in higher potentiostatic accuracy. Units
            of degrees. (default None)

        force_parallel_step_scale: float or None
            Scales the size of the potentiostat step. The potentiostat step is
            determined by linear extrapolation from the current potential energy
            to the target_energy with the current forces. A
            force_parallel_step_scale > 1.0 overcorrects and < 1.0
            undercorrects. By default, a simple hueristic is used to selected
            the valued based on the parallel_drift. (default None)

        remove_translation: boolean
            When True, the net momentum is removed at each step. Improves
            potentiostatic accuracy slightly for bulk systems but should not be
            used with constraints. (default False)

        use_fs: Bool
            Controls whether or not the Taylor expansion of the Frenet-Serret
            formulas for curved path extrapolation are used. Required for using
            angle_limit based step scalling. (default True)

        initialize_old: boolean
            When use_fs == True, a previous step is required for calculating the
            curvature. When initialize_old=True a small step is made to
            initialize the curvature. (default True)

        initialization_step_scale: float
            Controls the scale of the initial step as a multiple of maxstep.
            (default 1e-2)

        use_target_shift: boolean
            Enables shifting of the potentiostat target to compensate for
            systematic undercorrection or overcorrection by the potentiostat.
            Uses the average of the *target_shift_previous_steps* to prevent
            coupled occilations. (default True)

        target_shift_previous_steps: int
            The number of pevious steps to average when using use_target_shift.
            (default 10)

        use_tangent_curvature: boolean
            Use the velocity unit tangent rather than the contour normals from
            forces to compute the curvature. Usually not as accurate.
            (default False)

        rng: a random number generator
            Lets users control the random number generator for the
            parallel_drift vector. (default numpy.random)

         force_consistent: boolean
             (default None)

        trajectory: Trajectory object or str  (optional)
            Attach trajectory object.  If *trajectory* is a string a
            Trajectory will be constructed.  Default: None.

        logfile: file object or str (optional)
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.  Default: None.

        loginterval: int (optional)
            Only write a log line for every *loginterval* time steps.
            Default: 1

        append_trajectory: boolean
            Defaults to False, which causes the trajectory file to be
            overwriten each time the dynamics is restarted from scratch.
            If True, the new structures are appended to the trajectory
            file instead.
        """

        if force_parallel_step_scale is None:
            # a hureistic guess since most systems will overshoot when there is
            # drift
            FPSS = 1.1 + 0.6 * parallel_drift
        else:
            FPSS = force_parallel_step_scale

        self.atoms = atoms
        self.rng = rng
        self.remove_translation = remove_translation
        self.use_fs = use_fs
        self.force_consistent = force_consistent
        self.kappa = 0.0  # initializing a value so logging can work.
        self.use_tangent_curvature = use_tangent_curvature

        # these will be populated once self.step() is called, but could be set
        # after instantiating with ContourExploration(..., initialize_old=False)
        # to resume a previous contour trajectory
        self.T = None
        self.Told = None
        self.N = None
        self.Nold = None
        self.r_old = None
        self.r = None

        if energy_target is None:
            self.energy_target = atoms.get_potential_energy(
                force_consistent=self.force_consistent)
        else:
            self.energy_target = energy_target

        # Initizing the previous steps at the target energy slows
        # target_shifting untill the system has had
        # 'target_shift_previous_steps' steps to equilibrate and should prevent
        # occilations. These need to be initialized before the initialize_old
        # step to prevent a crash
        self.previous_energies = np.full(target_shift_previous_steps,
                                         self.energy_target)

        # we need velocities or NaNs will be produced,
        # if none are provided we make random ones
        velocities = atoms.get_velocities()
        if np.linalg.norm(velocities) < 1e-6:
            # we have to pass dimension since atoms are not yet stored
            atoms.set_velocities(self.rand_vect())

        if initialize_old:
            self.maxstep = maxstep * initialization_step_scale
            self.parallel_drift = 0.0
            # should force_parallel_step_scale be 0.0 for a better initial
            # curvature? Or at least smaller than 1.0?
            # Doesn't seem to matter much
            self.force_parallel_step_scale = 1.0
            self.use_target_shift = False
            #self.atoms = atoms
            self.step()

        self.maxstep = maxstep
        self.angle_limit = angle_limit
        self.parallel_drift = parallel_drift
        self.force_parallel_step_scale = FPSS
        self.use_target_shift = use_target_shift

        # loginterval exists for the MolecularDynamics class but not for
        # the more general Dynamics class
        Dynamics.__init__(self, atoms,
                          logfile, trajectory,  # loginterval,
                          append_trajectory=append_trajectory,
                          )

    # Required stuff for Dynamics
    def todict(self):
        return {'type': 'contour-exploration',
                'dyn-type': self.__class__.__name__,
                'stepsize': self.step_size}

    def run(self, steps=50):
        """ Call Dynamics.run and adjust max_steps """
        self.max_steps = steps + self.nsteps
        return Dynamics.run(self)

    def log(self):
        if self.logfile is not None:
            name = self.__class__.__name__
            if self.nsteps == 0:
                args = (
                    "Step",
                    "Energy_Target",
                    "Energy",
                    "Curvature",
                    "Step_Size",
                    "Energy_Deviation_per_atom")
                msg = "# %4s %15s %15s %12s %12s %15s\n" % args
                self.logfile.write(msg)
            e = self.atoms.get_potential_energy(
                force_consistent=self.force_consistent)
            dev_per_atom = (e - self.energy_target) / len(self.atoms)
            args = (
                self.nsteps,
                self.energy_target,
                e,
                self.kappa,
                self.step_size,
                dev_per_atom)
            msg = "%6d %15.6f %15.6f %12.6f %12.6f %24.9f\n" % args
            self.logfile.write(msg)

            self.logfile.flush()

    def vector_rejection(self, a, b):
        '''returns new vector that removes vector a's projection vector b'''
        aout = a - np.vdot(a, b) / np.vdot(b, b) * b
        return aout

    def normalize(self, a):
        '''Makes a unit vector out of a vector'''
        return a / np.linalg.norm(a)

    def rand_vect(self):
        vect = self.rng.rand(len(self.atoms), 3) - 0.5
        return vect

    def create_drift_unit_vector(self, N, T):
        '''Creates a random drift unit vector with no projection on N or T and
        with out a net translation so systems don't wander'''
        drift = self.rand_vect()
        drift = self.vector_rejection(drift, N)
        drift = self.vector_rejection(drift, T)
        # removes net translation, so systems don't wander
        drift = drift - drift.sum(axis=0) / len(self.atoms)
        D = self.normalize(drift)
        return D

    def _compute_update_without_fs(self, delta_s_perpendicular):

        dr_perpendicular = self.N * delta_s_perpendicular
        # Without the use of curvature there is no way to estimate the
        # limiting step size
        self.step_size = self.maxstep

        if abs(delta_s_perpendicular) < self.step_size:
            contour_step_size = np.sqrt(
                self.step_size**2 - delta_s_perpendicular**2)

            delta_s_parallel = np.sqrt(
                1 - self.parallel_drift**2) * contour_step_size
            dr_parallel = delta_s_parallel * self.T
            delta_s_drift = contour_step_size * self.parallel_drift

            D = self.create_drift_unit_vector(self.N, self.T)
            dr_drift = D * delta_s_drift

            dr = dr_parallel + dr_drift + dr_perpendicular
        else:
            delta_s_parallel = 0.0
            delta_s_drift = 0.0
            dr = self.step_size / \
                abs(delta_s_perpendicular) * dr_perpendicular
        return dr

    def _compute_update_with_fs(self, delta_s_perpendicular):
        '''Uses the Frenet–Serret formulas to perform curvature based
        extrapolation'''
        # this should keep the dr clear of the constraints
        delta_r = self.r - self.rold
        delta_s = np.linalg.norm(delta_r)
        # approximation of delta_s we use this incase an adaptive step_size
        # algo get used

        delta_T = self.T - self.Told
        delta_N = self.N - self.Nold
        dTds = delta_T / delta_s
        dNds = delta_N / delta_s
        if self.use_tangent_curvature:
            kappa = np.linalg.norm(dTds)
        else:
            # normals are better since they are fixed to the reality of
            # forces. I see smaller forces and energy errors in bulk systems
            # using the normals for kappa
            kappa = np.linalg.norm(dNds)
        self.kappa = kappa

        # on a perfect trajectory, the normal can be computed this way,
        # But the normal should always be tied to forces
        # N = dTds / kappa

        if self.angle_limit is not None:
            phi = np.pi / 180 * self.angle_limit
            self.step_size = np.sqrt(2 - 2 * np.cos(phi)) / kappa
            self.step_size = min(self.step_size, self.maxstep)

        # now we can compute a safe step
        if abs(delta_s_perpendicular) < self.step_size:
            contour_step_size = np.sqrt(
                self.step_size**2 - delta_s_perpendicular**2)
            delta_s_parallel = np.sqrt(
                1 - self.parallel_drift**2) * contour_step_size
            delta_s_drift = contour_step_size * self.parallel_drift

            N_guess = self.N + dNds * delta_s_parallel
            T_guess = self.T + dTds * delta_s_parallel
            # the extrapolation is good at keeping N_guess and T_guess
            # orthogonal but not normalized:
            N_guess = self.normalize(N_guess)
            T_guess = self.normalize(T_guess)

            dr_perpendicular = delta_s_perpendicular * (N_guess)

            dr_parallel = delta_s_parallel * self.T * (1 - (delta_s_parallel * kappa)**2 / 6.0) \
                + self.N * (kappa / 2.0) * delta_s_parallel**2

            D = self.create_drift_unit_vector(N_guess, T_guess)
            dr_drift = D * delta_s_drift

            # combine the components
            dr = dr_perpendicular + dr_parallel + dr_drift
            dr = self.step_size * self.normalize(dr)
            # because we guess our orthonormalization directions,
            # we renormalize to ensure a correct step size

        else:
            # in this case all priority goes to potentiostat terms
            delta_s_parallel = 0.0
            delta_s_drift = 0.0
            # if dr_perpendicular was used here alone, the step size could
            # be larger than the step_size which is adaptive to the
            # curvature the in this case.
            dr = self.step_size / \
                abs(delta_s_perpendicular) * dr_perpendicular

        return dr

    def step(self, f=None):
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()

        self.r = atoms.get_positions()
        current_energy = atoms.get_potential_energy(
            force_consistent=self.force_consistent)

        # Update our history of self.previous_energies to include our current
        # energy. np.roll shifts the values to keep nice sequential ordering.
        self.previous_energies = np.roll(self.previous_energies, 1)
        self.previous_energies[0] = current_energy

        if self.use_target_shift:
            target_shift = self.energy_target - np.mean(self.previous_energies)
        else:
            target_shift = 0.0

        # deltaU is the potential error that will be corrected for
        deltaU = current_energy - (self.energy_target + target_shift)

        # get the velocity vector and old kinetic energy for momentum rescaling
        velocities = atoms.get_velocities()
        KEold = atoms.get_kinetic_energy()

        # get force an and correction distance
        f_norm = np.linalg.norm(f)
        self.N = self.normalize(f)
        # can be positive or negative
        delta_s_perpendicular = (deltaU / f_norm) * \
            self.force_parallel_step_scale

        # remove velocity  projection on forces
        v_parallel = self.vector_rejection(velocities, self.N)
        self.T = self.normalize(v_parallel)

        if self.Nold is None or not self.use_fs:
            dr = self._compute_update_without_fs(delta_s_perpendicular)
        else:
            dr = self._compute_update_with_fs(delta_s_perpendicular)

        # now that dr is done, we check if there is translation
        if self.remove_translation:
            net_motion = dr.sum(axis=0) / len(atoms)
            # print(net_motion)
            dr = dr - net_motion
            dr_unit = dr / np.linalg.norm(dr)
            dr = dr_unit * self.step_size

        # save old positions before update
        self.Told = self.T
        self.Nold = self.N
        self.rold = self.r

        # if we have constraints then this will do the first part of the
        # RATTLE algorithm:
        # If we can avoid using momenta, this will be simpler.
        masses = atoms.get_masses()[:, np.newaxis]
        atoms.set_positions(self.r + dr)
        new_momenta = (atoms.get_positions() - self.r) * masses  # / self.dt

        # We need to store the momenta on the atoms before calculating
        # the forces, as in a parallel Asap calculation atoms may
        # migrate during force calculations, and the momenta need to
        # migrate along with the atoms.
        atoms.set_momenta(new_momenta, apply_constraint=False)

        # Now we get the new forces!
        f = atoms.get_forces(md=True)

        # I don't really know if removing md=True from above will break
        # compatibility with RATTLE, leaving it alone for now.
        f_constrained = atoms.get_forces()
        # but this projection needs the forces to be consistent with the
        # constraints. We have to set the new velocities perpendicular so they
        # get logged properly in the trajectory files.
        vnew = self.vector_rejection(atoms.get_velocities(), f_constrained)
        # using the md = True forces like this:
        #vnew = self.vector_rejection(atoms.get_velocities(), f)
        # will not work with constraints
        atoms.set_velocities(vnew)

        # rescaling momentum to maintain constant kinetic energy.
        KEnew = atoms.get_kinetic_energy()
        Ms = np.sqrt(KEold / KEnew)  # Ms = Momentum_scale
        atoms.set_momenta(Ms * atoms.get_momenta())

        # Normally this would be the second part of RATTLE will be done here like this:
        #  atoms.set_momenta(atoms.get_momenta() + 0.5 * self.dt * f)
        return f
