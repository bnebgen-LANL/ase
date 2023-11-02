import pytest
import numpy as np
from math import sqrt
from os.path import getsize
from pathlib import Path

from ase.io import read
from ase.io.trajectory import Trajectory
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.optimize import BFGS, BFGSLineSearch, CellAwareBFGS
from ase.filters import FrechetCellFilter


break_now = 2


def params(opt):
    run_params = dict(fmax=0.005)
    if opt is CellAwareBFGS:
        run_params.update(dict(smax=0.0005))

    return run_params


def opt_filter_atoms(opt, trajectory, restart):
    atoms = bulk("Au")
    atoms *= 2
    atoms.rattle(stdev=0.005, seed=1)

    # if restarting, we take atoms from traj
    if Path(trajectory).is_file():
        atoms = read(trajectory, -1)

    atoms.calc = EMT()

    opt_relax = opt(
        FrechetCellFilter(atoms, exp_cell_factor=1.0),
        alpha=70,
        trajectory=trajectory,
        append_trajectory=True,
        restart=restart,
    )

    return opt_relax


def fragile_optimizer(opt, trajectory, restart, kwargs):
    fragile_init = opt_filter_atoms(opt, trajectory, restart)
    for idx, _ in enumerate(fragile_init.irun(**kwargs)):
        if idx == break_now:
            break

    else:
        raise "Fragile Optimizer did not break. Check if nsteps is to large."

    # pick up where we left off, assert we have written the files, and they
    # contain data. We check this here since these files are required in
    # order to properly restart.
    assert fragile_init.nsteps == break_now
    assert Path(restart).is_file() and Path(trajectory).is_file()
    assert all(size != 0 for size in [getsize(restart), getsize(trajectory)])

    fragile_restart = opt_filter_atoms(opt=opt, trajectory=trajectory, restart=restart)
    fragile_restart.run(**kwargs)

    return fragile_init, fragile_restart


@pytest.mark.parametrize("opt", [BFGS])
def test_optimizers_restart(testdir, opt):
    restart_filename = f"restart_{opt.__name__}.dat"
    trajectory_filename = f"{opt.__name__}.traj"
    run_kwargs = params(opt)

    # single run
    print("\nrunning single relaxations")
    single = opt_filter_atoms(
        opt=opt,
        trajectory="single_" + trajectory_filename,
        restart="single_" + restart_filename,
    )
    single.run(**run_kwargs)

    # fragile restart
    print("\nrunning fragile relaxations")
    fragile_init, fragile_restart = fragile_optimizer(
        opt=opt,
        trajectory="fragile_" + trajectory_filename,
        restart="fragile_" + restart_filename,
        kwargs=run_kwargs,
    )
    fragile_traj = read_traj("fragile_" + trajectory_filename)
    single_traj = read_traj("single_" + trajectory_filename)

    assert single.nsteps == fragile_init.nsteps + fragile_restart.nsteps

    # last step of init == first step of restart == single run break_now step
    # last step of restart == last step of single run
    for f_traj, s_traj in zip(fragile_traj, single_traj):
        # print(f_traj)
        # print(s_traj)
        for f, s in zip(f_traj, s_traj):
            # print(f_traj[f], s_traj[s])
            # print(f_traj[f] == s_traj[s])
            # print(f == s)

            assert np.allclose(f_traj[f], s_traj[s])
            assert f == s


def read_traj(file: str):
    data = []

    traj = Trajectory(file, "r")
    for idx, atoms in enumerate(traj):

        pos = atoms.get_positions()
        forces = atoms.calc.results["forces"]
        stress = atoms.calc.results["stress"]
        energy = atoms.calc.results["energy"]

        fmax = sqrt((forces**2).sum(axis=1).max())
        smax = abs(stress).max()

        tmp = {
            "step": idx,
            "energy": energy,
            "position": pos,
            "forces": forces,
            "fmax": fmax,
            "smax": smax,
        }
        data.append(tmp)

    return data
