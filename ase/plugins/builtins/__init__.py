"""
This module contains the functionality, that is built-in into ASE
and which is supposed to be used via plugins mechanism.
"""
from ase.register import register_calculator

"""
note: ase.register.register_io_format is just a wrapper
to define_io_format. We call here define_io_format to
retain the old order of the parameters (since register
begins as all register_.... functions with implementation)
"""


def ase_register():
    register_calculators()
    register_formats()


def register_calculators():
    register_calculator("ase.calculators.demon.demon.Demon")
    register_calculator("ase.calculators.kim.kim.KIM")
    register_calculator("ase.calculators.openmx.openmx.OpenMX")
    register_calculator("ase.calculators.siesta.siesta.Siesta")
    register_calculator("ase.calculators.turbomole.turbomole.Turbomole")
    register_calculator("ase.calculators.vasp.vasp.Vasp")
    register_calculator("ase.calculators.abinit.Abinit")
    register_calculator("ase.calculators.acemolecule.ACE")
    register_calculator("ase.calculators.acn.ACN")
    register_calculator("ase.calculators.aims.Aims")
    register_calculator("ase.calculators.amber.Amber")
    register_calculator("ase.calculators.castep.Castep")
    register_calculator("ase.calculators.combine_mm.CombineMM")
    register_calculator("ase.calculators.counterions.AtomicCounterIon")
    register_calculator("ase.calculators.cp2k.CP2K")
    register_calculator("ase.calculators.demonnano.DemonNano")
    register_calculator("ase.calculators.dftb.Dftb")
    register_calculator("ase.calculators.dftd3.DFTD3")
    register_calculator("ase.calculators.dmol.DMol3")
    register_calculator("ase.calculators.eam.EAM")
    register_calculator("ase.calculators.elk.ELK")
    register_calculator("ase.calculators.emt.EMT")
    register_calculator("ase.calculators.espresso.Espresso")
    register_calculator("ase.calculators.ff.ForceField")
    register_calculator("ase.calculators.gamess_us.GAMESSUS")
    register_calculator("ase.calculators.gaussian.GaussianDynamics")
    register_calculator("ase.calculators.gromacs.Gromacs")
    register_calculator("ase.calculators.lammpslib.LAMMPSlib")
    register_calculator("ase.calculators.lammpsrun.LAMMPS")
    register_calculator("ase.calculators.lj.LennardJones")
    register_calculator("ase.calculators.mopac.MOPAC")
    register_calculator("ase.calculators.morse.MorsePotential")
    register_calculator("ase.calculators.nwchem.NWChem")
    register_calculator("ase.calculators.octopus.Octopus")
    register_calculator("ase.calculators.onetep.Onetep")
    register_calculator("ase.calculators.orca.ORCA")
    register_calculator("ase.calculators.plumed.Plumed")
    register_calculator("ase.calculators.psi4.Psi4")
    register_calculator("ase.calculators.qchem.QChem")
    register_calculator("ase.calculators.qmmm.SimpleQMMM")
    register_calculator("ase.calculators.qmmm.EIQMMM")
    register_calculator("ase.calculators.qmmm.RescaledCalculator")
    register_calculator("ase.calculators.qmmm.ForceConstantCalculator")
    register_calculator("ase.calculators.qmmm.ForceQMMM")
    register_calculator("ase.calculators.tip3p.TIP3P")
    register_calculator("ase.calculators.tip4p.TIP4P")
