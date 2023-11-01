import pytest

from ase.build import bulk
from ase.calculators.emt import EMT
from ase.optimize import BFGS, BFGSLineSearch, CellAwareBFGS, RestartError
from ase.filters import FrechetCellFilter


def relax(opt):
    atoms = bulk("Au")
    atoms.calc = EMT()
    relax = opt(FrechetCellFilter(atoms, exp_cell_factor=1.0),
                alpha=70, long_output=True)
    relax.run(fmax=0.005, smax=0.00005)
    return relax.nsteps

def fragile_optimizer():
    ...

def optimizter():
    ...

@pytest.mark.parametrize("opt", [CellAwareBFGS, BFGS, BFGSLineSearch])
def test_cellaware_bfgs_restart(testdir, opt):
    fname = "tmp.dat"

    with open(fname, "w") as fd:
        fd.write("hello world\n")

    with pytest.raises(RestartError, match="Could not decode"):
        opt(Atoms(), restart=fname)
