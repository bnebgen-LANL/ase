import os
import pytest
from pathlib import Path

from ase.io import read
from ase.build import bulk
from ase.atoms import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS, BFGSLineSearch, CellAwareBFGS
from ase.filters import FrechetCellFilter

# optclasses = [BFGS] #LineSearch, BFGS, CellAwareBFGS]


def params(opt):
    run_params = dict(fmax=0.005)
    if opt is CellAwareBFGS:
        run_params.update(dict(smax=0.0005))

    return run_params


def ref_atoms():
    atoms = Atoms("Au2", positions=[[0, 0, 0], [2, 0, 0]])
    atoms.center(vacuum=10)
    # atoms = bulk("Au")
    # atoms *= 2
    # atoms.rattle(stdev=0.005, seed=1)
    atoms.calc = EMT()
    return atoms


def init_relax(opt, trajectory, restart):
    atoms = ref_atoms()
    if Path(trajectory).is_file() and os.path.getsize(trajectory) != 0:
        print("retrieving atoms from trajectory")
        atoms = read(trajectory, -1)
        atoms.calc = EMT()

    opt_relax = opt(
        FrechetCellFilter(atoms, exp_cell_factor=1.0),
        alpha=70,
        trajectory=trajectory,
        restart=restart,
    )

    return opt_relax


def fragile_optimizer(opt, trajectory, restart, kwargs):
    opt_init = init_relax(opt, trajectory, restart)
    for idx, _ in enumerate(opt_init.irun(**kwargs)):
        if idx == 5:
            break

    else:
        raise "Fragile Optimizer did not break. Check if nsteps is to large."

    #### now we restart here!
    assert opt_init.nsteps + 1 == 5 + 1
    del opt_init
    assert Path(restart).is_file() and Path(trajectory).is_file()

    print("restarting fragile relaxation")
    fragile_opt = init_relax(opt=opt, trajectory=trajectory, restart=restart)
    fragile_opt.run(**kwargs)
    # print("nsteps after break in opt.run()", opt_init.nsteps)
    print("nsteps after opt.run() in fragile", fragile_opt.nsteps)
    return (fragile_opt.nsteps + 1) + (5 + 1)


@pytest.mark.parametrize("opt", [CellAwareBFGS, BFGS, BFGSLineSearch])
def test_cellaware_bfgs_restart(testdir, opt):
    restart_filename = f"restart_{opt.__name__}.dat"
    trajectory_filename = f"{opt.__name__}.traj"
    run_kwargs = params(opt)

    single = init_relax(
        opt=opt,
        trajectory="single_" + trajectory_filename,
        restart="single_" + restart_filename,
    )
    print("running fragile relaxations")
    fragile = fragile_optimizer(
        opt=opt,
        trajectory="fragile_" + trajectory_filename,
        restart="fragile_" + restart_filename,
        kwargs=run_kwargs,
    )

    print("running single relaxations")
    single.run(**run_kwargs)

    print("Number of steps taken by each routine")
    print("Single Optimization ", single.nsteps)
    # print("Fragile Optimization", fragile)
