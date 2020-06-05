""" Test to check that passing unknown keywords,
which are not processed by the interface to CP2K,
raises an error.
"""

from ase.test import must_raise
from ase.calculators.calculator import CalculatorSetupError


def test_unknown_keywords(cp2k_factory):
    unknown_keywords = ['kpts', 'smearing', 'maxiter', 'convergence',
                        'spinpol', 'occupations', 'mixer']

    for keyword in unknown_keywords:
        with must_raise(CalculatorSetupError):
            cp2k_factory.calc(key=None)

    print('passed test "unknown_keywords"')
