"""Temporary file while we deprecate this locaation."""

from ase.mep import NEB as RealNEB
from ase.mep import NEBTools as RealNEBTools
from ase.mep import idpp_interpolate as realidpp_interpolate
from ase.mep import interpolate as realinterpolate
from ase.utils import deprecated


class NEB(RealNEB):
    @deprecated('Please import NEB from ase.mep, not ase.neb.')
    def __init__(self, *args, **kwargs):
        """
        .. deprecated:: 3.23.0
            Please import ``NEB`` from :mod:`ase.mep`
        """
        super().__init__(*args, **kwargs)


class NEBTools(RealNEBTools):
    @deprecated('Please import NEBTools from ase.mep, not ase.neb.')
    def __init__(self, *args, **kwargs):
        """
        .. deprecated:: 3.23.0
            Please import ``NEBTools`` from :mod:`ase.mep`
        """
        super().__init__(*args, **kwargs)


@deprecated('Please import interpolate from ase.mep, not ase.neb.')
def interpolate(*args, **kwargs):
    """
    .. deprecated:: 3.23.0
            Please import ``interpolate`` from :mod:`ase.mep`
    """
    return realinterpolate(*args, **kwargs)


@deprecated('Please import idpp_interpolate from ase.mep, not ase.neb.')
def idpp_interpolate(*args, **kwargs):
    """
    .. deprecated:: 3.23.0
            Please import ``idpp_interpolate`` from :mod:`ase.mep`
    """
    return realidpp_interpolate(*args, **kwargs)
