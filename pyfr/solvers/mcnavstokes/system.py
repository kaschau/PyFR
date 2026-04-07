from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.mcnavstokes.elements import MCNavierStokesElements
from pyfr.solvers.mcnavstokes.inters import (MCNavierStokesBaseBCInters,
                                            MCNavierStokesIntInters,
                                            MCNavierStokesMPIInters)


class MCNavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'mcnavier-stokes'

    elementscls = MCNavierStokesElements
    intinterscls = MCNavierStokesIntInters
    mpiinterscls = MCNavierStokesMPIInters
    bbcinterscls = MCNavierStokesBaseBCInters