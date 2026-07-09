from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.mcnavstokes.elements import MCNavierStokesElements
from pyfr.solvers.mcnavstokes.inters import (MCNavierStokesIntInters,
                                             MCNavierStokesMPIInters,
                                             MCNavierStokesBaseBCInters)


class MCNavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'mcnavstokes'

    elementscls = MCNavierStokesElements
    intinterscls = MCNavierStokesIntInters
    mpiinterscls = MCNavierStokesMPIInters
    bbcinterscls = MCNavierStokesBaseBCInters
