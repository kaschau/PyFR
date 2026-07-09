from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.mcnavstokes.elements import MCNavierStokesElements
from pyfr.solvers.mcnavstokes.inters import (MCNavierStokesIntInters,
                                             MCNavierStokesMPIInters,
                                             MCNavierStokesBaseBCInters)


class MCNavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'mcnavstokes'
    ef_solver = 'mceuler'

    elementscls = MCNavierStokesElements
    intinterscls = MCNavierStokesIntInters
    mpiinterscls = MCNavierStokesMPIInters
    bbcinterscls = MCNavierStokesBaseBCInters
