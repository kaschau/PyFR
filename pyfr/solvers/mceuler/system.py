from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.mceuler.elements import MCEulerElements
from pyfr.solvers.mceuler.inters import (MCEulerIntInters,
                                         MCEulerMPIInters,
                                         MCEulerBaseBCInters)


class MCEulerSystem(BaseAdvectionSystem):
    name = 'mceuler'

    elementscls = MCEulerElements
    intinterscls = MCEulerIntInters
    mpiinterscls = MCEulerMPIInters
    bbcinterscls = MCEulerBaseBCInters
