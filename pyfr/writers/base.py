from pyfr.readers.native import NativeReader
from pyfr.util import subclass_where


class BaseWriter:
    def __init__(self, meshf, solnf):
        from pyfr.solvers.base import BaseSystem

        # Load the mesh and solution files
        self.reader = NativeReader(meshf, construct_con=False)
        self.mesh, self.soln = self.reader.load_subset_mesh_soln(solnf)

        # Load the configuration and stats files
        self.cfg = self.soln['config']
        self.stats = self.soln['stats']

        # Data file prefix
        self.dataprefix = self.stats.get('data', 'prefix')

        # System and elements classes
        self.systemscls = subclass_where(
            BaseSystem, name=self.cfg.get('solver', 'system')
        )
        self.elementscls = self.systemscls.elementscls

        # Dimensions
        self.ndims = self.mesh.ndims
