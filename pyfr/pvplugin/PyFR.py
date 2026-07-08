# ParaView reader plugin for native PyFR files (.pyfrs solutions and
# .pyfrm meshes)
#
# Load via Tools -> Manage Plugins -> Load New, or by adding the
# containing directory to PV_PLUGIN_PATH.  Requires the PyFR package to
# be importable; if it is not on ParaView's sys.path it is located via
# the PYFR_PYTHONPATH environment variable or, failing that, relative to
# this file (which works when the plugin is shipped inside the PyFR
# source tree or an installed package).

import os
import re
import sys
import types
import weakref

import numpy as np

# Directory containing the pyfr package; embedded by 'pyfr paraview
# install' when the plugin is copied out of the source tree
PYFR_ROOT = None


def _bootstrap_pyfr():
    try:
        import pyfr  # noqa
        return
    except ImportError:
        pass

    # Resolve symlinks so a dev-mode (symlinked) install locates the
    # source tree it points into
    here = os.path.dirname(os.path.realpath(__file__))

    candidates = [PYFR_ROOT] if PYFR_ROOT else []
    if 'PYFR_PYTHONPATH' in os.environ:
        candidates.append(os.environ['PYFR_PYTHONPATH'])

    # <root>/pyfr/pvplugin/PyFR.py -> <root>
    candidates.append(os.path.dirname(os.path.dirname(here)))

    for c in candidates:
        if os.path.isfile(os.path.join(c, 'pyfr', '__init__.py')):
            sys.path.insert(0, c)
            return

    raise ImportError(
        'Unable to locate the PyFR package; set PYFR_PYTHONPATH to the '
        'directory containing it'
    )


def _shim_platformdirs():
    # pyfr.cache imports platformdirs, which ParaView does not bundle;
    # only the in-memory memoize is used here so a stub suffices
    if 'platformdirs' not in sys.modules:
        try:
            import platformdirs  # noqa
        except ImportError:
            import tempfile

            mod = types.ModuleType('platformdirs')
            mod.user_cache_dir = (
                lambda *a, **k: os.path.join(tempfile.gettempdir(), 'pyfr')
            )
            sys.modules['platformdirs'] = mod


def _import_vtk_shapes():
    # Import pyfr/writers/vtk/shapes.py directly by path; importing it
    # as pyfr.writers.vtk.shapes would execute the package __init__
    # which drags in solver/plugin dependencies not present in ParaView
    import importlib.util

    import pyfr

    path = os.path.join(os.path.dirname(pyfr.__file__), 'writers', 'vtk',
                        'shapes.py')
    spec = importlib.util.spec_from_file_location('_pyfr_vtk_shapes', path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules['_pyfr_vtk_shapes'] = mod
    spec.loader.exec_module(mod)
    return mod


_bootstrap_pyfr()
_shim_platformdirs()

import h5py

from pyfr.inifile import Inifile
from pyfr.polys import get_polybasis
from pyfr.shapes import BaseShape, proj_pts
from pyfr.util import subclass_where

get_vtk_shape = _import_vtk_shapes().get_vtk_shape

from paraview.util.vtkAlgorithm import (VTKPythonAlgorithmBase, smdomain,
                                        smhint, smproperty, smproxy)
from vtkmodules.util import numpy_support as vnp
from vtkmodules.vtkCommonCore import vtkDataArraySelection, vtkPoints
from vtkmodules.vtkCommonDataModel import (vtkCellArray, vtkCompositeDataSet,
                                           vtkMultiBlockDataSet,
                                           vtkUnstructuredGrid)
from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline


def _np_to_vtk(arr, name=None):
    varr = vnp.numpy_to_vtk(np.ascontiguousarray(arr), deep=1)
    if name:
        varr.SetName(name)
    return varr


def _selection_observer(obj):
    ref = weakref.ref(obj)

    def callback(*args):
        if (o := ref()) is not None:
            o.Modified()

    return callback


def _clean_to_grid(pts, conn, pfields):
    """Merge coincident points and average the point data at them.

    The discontinuous element-wise representation duplicates points
    along element interfaces; here these are fused by quantising the
    coordinates, with the (generally multi-valued) point data averaged
    to give a C0 continuous field.
    """
    bmin = pts.min(axis=0)
    diag = max(float(np.max(pts.max(axis=0) - bmin)), 1e-300)

    # Quantise to ~2M ticks per axis; comfortably above float32 noise
    # while packing exactly into 21 bits per axis
    q = np.rint((pts - bmin)*(2_000_000/diag)).astype(np.int64)
    key = (q[:, 0] << 42) | (q[:, 1] << 21) | q[:, 2]

    ukey, inv, counts = np.unique(key, return_inverse=True,
                                  return_counts=True)
    n = len(ukey)

    def average(arr):
        arr = arr.reshape(len(arr), -1)
        out = np.empty((n, arr.shape[1]), dtype=np.float32)
        for c in range(arr.shape[1]):
            out[:, c] = np.bincount(inv, weights=arr[:, c], minlength=n)

        return out / counts[:, None]

    pts = average(pts)
    pfields = {fn: average(arr) for fn, arr in pfields.items()}

    return pts, inv[conn].astype(np.int64), pfields


def _build_grid(pieces, pfields, cfields, clean):
    """Build a vtkUnstructuredGrid from (pts, conn, clen, ctype) pieces.

    pfields and cfields map array names to pre-concatenated point and
    cell data; with clean=True coincident points are merged and the
    point data averaged.
    """
    poff = np.cumsum([0] + [len(p[0]) for p in pieces])

    pts = np.concatenate([p[0] for p in pieces])
    conn = np.concatenate([p[1] + poff[i] for i, p in enumerate(pieces)])
    clen = np.concatenate([p[2] for p in pieces])
    ctype = np.concatenate([p[3] for p in pieces])

    if clean:
        pts, conn, pfields = _clean_to_grid(pts, conn, pfields)

    offs = np.concatenate([[0], np.cumsum(clen)])

    grid = vtkUnstructuredGrid()

    vpts = vtkPoints()
    vpts.SetData(_np_to_vtk(pts))
    grid.SetPoints(vpts)

    cells = vtkCellArray()
    cells.SetData(_np_to_vtk(offs.astype(np.int64)),
                  _np_to_vtk(conn.astype(np.int64)))
    grid.SetCells(_np_to_vtk(ctype), cells)

    pdata = grid.GetPointData()
    for fn, arr in pfields.items():
        pdata.AddArray(_np_to_vtk(arr, fn))

    cdata = grid.GetCellData()
    for fn, arr in cfields.items():
        cdata.AddArray(_np_to_vtk(arr, fn))

    return grid


def _read_inifile(f, dset):
    return Inifile(f[dset][()].decode())


# Members whose field names are not prefixed with the member name
_PREFIX_EXEMPT = {'soln', 'aux'}

# Fields which are always linear vertex data, even when the vertex and
# solution point counts coincide (p1 on simplices)
_LIN_FIELDS = {'artvisc'}


def _field_plans(dtype, nupts, nverts, ndims):
    """Classify every (member, field) in a dataset by subfield shape.

    Returns a list of (member, field, outname, kind) tuples with kind
    one of:
      'point'     -- (nupts,); interpolated with the solution basis
      'point-lin' -- (nverts,); interpolated with the linear basis
                     (e.g. artificial viscosity vertex data)
      'grad'      -- (ndims, nupts); gradient vectors
      'cell'      -- anything else; per-element cell data (e.g.
                     part-id, ef-filter, dt-cfl)
    """
    plans = []
    for m in dtype.names:
        sub = dtype[m]
        if not sub.names:
            continue

        for fn in sub.names:
            shape = sub[fn].shape
            on = fn if m in _PREFIX_EXEMPT else f'{m}-{fn}'

            if shape == (nverts,) and fn in _LIN_FIELDS:
                kind = 'point-lin'
            elif shape == (nupts,):
                kind = 'point'
            elif shape == (nverts,):
                kind = 'point-lin'
            elif shape == (ndims, nupts):
                kind = 'grad'
            else:
                kind = 'cell'

            plans.append((m, fn, on, kind))

    return plans


def _ds_plans(name, ds):
    """Order, element type, and field plans for a solution dataset."""
    if (m := re.match(r'p(\d+)-([a-z]+)$', name)) is None:
        return None

    order, etype = int(m[1]), m[2]
    shapecls = subclass_where(BaseShape, name=etype)

    nupts = (len(ds.attrs['pts']) if 'pts' in ds.attrs
             else shapecls.npts_from_order(order))
    nverts = len(shapecls.std_ele(1))

    return order, etype, _field_plans(ds.dtype, nupts, nverts,
                                      shapecls.ndims)


def _h5py_site_dirs():
    vi = sys.version_info
    return [os.environ.get('PYFR_PV_SITE'),
            os.path.expanduser(f'~/.pyfr/pv-site-{vi[0]}.{vi[1]}')]


def _switch_h5py():
    """Swap in a user-provided h5py in place of ParaView's bundled copy.

    Modern h5py wheels carry their own private libhdf5, so a copy newer
    than ParaView's can coexist with the HDF5 that VTK itself links.
    """
    global h5py

    import importlib

    for c in _h5py_site_dirs():
        if c and os.path.isdir(os.path.join(c, 'h5py')):
            for m in [m for m in sys.modules
                      if m == 'h5py' or m.startswith('h5py.')]:
                del sys.modules[m]

            sys.path.insert(0, c)
            h5py = importlib.import_module('h5py')
            return True

    return False


def _h5probe(g):
    # Resolving each object touches its header; unreadable datatype
    # messages surface here as a KeyError (note that Group.items()
    # silently maps them to None, hence the explicit indexing)
    for k in g:
        obj = g[k]
        if isinstance(obj, h5py.Group):
            _h5probe(obj)


def _h5open(fname):
    """Open a PyFR file, falling back to a user-site h5py if needed.

    Files written by PyFR with libver='latest' under HDF5 >= 2.0 cannot
    be decoded by the HDF5 1.x that ParaView bundles.  ParaView's own
    h5py is tried first; if the file proves unreadable a newer h5py is
    loaded from PYFR_PV_SITE or ~/.pyfr/pv-site-X.Y when present.
    """
    f = h5py.File(fname, 'r')
    try:
        _h5probe(f)
        return f
    except (KeyError, OSError) as e:
        err = str(e)
        f.close()
        if 'bad version number' not in err:
            raise

    if not _switch_h5py():
        vi = sys.version_info
        raise RuntimeError(
            f'{fname} was written with a newer HDF5 than ParaView '
            f'bundles.  Run "pyfr paraview install" to provision a '
            f'compatible h5py, or manually:\n'
            f'  pip download h5py --only-binary=:all: --no-deps '
            f'--python-version {vi[0]}{vi[1]} -d /tmp/pvw\n'
            f'  unzip /tmp/pvw/h5py-*.whl -d '
            f'~/.pyfr/pv-site-{vi[0]}.{vi[1]}\n'
            f'(or set PYFR_PV_SITE to a directory containing h5py)'
        )

    f = h5py.File(fname, 'r')
    _h5probe(f)
    return f


class _PyFRMesh:
    """Topology and geometry of a .pyfrm file.

    Mesh topology (face codecs, curvature, partitionings) is read
    eagerly as it is small; node locations, shape points, and the
    volume visualisation geometry are read and built on demand so that
    boundary-only views need not load the volume.
    """

    def __init__(self, fname):
        self.fname = fname

        with _h5open(fname) as f:
            self.uuid = f['mesh-uuid'][()].decode()
            self.codec = [c.decode() for c in f['codec'][...]]
            self.ndims = f['nodes'].dtype['location'].shape[0]
            self.etypes = etypes = sorted(f['eles'])

            self.faces, self.curved = {}, {}
            self.nspts, self.neles = {}, {}
            for etype in etypes:
                ds = f[f'eles/{etype}']
                self.faces[etype] = ds.fields('faces')[...]['cidx']
                self.curved[etype] = ds.fields('curved')[...]
                self.nspts[etype] = ds.dtype['nodes'].shape[0]
                self.neles[etype] = ds.shape[0]

            # Invert each partitioning into per-element rank numbers
            self.parts = {}
            for p in f.get('partitionings', []):
                ds = f[f'partitionings/{p}/eles']
                regions, eles = ds.attrs['regions'], ds[...]

                parts = {et: np.zeros(self.neles[et], dtype=np.int32)
                         for et in etypes}
                for r, einfo in enumerate(regions):
                    chunks = np.split(eles[einfo[0]:einfo[-1]],
                                      einfo[1:-1] - einfo[0])
                    for et, ch in zip(etypes, chunks):
                        parts[et][ch] = r

                self.parts[p] = parts

        self.bcs = [c[3:] for c in self.codec if c.startswith('bc/')]

        self._nodes = None
        self._spts = {}
        self._vols = {}

    def nodes(self):
        if self._nodes is None:
            with _h5open(self.fname) as f:
                self._nodes = f['nodes'].fields('location')[...]

        return self._nodes

    def spts(self, etype):
        """Shape points of all elements; (nspts, neles, ndims)."""
        if etype not in self._spts:
            with _h5open(self.fname) as f:
                enodes = f[f'eles/{etype}'].fields('nodes')[...]

            self._spts[etype] = self.nodes()[enodes].swapaxes(0, 1)

        return self._spts[etype]

    def spts_subset(self, etype, idxs):
        """Shape points of selected elements only."""
        if etype in self._spts:
            return self._spts[etype][:, idxs]

        with _h5open(self.fname) as f:
            enodes = f[f'eles/{etype}'].fields('nodes')[idxs]

        return self.nodes()[enodes].swapaxes(0, 1)

    def volume(self, divisor, ho):
        """Volume visualisation geometry at a given subdivision level.

        Returns a per-etype dict with the vis points piece, the std
        element vis points, and the sub-cell count per element.
        """
        if (key := (divisor, ho)) in self._vols:
            return self._vols[key]

        vol = {}
        for etype in self.etypes:
            pts = self.spts(etype)
            shapecls = subclass_where(BaseShape, name=etype)
            nspts, neles = self.nspts[etype], self.neles[etype]

            div = divisor or shapecls.order_from_npts(nspts)
            vshape = get_vtk_shape(etype, div)

            svpts = np.array(shapecls.std_ele(div))
            nsv = len(svpts)

            # High-order Lagrange cells (unsupported for pyramids)
            if ho and etype != 'pyr':
                svpts = svpts[vshape.nodemaps[nsv]]
                conn1 = np.arange(nsv)
                clen = np.array([nsv])
                ctype1 = np.array([vshape.vtk_ho_type], dtype=np.uint8)
            else:
                conn1 = vshape.subnodes
                clen = np.diff(vshape.subcelloffs, prepend=0)
                ctype1 = vshape.subcelltypes.astype(np.uint8)

            # Interpolate the shape points to the visualisation points
            sord = shapecls.order_from_npts(nspts)
            sbasis = get_polybasis(etype, sord, shapecls.std_ele(sord))
            mesh_op = sbasis.nodal_basis_at(svpts).astype(np.float32)

            vpts = np.tensordot(mesh_op, pts.astype(np.float32), axes=1)
            vpts = vpts.swapaxes(0, 1).reshape(-1, self.ndims)
            if self.ndims == 2:
                vpts = np.pad(vpts, [(0, 0), (0, 1)])

            # Tile the single-element connectivity across all elements
            eoff = np.arange(neles, dtype=np.int64)*nsv
            piece = (vpts, (conn1[None, :] + eoff[:, None]).ravel(),
                     np.tile(clen, neles), np.tile(ctype1, neles))

            vol[etype] = dict(piece=piece, svpts=svpts, nsub=len(ctype1))

        self._vols[key] = vol
        return vol

    def cell_arrays(self, etype, idxs=None):
        """Mesh-derived per-element cell data."""
        out = {'curved': self.curved[etype].astype(np.uint8)}
        for p, parts in self.parts.items():
            out[f'partition-{p}'] = parts[etype]

        if idxs is not None:
            out = {k: v[idxs] for k, v in out.items()}

        return out


class _BasePyFRReader(VTKPythonAlgorithmBase):
    def __init__(self, outputType='vtkUnstructuredGrid'):
        super().__init__(nInputPorts=0, nOutputPorts=1,
                         outputType=outputType)

        self._divisor = 0
        self._ho = False
        self._clean = True
        self._mesh = None

    def _mesh_key(self, fname):
        return (fname, os.path.getmtime(fname))

    def _get_mesh(self, fname):
        if self._mesh is None or self._mesh[0] != self._mesh_key(fname):
            self._mesh = (self._mesh_key(fname), _PyFRMesh(fname))

        return self._mesh[1]

    # NB: smproperty decorators are not inherited by subclasses, so the
    # concrete readers wrap these with their own decorated methods
    def SetSubdivisionLevel(self, n):
        if n != self._divisor:
            self._divisor = n
            self.Modified()

    def SetHighOrderCells(self, ho):
        if bool(ho) != self._ho:
            self._ho = bool(ho)
            self.Modified()

    def SetCleanToGrid(self, c):
        if bool(c) != self._clean:
            self._clean = bool(c)
            self.Modified()


@smproxy.reader(name='PyFRMeshReader', label='PyFR Mesh Reader',
                extensions='pyfrm', file_description='PyFR mesh files')
class PyFRMeshReader(_BasePyFRReader):
    def __init__(self):
        super().__init__()
        self._filename = None

    @smproperty.stringvector(name='FileName')
    @smdomain.filelist()
    @smhint.filechooser(extensions='pyfrm',
                        file_description='PyFR mesh files')
    def SetFileName(self, fname):
        if fname != self._filename:
            self._filename = fname
            self.Modified()

    @smproperty.intvector(name='SubdivisionLevel', default_values=0)
    @smdomain.intrange(min=0, max=12)
    def SetSubdivisionLevel(self, n):
        super().SetSubdivisionLevel(n)

    @smproperty.intvector(name='HighOrderCells', default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetHighOrderCells(self, ho):
        super().SetHighOrderCells(ho)

    @smproperty.intvector(name='CleanToGrid', default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCleanToGrid(self, c):
        super().SetCleanToGrid(c)

    def RequestData(self, request, inInfo, outInfo):
        mesh = self._get_mesh(self._filename)
        vol = mesh.volume(self._divisor, self._ho)

        # Mesh-derived cell data, repeated for each sub-cell
        cfields = {}
        for et in mesh.etypes:
            for on, arr in mesh.cell_arrays(et).items():
                cfields.setdefault(on, []).append(
                    np.repeat(arr, vol[et]['nsub'], axis=0))
        cfields = {on: np.concatenate(v) for on, v in cfields.items()}

        grid = _build_grid([vol[et]['piece'] for et in mesh.etypes],
                           {}, cfields, self._clean)

        output = vtkUnstructuredGrid.GetData(outInfo)
        output.ShallowCopy(grid)
        return 1


@smproxy.reader(name='PyFRSolutionReader', label='PyFR Solution Reader',
                extensions='pyfrs', file_description='PyFR solution files',
                support_reload=True)
class PyFRSolutionReader(_BasePyFRReader):
    def __init__(self):
        super().__init__(outputType='vtkMultiBlockDataSet')

        self._filenames = []
        self._meshfname = ''
        self._times = None

        self._arrays = vtkDataArraySelection()
        self._arrays.AddObserver('ModifiedEvent',
                                 _selection_observer(self))

        self._cells = vtkDataArraySelection()
        self._cells.AddObserver('ModifiedEvent',
                                _selection_observer(self))

        self._regions = vtkDataArraySelection()
        self._regions.AddArray('Volume')
        self._regions.AddObserver('ModifiedEvent',
                                  _selection_observer(self))

    # File series of .pyfrs files; the GUI groups numbered series and
    # passes each file via AddFileName
    @smproperty.stringvector(name='FileNames', label='File Names',
                             animateable='1', repeat_command='1',
                             clean_command='RemoveAllFileNames',
                             command='AddFileName', number_of_elements='1',
                             panel_visibility='never')
    @smdomain.filelist()
    @smhint.filechooser(extensions='pyfrs',
                        file_description='PyFR solution files')
    def AddFileName(self, fname):
        if fname and fname not in self._filenames:
            self._filenames.append(fname)
            self._times = None
            self.Modified()

    def RemoveAllFileNames(self):
        if self._filenames:
            self._filenames = []
            self._times = None
            self.Modified()

    @smproperty.stringvector(name='MeshFileName', default_values='')
    @smdomain.filelist()
    @smhint.filechooser(extensions='pyfrm',
                        file_description='PyFR mesh files')
    def SetMeshFileName(self, fname):
        fname = fname or ''
        if fname != self._meshfname:
            self._meshfname = fname
            self.Modified()

    @smproperty.dataarrayselection(name='PointArrays')
    def GetPointArraySelection(self):
        return self._arrays

    @smproperty.dataarrayselection(name='CellArrays')
    def GetCellArraySelection(self):
        return self._cells

    @smproperty.dataarrayselection(name='Regions')
    def GetRegionSelection(self):
        return self._regions

    @smproperty.intvector(name='SubdivisionLevel', default_values=0)
    @smdomain.intrange(min=0, max=12)
    def SetSubdivisionLevel(self, n):
        super().SetSubdivisionLevel(n)

    @smproperty.intvector(name='HighOrderCells', default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetHighOrderCells(self, ho):
        super().SetHighOrderCells(ho)

    @smproperty.intvector(name='CleanToGrid', default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetCleanToGrid(self, c):
        super().SetCleanToGrid(c)

    @smproperty.doublevector(name='TimestepValues', information_only='1',
                             si_class='vtkSITimeStepsProperty')
    def GetTimestepValues(self):
        return self._get_times()

    def _get_times(self):
        if self._times is None:
            times = []
            for fname in self._filenames:
                with _h5open(fname) as f:
                    stats = _read_inifile(f, 'stats')
                    times.append(stats.getfloat('solver-time-integrator',
                                                'tcurr'))

            # Sort the files by time
            srt = np.argsort(times)
            self._filenames = [self._filenames[i] for i in srt]
            self._times = [times[i] for i in srt]

        return self._times

    def _find_mesh(self, fname):
        if self._meshfname:
            return self._meshfname

        with _h5open(fname) as f:
            uuid = f['mesh-uuid'][()].decode()

        sdir = os.path.dirname(os.path.abspath(fname))
        for g in sorted(os.listdir(sdir)):
            if g.endswith('.pyfrm'):
                mpath = os.path.join(sdir, g)
                with _h5open(mpath) as f:
                    if f['mesh-uuid'][()].decode() == uuid:
                        return mpath

        raise RuntimeError(
            f'Unable to locate a .pyfrm mesh with UUID {uuid} alongside '
            f'{fname}; set the MeshFileName property'
        )

    def _pri_info(self, cfg, fnames):
        """Momentum field names if the system has known primitives."""
        system = cfg.get('solver', 'system', '')
        momf = [fn for fn in ('rhou', 'rhov', 'rhow') if fn in fnames]

        if (system in ('euler', 'navier-stokes')
                and 'rho' in fnames and 'E' in fnames and momf):
            return momf

        return None

    def _out_field_names(self, cfg, fnames):
        """Output array names for the given raw solution field names."""
        if (momf := self._pri_info(cfg, fnames)):
            out = ['Density', 'Velocity', 'Pressure', 'Energy']
            out += [fn for fn in fnames if fn not in {'rho', 'E', *momf}]
        else:
            out = list(fnames)

        return out

    def _populate_selections(self):
        if not self._filenames:
            return

        fname = self._filenames[0]
        with _h5open(fname) as f:
            cfg = _read_inifile(f, 'config')
            stats = _read_inifile(f, 'stats')
            prefix = stats.get('data', 'prefix', 'soln')

            plans = []
            for name, ds in f[prefix].items():
                if (dsp := _ds_plans(name, ds)) is not None:
                    plans = dsp[2]
                    break

        # Point fields from the solution member may map to primitives
        sfields = [fn for m, fn, on, kind in plans
                   if m == 'soln' and kind == 'point']
        for fn in self._out_field_names(cfg, sfields):
            self._arrays.AddArray(fn)

        for m, fn, on, kind in plans:
            match kind:
                case 'point' | 'point-lin' if m != 'soln':
                    self._arrays.AddArray(on)
                case 'point-lin':
                    self._arrays.AddArray(on)
                case 'grad':
                    # Gradients default to off; loaded only when enabled
                    if self._arrays.AddArray(on):
                        self._arrays.DisableArray(on)
                case 'cell':
                    self._cells.AddArray(on)

        # Boundary regions default to off
        try:
            mfname = self._find_mesh(fname)
        except RuntimeError:
            return

        with _h5open(mfname) as f:
            codec = [c.decode() for c in f['codec'][...]]
            ndims = f['nodes'].dtype['location'].shape[0]
            pnames = list(f.get('partitionings', []))

        # Mesh-derived cell data
        self._cells.AddArray('curved')
        for p in pnames:
            self._cells.AddArray(f'partition-{p}')

        if ndims == 3:
            for c in codec:
                if c.startswith('bc/') and self._regions.AddArray(c[3:]):
                    self._regions.DisableArray(c[3:])

    def RequestInformation(self, request, inInfo, outInfo):
        info = outInfo.GetInformationObject(0)
        times = self._get_times()

        info.Remove(vtkStreamingDemandDrivenPipeline.TIME_STEPS())
        info.Remove(vtkStreamingDemandDrivenPipeline.TIME_RANGE())

        if times:
            for t in times:
                info.Append(vtkStreamingDemandDrivenPipeline.TIME_STEPS(), t)
            info.Append(vtkStreamingDemandDrivenPipeline.TIME_RANGE(),
                        times[0])
            info.Append(vtkStreamingDemandDrivenPipeline.TIME_RANGE(),
                        times[-1])

        self._populate_selections()

        return 1

    def _active_file(self, outInfo):
        times = self._get_times()
        info = outInfo.GetInformationObject(0)

        idx = 0
        if (len(times) > 1 and
            info.Has(vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())):
            t = info.Get(vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())
            idx = int(np.argmin(np.abs(np.array(times) - t)))

        return self._filenames[idx], times[idx] if times else 0.0

    def _load_soln(self, fname, esubs=None):
        """Read the raw solution data from a .pyfrs file.

        If esubs is given then only the listed elements of each type
        are read; types with an empty selection are omitted entirely.
        """
        with _h5open(fname) as f:
            cfg = _read_inifile(f, 'config')
            stats = _read_inifile(f, 'stats')
            prefix = stats.get('data', 'prefix', 'soln')

            sdata, etypes = {}, []
            for name, ds in f[prefix].items():
                if (dsp := _ds_plans(name, ds)) is None:
                    continue

                order, etype, plans = dsp
                if f'{name}-idxs' in f[prefix]:
                    raise RuntimeError('Subset solution files are not '
                                       'currently supported')

                idxs = esubs.get(etype) if esubs is not None else None
                if idxs is not None and not len(idxs):
                    continue

                etypes.append(etype)

                # Read each member, deferring members which consist
                # solely of gradient fields until one is enabled
                arrs = {}
                for m in {p[0] for p in plans}:
                    fplans = [p for p in plans if p[0] == m]
                    if (all(kind == 'grad' for *_, kind in fplans)
                        and not any(self._arrays.ArrayIsEnabled(on)
                                    for _, _, on, _ in fplans)):
                        continue

                    arrs[m] = (ds.fields(m)[...] if idxs is None
                               else ds.fields(m)[idxs])

                sdata[etype] = (order, plans, arrs)

        return cfg, sdata, etypes

    def _interp_fields(self, plans, arrs, soln_op, lin_op, ndims, key,
                       fields, idxs=None):
        """Interpolate the pointwise fields of one element type.

        Results are accumulated into fields[outname][key]; gradient
        fields are only interpolated if enabled.
        """
        def sel(arr, fn):
            return arr[fn] if idxs is None else arr[fn][idxs]

        for m, arr in arrs.items():
            for kind, op in [('point', soln_op), ('point-lin', lin_op)]:
                pairs = [(fn, on) for m2, fn, on, k in plans
                         if m2 == m and k == kind]
                if not pairs:
                    continue

                sarr = np.stack([sel(arr, fn) for fn, _ in pairs], axis=-1)
                varr = np.einsum('vu,euf->evf', op, sarr,
                                 optimize=True).reshape(-1, len(pairs))
                for i, (_, on) in enumerate(pairs):
                    fields.setdefault(on, {})[key] = varr[:, i]

            pairs = [(fn, on) for m2, fn, on, k in plans
                     if m2 == m and k == 'grad'
                     and self._arrays.ArrayIsEnabled(on)]
            if pairs:
                garr = np.stack([sel(arr, fn) for fn, _ in pairs], axis=-1)
                varr = np.einsum('vu,eduf->evdf', soln_op, garr,
                                 optimize=True)
                varr = varr.reshape(-1, ndims, len(pairs))
                for i, (_, on) in enumerate(pairs):
                    fields.setdefault(on, {})[key] = varr[..., i]

    def _cell_fields(self, plans, arrs, idxs=None):
        """Per-element cell data fields of one element type."""
        out = {}
        for m, fn, on, kind in plans:
            if kind != 'cell' or m not in arrs:
                continue

            arr = arrs[m][fn] if idxs is None else arrs[m][fn][idxs]
            out[on] = arr.reshape(len(arr), -1)

        return out

    def _interp_soln(self, cfg, sdata, etypes, mesh, vol):
        """Interpolate the solution fields to the vis points.

        Returns a map from field names to per-etype (nvpts, ncomps)
        arrays.
        """
        fields = {}
        for etype in etypes:
            order, plans, arrs = sdata[etype]
            svpts = vol[etype]['svpts']

            # Solution and linear bases at the vis points
            shapecls = subclass_where(BaseShape, name=etype)
            basis = shapecls(mesh.nspts[etype], cfg)
            soln_op = basis.ubasis.nodal_basis_at(svpts).astype(np.float32)

            lbasis = get_polybasis(etype, 1, shapecls.std_ele(1))
            lin_op = lbasis.nodal_basis_at(svpts).astype(np.float32)

            self._interp_fields(plans, arrs, soln_op, lin_op, mesh.ndims,
                                etype, fields)

        return fields

    def _con_to_pri(self, cfg, fields, etypes):
        """Convert conservative fields to primitives where recognised."""
        fnames = [fn for fn in fields if not fn.startswith('grad-')]

        if (momf := self._pri_info(cfg, fnames)) is None:
            return fields

        gamma = cfg.getfloat('constants', 'gamma')
        out = {}

        for et in etypes:
            rho = fields['rho'][et]
            rhov = np.stack([fields[fn][et] for fn in momf], axis=-1)

            v = rhov / rho[:, None]
            p = (gamma - 1)*(fields['E'][et] - 0.5*np.sum(rhov*v, axis=-1))

            if v.shape[1] == 2:
                v = np.pad(v, [(0, 0), (0, 1)])

            out.setdefault('Density', {})[et] = rho
            out.setdefault('Velocity', {})[et] = v
            out.setdefault('Pressure', {})[et] = p
            out.setdefault('Energy', {})[et] = fields['E'][et]

        # Preserve any remaining (e.g. gradient) fields
        consumed = {'rho', 'E', *momf}
        for fn, data in fields.items():
            if fn not in consumed:
                out[fn] = data

        return out

    @staticmethod
    def _gather(fields, keys, sel):
        """Concatenate enabled per-key fields into named arrays."""
        return {fn: np.concatenate([data[k] for k in keys])
                for fn, data in fields.items()
                if sel.ArrayIsEnabled(fn) and all(k in data for k in keys)}

    def _boundary_grid(self, bname, cfg, sdata, mesh, esubs=None):
        """Extract a boundary surface with the solution at its points."""
        cidx = mesh.codec.index(f'bc/{bname}')
        div = self._divisor or cfg.getint('solver', 'order')

        pieces, fieldp, cellp, keys = [], {}, {}, []
        for etype in sdata:
            order, plans, arrs = sdata[etype]

            # Find the elements and face numbers on this boundary
            eoff, fidx = (mesh.faces[etype] == cidx).nonzero()
            if not len(eoff):
                continue

            shapecls = subclass_where(BaseShape, name=etype)
            basis = shapecls(mesh.nspts[etype], cfg)
            lbasis = get_polybasis(etype, 1, shapecls.std_ele(1))

            # Map global element numbers into the loaded subset
            esub = esubs.get(etype) if esubs is not None else None

            for f in np.unique(fidx):
                idxs = eoff[fidx == f]
                lidxs = idxs if esub is None else np.searchsorted(esub, idxs)
                itype, proj, _ = shapecls.faces[f]

                # Face visualisation points and connectivity
                ishapecls = subclass_where(BaseShape, name=itype)
                fsvpts = np.array(ishapecls.std_ele(div))
                vshape = get_vtk_shape(itype, div)

                if self._ho:
                    fsvpts = fsvpts[vshape.nodemaps[len(fsvpts)]]
                    conn1 = np.arange(len(fsvpts))
                    clen = np.array([len(fsvpts)])
                    ctype1 = np.array([vshape.vtk_ho_type], dtype=np.uint8)
                else:
                    conn1 = vshape.subnodes
                    clen = np.diff(vshape.subcelloffs, prepend=0)
                    ctype1 = vshape.subcelltypes.astype(np.uint8)

                # Project onto the face and build the operators
                pvpts = proj_pts(proj, fsvpts)
                mesh_op = basis.sbasis.nodal_basis_at(pvpts)
                soln_op = basis.ubasis.nodal_basis_at(pvpts)
                lin_op = lbasis.nodal_basis_at(pvpts)
                mesh_op = mesh_op.astype(np.float32)
                soln_op = soln_op.astype(np.float32)
                lin_op = lin_op.astype(np.float32)

                spts = mesh.spts_subset(etype, idxs).astype(np.float32)
                vpts = np.tensordot(mesh_op, spts, axes=1)
                vpts = vpts.swapaxes(0, 1).reshape(-1, 3)

                nsv, nsel = len(fsvpts), len(idxs)
                seoff = np.arange(nsel, dtype=np.int64)*nsv
                pieces.append((vpts,
                               (conn1[None, :] + seoff[:, None]).ravel(),
                               np.tile(clen, nsel), np.tile(ctype1, nsel)))

                # Interpolate the solution to the face points
                key = (etype, int(f))
                keys.append(key)

                self._interp_fields(plans, arrs, soln_op, lin_op,
                                    mesh.ndims, key, fieldp, idxs=lidxs)

                # Per-element cell data, repeated for each sub-cell
                cf = self._cell_fields(plans, arrs, idxs=lidxs)
                cf |= mesh.cell_arrays(etype, idxs=idxs)
                for on, arr in cf.items():
                    cellp.setdefault(on, {})[key] = np.repeat(
                        arr, len(ctype1), axis=0)

        if not pieces:
            return None

        pfields = self._gather(self._con_to_pri(cfg, fieldp, keys), keys,
                               self._arrays)
        cfields = self._gather(cellp, keys, self._cells)

        return _build_grid(pieces, pfields, cfields, self._clean)

    def RequestData(self, request, inInfo, outInfo):
        fname, time = self._active_file(outInfo)
        mesh = self._get_mesh(self._find_mesh(fname))

        vol_on = self._regions.ArrayIsEnabled('Volume')
        bnames = ([b for b in mesh.bcs if self._regions.ArrayIsEnabled(b)]
                  if mesh.ndims == 3 else [])

        # With only boundary regions enabled, restrict the solution
        # read to the elements on the selected boundaries
        esubs = None
        if not vol_on:
            cidxs = [mesh.codec.index(f'bc/{b}') for b in bnames]
            esubs = {et: np.unique(np.isin(mesh.faces[et],
                                           cidxs).nonzero()[0])
                     for et in mesh.etypes}

        cfg, sdata, etypes = self._load_soln(fname, esubs)

        # Default the subdivision level to the solver order
        div = self._divisor or cfg.getint('solver', 'order')

        mb = vtkMultiBlockDataSet()
        nblk = 0

        if vol_on:
            vol = mesh.volume(div, self._ho)

            fields = self._interp_soln(cfg, sdata, etypes, mesh, vol)
            fields = self._con_to_pri(cfg, fields, etypes)
            pfields = self._gather(fields, etypes, self._arrays)

            # Per-element cell data, repeated for each sub-cell
            cfields = {}
            for et in etypes:
                order, plans, arrs = sdata[et]
                cf = self._cell_fields(plans, arrs)
                cf |= mesh.cell_arrays(et)
                for on, arr in cf.items():
                    cfields.setdefault(on, {})[et] = np.repeat(
                        arr, vol[et]['nsub'], axis=0)
            cfields = self._gather(cfields, etypes, self._cells)

            grid = _build_grid([vol[et]['piece'] for et in etypes],
                               pfields, cfields, self._clean)

            mb.SetBlock(nblk, grid)
            mb.GetMetaData(nblk).Set(vtkCompositeDataSet.NAME(), 'Volume')
            nblk += 1

        for bname in bnames:
            grid = self._boundary_grid(bname, cfg, sdata, mesh, esubs)
            if grid is not None:
                mb.SetBlock(nblk, grid)
                mb.GetMetaData(nblk).Set(vtkCompositeDataSet.NAME(), bname)
                nblk += 1

        mb.GetInformation().Set(vtkMultiBlockDataSet.DATA_TIME_STEP(), time)

        output = vtkMultiBlockDataSet.GetData(outInfo)
        output.ShallowCopy(mb)
        return 1
