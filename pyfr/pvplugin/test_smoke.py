# Headless smoke test for the PyFR ParaView reader plugin
# Usage: pvpython test_smoke.py <plugin.py> <soln.pyfrs> [more.pyfrs ...]
import sys

from paraview.simple import LoadPlugin, UpdatePipeline

plugin, *files = sys.argv[1:]
LoadPlugin(plugin, remote=False, ns=globals())

r = PyFRSolutionReader(registrationName='soln', FileNames=files)  # noqa
r.UpdatePipelineInformation()

tv = r.TimestepValues if hasattr(r, 'TimestepValues') else None
print('timesteps:', list(tv) if hasattr(tv, '__len__') else tv)
print('point arrays available:', list(r.PointArrays.Available))
print('point arrays enabled:', list(r.PointArrays))
print('cell arrays available:', list(r.CellArrays.Available))
print('regions available:', list(r.Regions.Available))
print('regions enabled:', list(r.Regions))

UpdatePipeline(proxy=r)


def report(tag):
    di = r.GetDataInformation()
    print(f'[{tag}] npoints={di.GetNumberOfPoints()} '
          f'ncells={di.GetNumberOfCells()}')

    if (hier := di.GetHierarchy()):
        blocks = [hier.GetNodeName(hier.GetChild(0, i))
                  for i in range(hier.GetNumberOfChildren(0))]
        print(f'  blocks: {blocks}')

    for kind, ddi in [('point', di.GetPointDataInformation()),
                      ('cell', di.GetCellDataInformation())]:
        for i in range(ddi.GetNumberOfArrays()):
            ai = ddi.GetArrayInformation(i)
            ncomp = ai.GetNumberOfComponents()
            rng = [ai.GetComponentRange(c) for c in range(ncomp)]
            print(f'  {kind} array {ai.GetName()}: ncomp={ncomp} '
                  f'range={rng}')


report('default')

# Enable everything
r.PointArrays = list(r.PointArrays.Available)
r.Regions = list(r.Regions.Available)
UpdatePipeline(proxy=r)
report('all regions + all arrays')

# Boundaries only, first point array only
r.PointArrays = list(r.PointArrays.Available)[:1]
r.Regions = [n for n in r.Regions.Available if n != 'Volume']
UpdatePipeline(proxy=r)
report('boundaries only, first array only')

print('SMOKE OK')
