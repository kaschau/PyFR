# Headless render test for the PyFR ParaView reader plugin
# Usage: pvbatch test_render.py <plugin.py> <out.png> <soln.pyfrs> [...]
import sys

from paraview.simple import (ColorBy, GetActiveViewOrCreate,
                             GetAnimationScene, GetColorTransferFunction,
                             GetDisplayProperties, GetScalarBar, Hide,
                             LoadPlugin, Render, SaveScreenshot, Show, Slice,
                             UpdatePipeline)

plugin, outpng, *files = sys.argv[1:]
LoadPlugin(plugin, remote=False, ns=globals())

r = PyFRSolutionReader(registrationName='soln', FileNames=files)  # noqa
r.UpdatePipelineInformation()

tv = r.TimestepValues
tv = list(tv) if hasattr(tv, '__len__') else [tv]
print('timesteps:', tv)

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.AnimationTime = tv[-1]

sl = Slice(registrationName='slice', Input=r)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [5.0, 0.0, 0.0]
sl.SliceType.Normal = [0.0, 0.0, 1.0]

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1200, 800]

disp = Show(sl, view)
ColorBy(disp, ('POINTS', 'Pressure'))
disp.RescaleTransferFunctionToDataRange(True)
disp.SetScalarBarVisibility(view, True)

view.ResetCamera()
view.CameraPosition = [5.0, 0.0, 25.0]
view.CameraFocalPoint = [5.0, 0.0, 0.0]
view.CameraViewUp = [0.0, 1.0, 0.0]

Render()
SaveScreenshot(outpng, view)
print('saved', outpng)

di = sl.GetDataInformation()
print('slice npoints:', di.GetNumberOfPoints(), 'time:', view.ViewTime)
print('RENDER OK')
