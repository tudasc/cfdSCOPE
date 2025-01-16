# trace generated using paraview version 5.13.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
import argparse
import glob


# Input parameters
parser = argparse.ArgumentParser("Visualizes CFD data")
parser.add_argument("input", type=str, action="store", help="input file wildcard. ")
parser.add_argument("--frame", type=int, action="store", default=1, help="frame to render.")
parser.add_argument("--animate", action="store_true", help="save animation.")
parser.add_argument("--size", type=int, action="store", help="grid size")

argvals = parser.parse_args()

file_list = sorted(glob.glob(argvals.input))

print(f"Loaded {len(file_list)} files. First filename: '{file_list[0]}'")

frame_to_render = argvals.frame
animate = argvals.animate

grid_size = argvals.size

num_frames = len(file_list)-1


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

fieldscsv0 = CSVReader(registrationName='fields.csv.0*', FileName=file_list)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
fieldscsv0Display = Show(fieldscsv0, spreadSheetView1, 'SpreadSheetRepresentation')

# trace defaults for the display properties.
fieldscsv0Display.Assembly = ''

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=0)

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Table To Structured Grid'
tableToStructuredGrid1 = TableToStructuredGrid(registrationName='TableToStructuredGrid1', Input=fieldscsv0)
tableToStructuredGrid1.XColumn = 'p'
tableToStructuredGrid1.YColumn = 'p'
tableToStructuredGrid1.ZColumn = 'p'

# Properties modified on tableToStructuredGrid1
tableToStructuredGrid1.WholeExtent = [0, grid_size-1, 0, grid_size-1, 0, grid_size-1]
tableToStructuredGrid1.XColumn = 'x'
tableToStructuredGrid1.YColumn = 'y'
tableToStructuredGrid1.ZColumn = 'z'

# show data in view
tableToStructuredGrid1Display = Show(tableToStructuredGrid1, spreadSheetView1, 'SpreadSheetRepresentation')

# trace defaults for the display properties.
tableToStructuredGrid1Display.Assembly = ''

# hide data in view
Hide(fieldscsv0, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=tableToStructuredGrid1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'Velocity'
calculator1.Function = 'u*iHat + v*jHat + w*kHat'

# show data in view
calculator1Display = Show(calculator1, spreadSheetView1, 'SpreadSheetRepresentation')

# trace defaults for the display properties.
calculator1Display.Assembly = ''

# hide data in view
Hide(tableToStructuredGrid1, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)
calculator2.Function = ''

# Properties modified on calculator2
calculator2.ResultArrayName = 'Velocity2'
calculator2.Function = 'u*iHat + v*jHat'

# show data in view
calculator2Display = Show(calculator2, spreadSheetView1, 'SpreadSheetRepresentation')

# trace defaults for the display properties.
calculator2Display.Assembly = ''

# hide data in view
Hide(calculator1, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=calculator2)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [grid_size/2.0, grid_size/2.0, grid_size/2.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [grid_size/2.0, grid_size/2.0, grid_size/2.0]

# set active source
SetActiveSource(slice1)

# show data in view
slice1Display = Show(slice1, spreadSheetView1, 'SpreadSheetRepresentation')

# trace defaults for the display properties.
slice1Display.Assembly = ''

# set active view
SetActiveView(renderView1)

# show data in view
slice1Display_1 = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display_1.Representation = 'Surface'
slice1Display_1.ColorArrayName = [None, '']
slice1Display_1.SelectNormalArray = 'None'
slice1Display_1.SelectTangentArray = 'None'
slice1Display_1.SelectTCoordArray = 'None'
slice1Display_1.TextureTransform = 'Transform2'
slice1Display_1.OSPRayScaleArray = 'Velocity'
slice1Display_1.OSPRayScaleFunction = 'Piecewise Function'
slice1Display_1.Assembly = ''
slice1Display_1.SelectedBlockSelectors = ['']
slice1Display_1.SelectOrientationVectors = 'Velocity2'
slice1Display_1.ScaleFactor = 6.300000000000001
slice1Display_1.SelectScaleArray = 'None'
slice1Display_1.GlyphType = 'Arrow'
slice1Display_1.GlyphTableIndexArray = 'None'
slice1Display_1.GaussianRadius = 0.315
slice1Display_1.SetScaleArray = ['POINTS', 'Velocity']
slice1Display_1.ScaleTransferFunction = 'Piecewise Function'
slice1Display_1.OpacityArray = ['POINTS', 'Velocity']
slice1Display_1.OpacityTransferFunction = 'Piecewise Function'
slice1Display_1.DataAxesGrid = 'Grid Axes Representation'
slice1Display_1.PolarAxes = 'Polar Axes Representation'
slice1Display_1.SelectInputVectors = ['POINTS', 'Velocity2']
slice1Display_1.WriteLog = ''

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display_1.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display_1.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [242.55, 31.5, 31.5]
renderView1.CameraFocalPoint = [grid_size/2.0, grid_size/2.0, grid_size/2.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

renderView1.ResetActiveCameraToPositiveZ()

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

renderView1.AdjustRoll(-90.0)

renderView1.AdjustRoll(-90.0)

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# show data in view
slice1Display_1 = Show(slice1, renderView1, 'GeometryRepresentation')

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

#changing interaction mode based on data extents
renderView1.CameraPosition = [grid_size/2.0, grid_size/2.0, 242.55]
renderView1.CameraViewUp = [0.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# set scalar coloring
ColorBy(slice1Display_1, ('POINTS', 'Velocity', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
slice1Display_1.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
slice1Display_1.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Velocity'
velocityLUT = GetColorTransferFunction('Velocity')

# get opacity transfer function/opacity map for 'Velocity'
velocityPWF = GetOpacityTransferFunction('Velocity')

# get 2D transfer function for 'Velocity'
velocityTF2D = GetTransferFunction2D('Velocity')

# Rescale transfer function
velocityLUT.RescaleTransferFunction(0.0, 200)

# Rescale transfer function
velocityPWF.RescaleTransferFunction(0.0, 200)

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=slice1,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'Velocity2']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 6.300000000000001
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.GlyphType = '2D Glyph'
glyph1.ScaleFactor = 2.2

# show data in view
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = ['POINTS', 'Velocity']
glyph1Display.LookupTable = velocityLUT
glyph1Display.SelectNormalArray = 'None'
glyph1Display.SelectTangentArray = 'None'
glyph1Display.SelectTCoordArray = 'None'
glyph1Display.TextureTransform = 'Transform2'
glyph1Display.OSPRayScaleArray = 'Velocity'
glyph1Display.OSPRayScaleFunction = 'Piecewise Function'
glyph1Display.Assembly = ''
glyph1Display.SelectedBlockSelectors = ['']
glyph1Display.SelectOrientationVectors = 'Velocity2'
glyph1Display.ScaleFactor = 6.455563426017761
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 0.3227781713008881
glyph1Display.SetScaleArray = ['POINTS', 'Velocity']
glyph1Display.ScaleTransferFunction = 'Piecewise Function'
glyph1Display.OpacityArray = ['POINTS', 'Velocity']
glyph1Display.OpacityTransferFunction = 'Piecewise Function'
glyph1Display.DataAxesGrid = 'Grid Axes Representation'
glyph1Display.PolarAxes = 'Polar Axes Representation'
glyph1Display.SelectInputVectors = ['POINTS', 'Velocity2']
glyph1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-51.6118, 0.0, 0.5, 0.0, 199.875, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-51.6118, 0.0, 0.5, 0.0, 199.875, 1.0, 0.5, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

print("Glyphs added!")


## Final camera

renderView1.ResetActiveCameraToPositiveZ()

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

renderView1.AdjustRoll(-90.0)

renderView1.AdjustRoll(-90.0)

# layout/tab size in pixels
layout1 = GetLayout()
layout1.SetSize(2200, 900)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [29.38346951893936, 28.625569272380417, -221.46253025031308]
renderView1.CameraFocalPoint = [29.38346951893936, 28.625569272380417, 31.5]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 36.95696649796621

renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [grid_size/2.0, grid_size/2.0, -141.27681186415748]
renderView1.CameraFocalPoint = [grid_size/2.0, grid_size/2.0, grid_size/2.0]
renderView1.CameraViewUp = [4.440892098500626e-16, -1.0, 0.0]
renderView1.CameraParallelScale = 38.581574850542545

# get color transfer function/color map for 'Velocity'
velocityLUT = GetColorTransferFunction('Velocity')

# get color legend/bar for velocityLUT in view renderView1
velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)

# change scalar bar placement
velocityLUTColorBar.WindowLocation = 'Any Location'
velocityLUTColorBar.Position = [0.8661909009812666, 0.09999999999999998]
velocityLUTColorBar.ScalarBarLength = 0.3299999999999998

## Set frame

animationScene1.AnimationTime = frame_to_render

print(f"Saving frame {frame_to_render}.")

SaveScreenshot(f"frame_{frame_to_render}.png", viewOrLayout=renderView1, location=16, ImageResolution=[761, 900])

if animate:
    # save animation
    print(f"Saving animation with {num_frames} frames")
    SaveAnimation(filename='anim.avi', viewOrLayout=renderView1, location=16, ImageResolution=[1120, 900],
FrameRate=4, FrameWindow=[0, num_frames])
    print("Done!")
    
Interact()


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
