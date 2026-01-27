# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    renderBeamProfiles=ON, renderShellThickness=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
listKey = mdb.models.keys()
print(listKey)
if listKey[0] == 'Model-1':
	keyValue = 1
else:
	keyValue = 0
p = mdb.models[listKey[keyValue]].parts['L4-L5-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap=session.viewports['Viewport: 1'].colorMappings['Section']
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap = session.viewports['Viewport: 1'].colorMappings['Section']
cmap.updateOverrides(overrides={'Section-1-AF':(True, '#FF00BA', 'Default', 
    '#FF00BA'), 'Section-2-AF-Z2':(True, '#8510AA', 'Default', '#8510AA'), 
    'Section-3-AF-Z1': (True, '#FFB3D2', 'Default', '#FFB3D2'), 
    'Section-4-NP-Z5': (True, '#FF0000', 'Default', '#FF0000'), 
    'Section-5-NP-Z4': (True, '#FF7F00', 'Default', '#FF7F00'), 
    'Section-6-NP-Z3': (True, '#FFB200', 'Default', '#FFB200'), 
    'Section-7-NP-Z2': (True, '#FFD700', 'Default', '#FFD700'), 
    'Section-8-NP-Z1': (True, '#FFFF00', 'Default', '#FFFF00'), 
    'Section-9-NP': (True, '#0095D0', 'Default', '#0095D0'), 
    'Section-10-CEP': (True, '#00B094', 'Default', '#00B094'), 
    'Section-11-REBARS_PI': (True, '#0000FF', 'Default', '#0000FF'), 
    'Section-12-REBARS_AI': (True, '#00FFFF', 'Default', '#00FFFF'), 
    'Section-13-REBARS_PO': (True, '#800080', 'Default', '#800080'), 
    'Section-14-REBARS_AO': (True, '#C36F80', 'Default', '#C36F80'), 
    'Section-15-BEPTOP': (True, '#C8B89B', 'Default', '#C8B89B'), 
    'Section-16-BEPBOTTOM': (True, '#C8B89B', 'Default', '#C8B89B')})
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap = session.viewports['Viewport: 1'].colorMappings['Section']
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.graphicsOptions.setValues(backgroundStyle=SOLID, 
    backgroundColor='#FFFFFF', translucencyMode=2)
session.viewports['Viewport: 1'].partDisplay.setValues(renderBeamProfiles=ON, 
    renderShellThickness=ON)
e = p.elements
elements = e.getSequenceFromMask(mask=(
    '[#ffffffff:18 #0:18 #ffffffff:12 #0:12 #ffffffff:18 #0:18 #ffffffff:12', 
    ' #0:12 #ffffffff #f #0 #ff000000 #ffffffff:2 #0', 
    ' #f0000000 #ffffffff:2 #f #0 #fffff000 #ffffffff:2 #0:2', 
    ' #ffffff00 #0:9 #ffffffff:21 #0:24 #ffffffff:20 #0:8 #ffffffff:17', 
    ' #ffff #0:17 #fffff #f0000000 #ffffffff #f #fffff000', 
    ' #fffff #fffff00 #0 #fffffff0 #fff #fffff #f0000000', 
    ' #ffff #fffffff0 #fff #fff00000 #f00000ff #ffffffff:13 #0:13', 
    ' #fffffff0 #fff #fff00000 #fffffff #ffff0000 #f #fffff000', 
    ' #fffff #fffff00 #ffff0000 #f #fffff000 #fff00000 #fffffff', 
    ' #ffff0000 #ff00000f #fff #fffff #f0000000 #ffffffff #f', 
    ' #fffff000 #f0f0f00f #f00f0f0f #f00ff00f:2 #f0f00f0f #f00f0f0f #f0f00f0f', 
    ' #f0f0f00f #ff00f0f #f00ff0f0 #f00f0ff0 #f00ff00f:2 #f0ff00f #f00ff0f0', 
    ' #fff00ff0 #f0f0f00 #f0ff00f #ff0f0f0 #f0f0f0f0 #ffff:4 #ff00ff', 
    ' #ffff #ff00ff #ff00ff00:3 #ff00ff #ffff00 #ff00ff00 #ffff:4', 
    ' #ff00ff #ffff #ff00ff #ff00ff00:3 #ff00ff:2 #ff00ff00 #ff00ff:2', 
    ' #ffff #ff00ff #ff00ff00 #ff00ff #ffff00 #ff00ff:4 #ffff:2', 
    ' #ff00ff:2 #ffff #ff00ff #ff00ff00 #ff00ff:5 #ffff #ff00ff:4', 
    ' #ffff #ff00ff:2 #ff00ff00:2 #ff00ff:8 #ffff:2 #ff00ff #ff00ff00:2', 
    ' #ff00ff:6 #ffff:3 #ffff0000 #ff00ff #ff00ff00 #ff00ff:7 #ffff:3', 
    ' #ffff0000 #ff00ff #ff00ff00 #ff00ff:7 #ffffff #ffff0000 #ff', 
    ' #ffffff #ffff0000 #ff #ffffff #ffff0000 #ff #ffffff', 
    ' #ffff0000 #ff #ffffff #ffff0000 #ff #ff0000ff #ff00ff:2', 
    ' #ffffffff #ff #0 #ff000000 #ffffffff:2 #ff #ffff0000', 
    ' #ffffff #0 #ffffffff #ff #ffff0000 #ffffff #0', 
    ' #ffff:10 #ff0000ff #ffff00 #ffffffff #ff #0 #ff000000', 
    ' #ffffffff #0 #ffffff00 #ffffffff #ffff #ff000000 #ff0000ff:10', 
    ' #ffffff #0 #ff000000 #ffffff #ffff00 #ffffff00 #ff0000ff:3', 
    ' #ff #ffff00 #ffffff00 #ff0000ff #ff #ffff00 #ffffff00', 
    ' #ffff #ffffff00 #ff000000 #ffff #ffffff00 #ff000000 #ffff', 
    ' #ffffff00 #ff000000 #ffff #ffffff00 #ff000000 #ffff #ffffff00', 
    ' #ff000000 #ff0000 #ffff00ff #ffffffff #0:2 #ffff0000 #ffffff', 
    ' #0 #ffffffff #ff #ffff0000 #ffffffff #ff0000ff:10 #ff0000', 
    ' #ff00ffff #ffff0000 #ffffffff #0:2 #ffffc000 #3fffff #f0000000', 
    ' #ffffffff:2 #3ff #0 #ff000000 #ffffffff #0 #ffffff00', 
    ' #ffffffff #ffff0000:7 #0 #ffff #ffffffff #ffff0000:15 #ff300000', 
    ' #ffffffff:18 #0:18 #ffffffff:12 #0:12 #ffffffff:18 #0:18 #ffffffff:12', 
    ' #0:12 #ffffffff #f #0 #ff000000 #ffffffff:2 #0', 
    ' #f0000000 #ffffffff:2 #f #0 #fffff000 #ffffffff:2 #0:2', 
    ' #ffffff00 #0:9 #ffffffff:21 #0:24 #ffffffff:20 #0:8 #ffffffff:3', 
    ' #ffff #0:3 #f00ff00f #ff00f0f #ff0f00f #fffff0f0 #ffffffff:2', 
    ' #0:2 #ff00000 #f00f0ff0 #f00f0f0f #f0f0ff0 #f00ff00f #ffffffff:2', 
    ' #ff #0 #ffff0000 #ffffffff #0 #fff0000 #0', 
    ' #ffffff #f0000000 #ffffffff:2 #ffff #0 #ff000fff #fff000f', 
    ' #fff00 #ff000fff #fff000f #fff00 #ff000fff #f00f000f #f0f0f0f', 
    ' #fffff #f0000000 #ffffffff #ff00000f #fff #fffff #fffff00', 
    ' #ff0000 #ff00ff:4 #f00f00ff #ffff0ff0 #f #fffff000 #fff00000', 
    ' #fffffff #ff0000 #ff00ff:4 #ffff00ff #0 #ff00ffff #ffff00', 
    ' #ff00ff #ff0000ff #ffff00 #ff0000ff #ffffff00 #ffffffff #ffffff', 
    ' #0:2 #ffffffff #ffff #0 #ffffffff #f #0', 
    ' #fffff000 #ffffffff:2 #0 #fff0000 #fff00 #ff000fff #fff000f', 
    ' #fff00 #ff000fff #fff000f #fff00 #f0f0f00f #fffff #f0000000', 
    ' #ffff #fffff0 #fffff000 #ff00ff:5 #ff0f00f #fffff00f #f', 
    ' #3ffffc00 #ffff0000 #3fffff #f0000000 #ffff #fffffff0 #ff00ff00:3', 
    ' #ff00 #ffff00ff #ff00ff00:7 #3000ff00 ]', ), )
p.deleteElement(elements=elements, deleteUnreferencedNodes=ON)


