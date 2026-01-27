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


