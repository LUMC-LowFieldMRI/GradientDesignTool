from gradientCalculation import generateCoil, printGridding
from numpy import isclose
import datetime
"""
Script to perform "unit" tests of the gradientCalculation script
Checks consistency across gridding (resolution and max grid size)
as well as numerical stability of the code for various gradient 
types.
"""

designParameters = {}
designParameters['length'] = 0.140
designParameters['linearity'] = 16
designParameters['gradDir'] = "z"
designParameters['B0'] = "x"
designParameters['coilRad'] = 0.135
designParameters['fieldRad'] = 0.001*0.135
designParameters['apo'] = 0.05
designParameters['nrModes'] = 10
designParameters['nrWires'] = 15

coilTypes = [{'gradDir':'z', 'gradShape':'longitudinal'},
             {'gradDir':'x', 'gradShape':'transverse'}]

grids = [{'dz':1e-3, 'gridSize':25*designParameters['length']},
         {'dz':1e-3, 'gridSize':35*designParameters['length']},
         {'dz':1e-3, 'gridSize':45*designParameters['length']},
         {'dz':1e-1, 'gridSize':25*designParameters['length']},
         {'dz':1e-1, 'gridSize':35*designParameters['length']},
         {'dz':1e-1, 'gridSize':45*designParameters['length']}]
nrGrids = len(grids)

print("----")
print("Running Gradient Tool test at {}".format(datetime.datetime.now()))
print("----")
nrWarnings = 0
nrErrors = 0

print('')
print("Comparing {} different discretisations for consistency.".format(nrGrids))
# Test resolution and domain consistencies
iCoil = 1
for coilType in coilTypes:
    # Choose gradient type
    designParameters['type'] = coilType['gradShape']
    designParameters['gradDir'] = coilType['gradDir']
    
    # "Ground truth" for comparisons
    designParameters['resolution'] = 1e-2  # metre
    designParameters['gridSize'] = 35*designParameters['length']

    print('Testing {} coil shape in the {} direction'.format(coilType['gradShape'],coilType['gradDir']))
    print('')
    print('Base grid is')
    printGridding(designParameters)
    baseCoil = generateCoil(designParameters)
    iGrid = 1
    for discretisation in grids:
        print("Testing grid {} of {}.".format(iGrid, nrGrids))
        # Select different griddings
        designParameters['resolution'] = discretisation['dz']
        designParameters['gridSize'] = discretisation['gridSize']

        printGridding(designParameters)
        compCoil = generateCoil(designParameters)

        # Check inductance and resistance of the generated coil
        if not(isclose(compCoil['R'], baseCoil['R'],1e-3)):
            print('\tWARNING: Resistances are not equal.')
            nrWarnings += 1
        if not(isclose(compCoil['L'], baseCoil['L'],1e-3)):
            print('\tWARNING: Inductances are not equal.')
            nrWarnings += 1
        # TODO: Check wire segmets or Biot-savard consitency?
        iGrid += 1
    print('')
    iCoil += 1
        
print("----")
print("All tests finished")
print("{} warnings and {} errors observed.".format(nrWarnings, nrErrors))
print("----")


