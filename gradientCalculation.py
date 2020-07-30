# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 14:33:07 2020

@author: pfuchs
"""
import numpy as np

from numpy import pi, sqrt, exp, cos, sin

from scipy import special
from scipy import integrate

from numpy.fft import fft, ifft, fftfreq, fftshift

import warnings

import matplotlib.pyplot as plt
"""
==============================================================================
    Delft University of Technology
==============================================================================

    Source Name   : gradientCalculation.py

    Author        : Patrick Fuchs
    Date          : 14/04/2020

==============================================================================
"""

"""
Script to perform gradient coil design calculations based on the target field
method. Based on the work of Bart de Vos and previous python implementation by
Tom O'Reilly of the Leiden University Medical Center

Test with (runs at bottom page)

designParameters = {}
designParameters['length'] = 0.140
designParameters['linearity'] = 16
designParameters['type'] = "longitudinal"
designParameters['gradDir'] = "z"
designParameters['B0'] = "x"
designParameters['coilRad'] = 0.135
designParameters['fieldRad'] = 0.001*0.135
designParameters['apo'] = 0.05
designParameters['nrModes'] = 10
designParameters['nrWires'] = 15
coil = generateCoil(designParameters)

or

designParameters = {}
designParameters['length'] = 0.155
designParameters['linearity'] = 30
designParameters['type'] = "transverse"
designParameters['gradDir'] = "y"
designParameters['B0'] = "x"
designParameters['coilRad'] = 0.137
designParameters['fieldRad'] = 0.0001*0.135
designParameters['apo'] = 0.05
designParameters['nrModes'] = 10
designParameters['nrWires'] = 12
coil = generateCoil(designParameters)
"""

INFO = False
DEBUG = False


def generateCoil(designParameters, **kwargs):
    """
    Wrapper function to call on the targetFunction and computeCurrent functions
    to perform the whole coil design process as well as the computeWirePaths
    function to end up with discrete wire paths that approximate the current
    distribution.

    This function also computes resistance of the wire paths, efficiency and
    power efficiency of the coil as well as inductance of the given current
    distribution and outputs these.

    Input: design Parameters dictionary with fields:
                ['length']    - shape function length parameter (d) [metre]
                ['linearity'] - shape function linearity parameter (n)
                ['type']      - shape type "longitudinal", "transverse",
                                           "mexican hat", "sawtooth" etc.
                ['gradDir']   - gradient direction "x","y","z"
                ['B0']        - background field direction "x","y","z"
                ['coilRad']   - coil radius (a) [metre]
                ['fieldRad']  - field radius (b) [metre]
                ['apo']       - apodisation factor (h)
                ['nrWires']   - number of wire paths to generate
            optional keyword arguments:
                nrModes      - number of modes to compute (10)
                wireDiameter - wire diameter (for resistance) (1.5e-3) [metre]
                current      - current for inductance computation(1) [ampere]
                gradStrength - strength of gradient (scaling) [tesla]
                resolution   - discretisation step size for z axis [m]
                gridSize     - limits of the grid z in metre leads to a grid of
                               -(gridSize/2):(gridsize/2)


    Output coil dictionary with fields:
                        ['j']     - current density (spatial domain)
                        ['j_f']   - current density (frequency domain)
                        ['b']     - target magnetic flux density (along z axis)
                        ['z']     - z-axis coordinates [metre]
                        ['phi']   - angular coordinates [rad]
                        ['k']     - frequency domain coordinates [rad]
                        ['m']     - modes of the nonzero current densities
                        ['S']     - stream function output
                        ['wires'] - contour plot wire path output
                        ['R']     - coil resistance [Ohm]
                        ['L']     - coil inductance [micro Henry]
    """

    # TODO: kwargs and defaults should be checked and made consistent as well
    # as where something is defined and which script is the master / slave

   
    lengthTerm = designParameters['length']
    linearityTerm = designParameters.get('linearity',8)
    gradType = designParameters['type']
    
    resolution = designParameters.get('resolution', 1e-3)
    # Default is so much larger than linear length for fourier resolution and
    # numerical stability
    gridSize = designParameters.get('gridSize', 35*lengthTerm)

    targetField = targetFunction(lengthTerm, linearityTerm, gradType, 
                                 resolution, gridSize)

    targetDir = designParameters['gradDir']
    currRadius = designParameters['coilRad']
    
    backgroundDir = designParameters.get('B0')
    fieldRadius = designParameters.get('fieldRad', 1e-3*currRadius)
    
    apodisation = designParameters.get('apo', 0.05)
    nrModes = kwargs.get('nrModes', 10)
    gradStrength = kwargs.get('gradStrength', 1e-3)

    currentDens = computeCurrents(targetField, targetDir, backgroundDir,
                                  currRadius, fieldRadius, apodisation,
                                  nrModes, gradStrength)

    streamFun = computeStreamFunctions(currentDens)

    nrWires = designParameters['nrWires']
    wires = calculateContour(streamFun, lengthTerm, nrWires)

    wireDiameter = kwargs.get('wireDiameter', 1.5e-3)  # metre
    wireLength = calculateWireLength(wires, currRadius)
    resistance = calculateResistance(wireLength, wireDiameter)
    if INFO:
        print("Resistance is {:.4g} \u03A9.".format(resistance))

    current = kwargs.get('current', 1)  # ampere
    inductance = calculateInductance(currentDens, currRadius, fieldRadius,
                                     current)
    if INFO:
        print("Inductance is {:.4g} \u03BCH.".format(inductance*1e6))

    coilDict = {}
    coilDict['j'] = currentDens['j']
    coilDict['j_f'] = currentDens['j_f']
    coilDict['b'] = targetField['shape']
    coilDict['b_f'] = targetField['shape-f']
    coilDict['z'] = currentDens['z']
    coilDict['phi'] = currentDens['phi']
    coilDict['k'] = currentDens['k']
    coilDict['m'] = currentDens['m']
    coilDict['S'] = streamFun['S']
    coilDict['wires'] = wires
    coilDict['length'] = wireLength
    coilDict['R'] = resistance
    coilDict['L'] = inductance
    return coilDict


def targetFunction(lengthTerm, linearityTerm, gradType, resolution, gridSize,
                    **kwargs):
    """
    targetFunction(lengthTerm, linearityTerm, gradType,
                   sampPoints=[], freqPoints=[])

    Function to return the gradient shape function (Gamma) for a given gradient
    and tuning parameters. Shape can be returned both in the spatial as well as
    the frequency domain. If sample points are provided then a DFT will be used
    to compute the frequency domain gradient shape, otherwise an FFT operation
    will be performed and uniform spatial sampling is assumed.

    Input:
        tuning parameters:   lengthTerm    - d
                             linearityTerm - n
                             radius        - b
        gradient shape type: gradType      - is "transverse", "longitudinal",
                                                "sawtooth", "mexicanhat"
        resolution                         - dz [metres]
        gridSize                           - zmax [metres] 
                                             (default grid is -zmax/2 : zmax/2)
        spatial samples    : sampPoints    - z-axis samples in [metre] array
        fourier samples    : freqPoints    - k-axis samples in [rad] array

    Output:
        gradient shape dictionary: z-domain      - 'shape'
                                   sample points - 'samp'
                                   k-domain      - 'shape-f'
                                   freq. samples - 'freq'


    The fourier domain output is using a unitary fourier transform so scaled by
    the square-root number of spatial samples used, so gives the actual scaled
    discrete fourier domain values (important for inductance computation later
    on)!

    This subFunction test.pyensures that the gridding following (for computing
    the current densities is consistent with the analytic (continuous) fourier
    transform as it is implemented in discrete fashion in this code, all
    following scripts will not check gridding. Except for the dc frequency
    which is removed when the actual current densities are computed.
    """
    # Start with the gridding
    # Default gridding is -1 to 1 metre with 201 points (cm resolution)
    # computed with numpy linspace, and corresponding frequency domain can be
    # computed with fftfreq
    step = resolution  # [m]
    Range = gridSize  # [m] 
    sampPoints = kwargs.get('sampPoints', np.arange(-Range/2, Range/2, step))
    Samples = len(sampPoints)
    freqPoints = kwargs.get('freqPoints', (fftfreq(Samples, d=step)))
    # Make sure 0th frequency is just a little off 
    # Or move this into bessel function territory?
    freqPoints[freqPoints == 0] = 1e-6

    # Check if provided gridding is consistent with FFT, or whether we use DFT
    # TODO: compute check possibly with fftshift consistency? / ifftshift!
    fftFlag = True

    # TODO: remove asserts -> move to GUI or input validation from wrapper!
    # Compute gradient shape function on the spatial grid
    
    def shape(x,d,n):
        return 1/(1+(x/d)**n)

    
    
    if gradType.lower() == "transverse":
        # Generate transverse waveform from Turner and others
        #    Shape = 1/(1+(z/d)^n)
        #gz-1./(1+((z+3.5*d)./(0.5*d)).^n)-1./(1+((z-3.5*d)./(0.5*d)).^n);
       
        gradShape = shape(sampPoints,lengthTerm,linearityTerm) 
       
        #gradShape = gradShape - 0.5*shape(sampPoints-lengthTerm*3.5,lengthTerm,linearityTerm) - 0.5*shape(sampPoints+lengthTerm*3.5,lengthTerm,linearityTerm)
       
        #gradShape = np.divide(1, 1 + (sampPoints / lengthTerm)**linearityTerm)
       
        
    elif gradType.lower() == "longitudinal":
        # Generate longitudinal waveform from Turner and others
        #    Shape = z/(1+(z/d)^n)
        gradShape = np.divide(sampPoints,
                              1 + (sampPoints / lengthTerm)**linearityTerm)
         
       
    elif gradType.lower() == "MeixcanHat":
        # Generate mexican hat or ricker wavelet
        #    Shape = A * (1 - (x/a)**2) * exp(-0.5*(x/a)**2),
        # with A = 2/(sqrt(3*a)*(pi**0.25))
        gradShape = mexicanHat(sampPoints, lengthTerm, linearityTerm)
    else:
        warnings.warn("Sorry, there is no gradient shape of type "
                      + "{} implemented yet".format(gradType))
        gradShape = np.zeros(sampPoints.shape)

    # Fourier transform gradient shape (and scale! (discrete integral - DFT))
    if fftFlag:
        Samples = len(sampPoints)
        # scale by 1/N to remain consistent -> inductance computation
        gradShape_f = (fft(gradShape))*(Range/Samples)
    else:
        # perform DFT to compute gradient shape (slower!)
        # TODO: implemenet DFT
        # gradShape_f = DFT(gradShape, sampPoints, freqPoints)
        return 0

    shapeDict = {}
    shapeDict['shape'] = gradShape
    shapeDict['shape-f'] = gradShape_f
    shapeDict['samp'] = sampPoints
    shapeDict['freq'] = freqPoints

    return shapeDict


def mexicanHat(points, length, linearity):
    """
    Generate mexican hat or ricker wavelet shape defined by
       A * (1 - (xla)**n) * exp(-0.5*(x/l)**n),
    with A = 2/(sqrt(3*l)*(pi**0.25)), l the length term, and
    n the linearity exponent.
    """
    # Generate wavelet
    scaling = 2/(sqrt(3*length)*(pi**0.25)) * (1 - (points/length)**2)
    shape = exp(-0.5*(points/length)**linearity)
    return scaling*shape


def computeCurrents(targetFunction, targetDir, backgroundDir, currRadius,
                    fieldRadius, apodisation, nrModes=10, gradStrength=1e-3):
    """
    computeCurrents(targetFunction, targetDir, backgroundDir, currRadius,
                    fieldRadius, apodisation, nrModes=10, gradStrength=1)

    Function to compute the current distribution required to generate a target
    field shape in the given direction, and for the given background field
    direction. The target direction and background direction determine the
    transfer function (combination of Bessel functions) that relate the current
    density to the given target function.

    Additionally this function will perform apodisation on the (freq. domain)
    target function and remove sidelobes from the computed current distribution
    when necessary.

    Input:
        targetFunction - should be gradient shape dict. from targetFunction
        targetDir      - "x", "y", or "z" direction, (alt. transv1/2 or longit)
  
    TODOL Righthanded vs Lefthanded axis!!
        currRadius     - radius of the cylinder of the currents density (a)
        fieldRadius    - radius of the cylinder of the target fields (b) < a!

    Output: current density J_f - frequency m/k domain
            current density J - spatial z/phi domain
    """
    gradShape = targetFunction['shape-f']
    freqPoints = targetFunction['freq']
    sampPoints = targetFunction['samp']

    nrFreqSamp = len(freqPoints)
    nrPhiSamp = 256
    phiPoints = np.linspace(-1.5*pi, 0.5*pi, nrPhiSamp)

    # Check b < a or field is defined in region I c.g. inside the coil.
    assert currRadius > fieldRadius, \
        "Oh no, fields are given outside of the cylinder coil!"

    # Compute apodisation
    apodScaling = exp(-2*(apodisation*freqPoints)**2)
    
    #permeability of free space
    mu0 = 4 * np.pi * 10**(-7)


    # initialise current density (in the frequency domain) as complex zeros
    currDen_f = np.zeros((nrModes, nrFreqSamp)) + \
        1j*np.zeros((nrModes, nrFreqSamp))

    if backgroundDir == "Conventional":
        # Ordinary target field method: see Robert Turner et alii
        
        if targetDir == "z":
        
            currDen_f[0]=-1j*gradStrength / (2 * pi * mu0 * currRadius )* gradShape * apodScaling / (freqPoints * special.kvp(0, abs(freqPoints) * currRadius) * special.iv(0, abs(freqPoints) * fieldRadius))
            
           #s=ifft(1j.*Gz./(2*pi.*mu0.*a).*(gz_fft.*exp(-2*(k.*h).^2))./(k.*abs(k).*dBesselk_m.*besseli(0,abs(k).*b)));

            
            modes = [1]
        elif (targetDir == "x") or (targetDir == "y"):
            # across the bore directed gradients

            currDen_f[0]=-1*gradStrength*fieldRadius/(2*pi*currRadius*mu0) * \
                gradShape*apodScaling / ( special.kvp(1, abs(freqPoints) * \
                currRadius) * special.iv(1, abs(freqPoints) * fieldRadius))
           
           
            
            modes = [1]
        if targetDir.lower() == "y":
             currDen_f[0] = currDen_f[0]*1j

    elif backgroundDir == "Halbach":
        # Transverse target field method: see Bart de Vos et alii
        if targetDir == "z":
            # along the bore directed gradient, only contains odd modes!
            
            # First order mode 1/2 for consistency with inverse transform (x/y)
            currDen_f[0] = -(1j)*(2*pi*gradStrength*gradShape*apodScaling) * \
                                np.divide(1, (calcP(1, currRadius, fieldRadius, freqPoints) +
                                              calcQ(1, currRadius, fieldRadius, freqPoints)))
         
       
            
            modes = [1]
            # Higher order modes
            for nn in range(1, nrModes):
                currDen_f[nn] = -np.divide(calcP(2*nn-1, currRadius, fieldRadius, freqPoints) -
                                           calcQ(2*nn-1, currRadius, fieldRadius, freqPoints),
                                           calcP(2*nn+1, currRadius, fieldRadius, freqPoints) +
                                           calcQ(2*nn+1, currRadius, fieldRadius, freqPoints))*currDen_f[nn-1]
                modes.append(2*nn+1)
                
                
                
                
        elif (targetDir.lower() == "x") or (targetDir.lower() == "y"):
            # across the bore directed gradients

            # TODO: also depends on bgDir (x/y)!
            # TODO: This is defined for y, check x consistency and recalc for x (efficiently)

            # First order mode
            currDen_f[0] = -2j*(pi*gradStrength*gradShape*fieldRadius*currRadius*apodScaling) * \
                                np.divide(1, (calcP(2, currRadius, fieldRadius, freqPoints) +
                                              calcQ(2, currRadius, fieldRadius, freqPoints)))
            modes = [2]
            
            
           
            
            if targetDir.lower() == "y":
                currDen_f[0] = -currDen_f[0]*1j
            
            # Higher order modes
            for nn in range(1, nrModes):
                currDen_f[nn] = -np.divide(calcP(2*nn, currRadius, fieldRadius, freqPoints) -
                                           calcQ(2*nn, currRadius, fieldRadius, freqPoints),
                                           calcP(2*nn+2, currRadius, fieldRadius, freqPoints) +
                                           calcQ(2*nn+2, currRadius, fieldRadius, freqPoints))*currDen_f[nn-1]
                modes.append(2*nn+2)

    # Now on to the spatial domain, inverse fourier and discrete cosine trans-
    # forms. First inverse fourier is taken along the k dimension, then sum
    # over modes Inverse transform term (whether cosine or sine has to be used)
    if backgroundDir == "Conventional":
        if targetDir.lower() == "y":
            iDT = 1j*sin(np.outer(modes, phiPoints))
        elif targetDir.lower() == "x":
            iDT = cos(np.outer(modes, phiPoints))
        else:
            iDT=phiPoints
        currDen = np.zeros((nrPhiSamp, nrFreqSamp)) +\
            1j*np.zeros((nrPhiSamp, nrFreqSamp))
        
        
            # TODO add Range to this normalisation (N/zmax)
           
        
        currDen = (1/(2*pi**2))*np.outer(iDT[0],nrFreqSamp*ifft(currDen_f[0]))
      
      
        
    elif backgroundDir == "Halbach":
        if targetDir.lower() == "y":
            iDT = 1j*sin(np.outer(modes, phiPoints))
        else:
            iDT = cos(np.outer(modes, phiPoints))
    
        currDen = np.zeros((nrPhiSamp, nrFreqSamp)) +\
            1j*np.zeros((nrPhiSamp, nrFreqSamp))
        for nn in range(nrModes):
            # TODO add Range to this normalisation (N/zmax)
            currDen += (1/(2*pi**2))*np.outer(iDT[nn],
                                              nrFreqSamp*ifft(currDen_f[nn]))
            
            
    
    currentDict = {}
    currentDict['j'] = currDen  # Spatial domain
    currentDict['j_f'] = currDen_f  # Frequency domain
    currentDict['z'] = sampPoints
    currentDict['phi'] = phiPoints
    currentDict['k'] = freqPoints
    currentDict['m'] = modes
    currentDict['iDT'] = iDT
    return currentDict


def computeStreamFunctions(currentDens):
    """
    Function to compute stream functions from given current densities (in the
    frequency domain). The stream functions are computed from the integrated
    current densities, which in the frequency domain ends up being a multipli-
    cation by 1/jk.

    Since the stream functions are not required in the frequency domain but
    rather are used in the spatial domain to find the wire paths only the
    inverse transformed stream functions are returned, the frequency trans-
    formed stream functions are used directly in the computation.

    Input: currentDensity dictionary with:
                ['phi'] - phase sample points
                ['k']   - frequency sample points
                ['m']   - modes used in the frequency domain
                ['z']   - spatial samples corresponding to freq domain shape
                ['j_f'] - frequency domain current densities

    Output: streamFun dictionary with:
                ['S']   - stream function surface
                ['z']   - z-axis samples [metre]
                ['phi'] - angle samples [rad]
    """
    phi = currentDens['phi']
    k = currentDens['k']
    modes = currentDens['m']

    nrModes = len(modes)
    nrPhiSamp = len(phi)
    nrFreq = len(k)

    currDen_f = currentDens['j_f']
    iDT = currentDens['iDT']

    streams = np.zeros((nrPhiSamp, nrFreq)) + \
        1j*np.zeros((nrPhiSamp, nrFreq))
    for nn in range(nrModes):
        streams += (1/(2*pi**2)) * \
             np.outer(iDT[nn], nrFreq*ifft(np.divide(currDen_f[nn], 1j*k)))
               
    
               
    
    streamFun = {}
    streamFun['S'] = streams
    streamFun['z'] = currentDens['z']
    streamFun['phi'] = phi

    return streamFun


def calculateContour(streamFun, length, numWires=10):
    """
    Function to compute the wire paths from a given stream function. It does
    this by computing the contours through a contour plotting routine where the
    midpoints are chosen to correspond with the given number of wires.

    Input: streamFun dictionary with:
                ['S']   - stream function surface
                ['z']   - z-axis samples [metre]
                ['phi'] - angle samples [rad]
           numWires     - scalar value, number of wires ?per quadrant?

    Output: wires  - contour plot outputs as matplotlib.contour.QuadContourSet
    """

    Z, Phi = np.meshgrid(streamFun['z'], streamFun['phi'])
    S = streamFun['S'].real

    # TODO: Why +4?
   # levels = np.linspace(S.min(), S.max(), numWires*2 + 4)
    levels = np.linspace(S.min(), S.max(), numWires*2 + 4)
    levels = levels[1:-1]

    # Wire should be laid along contours between the isolines, calculate mid-
    # point between isolines
    midPointLevels = [(levels[i]+levels[i+1])/2 for i in
                      range(np.size(levels)-1)]
    # Remove zeros, account for floating point error
    midPointLevels = np.array(midPointLevels)[np.abs(midPointLevels) >= 1e-12]

    plt.ioff()
    plt.figure()
    contWires = plt.contour(Phi, Z, S, levels=midPointLevels)
    plt.close()
    return contWires


def calculateResistance(length, wireDiameter, **kwargs):
    """
    TODO: Document
    """
    resistivity = kwargs.get('resistivity', 1.68*1e-8)  # default is copper

    resistance = resistivity*length/(pi*(wireDiameter / 2)**2)
    return resistance


def calculateWireLength(wires, coilRadius):
    """
    TODO: Document and understand
    """
    length = 0
    wireLevels = wires.allsegs
    for idx, wireLevel in enumerate(wireLevels):
        for wire in range(np.size(wireLevel, 0)):
            wirePath3D = np.stack((cos(wireLevel[wire][:, 0])*coilRadius,
                                   sin(wireLevel[wire][:, 0])*coilRadius,
                                   wireLevel[wire][:, 1]), axis=1)
            temp, indices = np.unique(wirePath3D.round(decimals=6),
                                      axis=0, return_index=True)
            wirePath3D = wirePath3D[sorted(indices)]
            delta = wirePath3D[1:, :] - wirePath3D[:-1, :]
            segLength = sqrt(np.sum(np.square(delta), axis=1))
            length += np.sum(segLength)
    return length


def calculateInductance(currentDens, currRad, fieldRad, current=1):
    """ 
    Compute the inductance for a given cylindrical current distribution
    TODO: add reference and elaborate
    """
    j_f = currentDens['j_f']
    k = currentDens['k']
    m = currentDens['m']

    dBesselk = [special.kvp(m[i], abs(k)*currRad, 1) for i in range(len(m))]
    dBesseli = [special.ivp(m[i], abs(k)*fieldRad, 1) for i in range(len(m))]

    # Weakened form for numerical stability, not necessary with apodisation!
    # dBesselk = [-0.5*(special.kve(m[i]-1,abs(k)*currRad) + \
    #               special.kve(m[i]+1,abs(k)*currRad)) for i in range(len(m))]
    # dBesseli = [0.5*(special.ive(m[i]-1,abs(k)*currRad) +  \
    #               special.ive(m[i]+1,abs(k)*currRad)) for i in range(len(m))]

    mu0 = 4*pi*1e-7
    inductance = -((2*mu0*currRad**2)/current**2) * \
        (sum(integrate.simps(abs(j_f)**2 * dBesselk*dBesseli, k)))
    return inductance.real


def calcP(mm, aa, bb, kk):
    mu0 = 4 * np.pi * 10**(-7)
    dBesselk_m = special.kvp(mm, abs(kk) * aa, 1)
    dBesseli_m = special.ivp(mm, abs(kk) * bb, 1)
    return aa * mu0 * kk * dBesseli_m * dBesselk_m


def calcQ(mm, aa, bb, kk):
    mu0 = 4 * np.pi * 10**(-7)
    dBesselk_m = special.kvp(mm, abs(kk)*aa, 1)
    Besseli_m = special.iv(mm, abs(kk) * bb)
    return mm * aa * mu0 / bb * (abs(kk) / kk) * Besseli_m * dBesselk_m


#--------------------- Old functions from Tom! -----------------
def calculateBfield(contours, DSV, resolution, coilRad, direction, backgroundDir):
    """
    TODO: Write documentation for this function and its working and
    tune inputs to fit dict or struct style of rest of code
    """
    import numexpr as ne
    
    
    radius = np.float32(DSV/2)
   
    # Removed N =  2*np.pi*radius/resolution and set to 256
    # and N= 1*np.pi*radius/resolution

    # Create grid for sphere (spherical coordinates)
    phi = np.linspace(0, 2*np.pi, 128, dtype = np.float32)
    theta = np.linspace(0, np.pi, 128, dtype = np.float32)
    phiGrid, thetaGrid = np.meshgrid(phi,theta)
    xSphere = radius*np.multiply(np.sin(thetaGrid), np.cos(phiGrid))
    ySphere = radius*np.multiply(np.sin(thetaGrid), np.sin(phiGrid))
    zSphere = radius*np.cos(thetaGrid)
    
    spherePoints = np.stack((np.ravel(xSphere), 
                             np.ravel(ySphere), 
                             np.ravel(zSphere)),axis=1)
    
    wireLevels = contours.allsegs
    
    #gradCurrent = np.float32(contours.levels[1] - contours.levels[0])
    #todo: in a more stable version go back to line above
    gradCurrent=1
    
    bField = np.zeros(np.shape(spherePoints)[0], dtype = np.float32)
    
    import time
    startTime = time.time()
    wireCounter = 1
    for wireLevel in wireLevels:
        for wire in wireLevel:
            print("Simulating wire {:.0f} of {:.0f}".format(wireCounter, np.size(wireLevels)))
            wire = np.array(wire, dtype = np.float32)
            wirePath3D =  np.stack((np.cos(wire[:,0])*np.float32(coilRad),
                                    np.sin(wire[:,0])*np.float32(coilRad),
                                    wire[:,1]),axis=1)
            # Take step along the wire path, each dl pointing to the next segment
            # TODO: Why 1e-7?
            idS = 1e-7*gradCurrent*np.diff(wirePath3D,1,0)
            
            # Promote dimension of points on sphere to take r as vector pointing
            # from wire segment (start) to sphere points
            r = spherePoints[:,np.newaxis] - wirePath3D[:-1,:]
            
            # Use numexpr to evaluate the length of the vectors for speedup
            r3 = ne.evaluate("sum(r**2, axis = 2)")[:,:,np.newaxis]
            # And compute the unit vectors using numexpr to speedup array processing
            rNorm = ne.evaluate("r/r3**1.5")
            # Field is given by superposition of fields from each wire segment by biot-savard
            if backgroundDir == "Halbach": 
                bField += np.matmul(rNorm[:,:,1], idS[:,2]) - np.matmul(rNorm[:,:,2],idS[:,1])
            elif backgroundDir == "Conventional":
                bField += np.matmul(rNorm[:,:,0], idS[:,1]) - np.matmul(rNorm[:,:,1],idS[:,0])
            wireCounter += 1
    print("Execution time: {:.2f} seconds".format(time.time()-startTime))
    [error,efficiency] = calculateError(spherePoints, bField, direction, backgroundDir)
    print(efficiency)
    return [xSphere*1e3, ySphere*1e3, zSphere*1e3], \
           np.reshape(bField, (np.size(theta), np.size(phi))), efficiency,\
           error


def calculateError(coords, bField, direction, backgroundDir):
    """
    TODO: write documentation for this function and structure input/output to fit rest
    of code.
    """
    
    if backgroundDir == "Halbach":
        if(direction == 0):
            coordAxis = coords[:,2]
        elif(direction == 1):
            coordAxis = coords[:,1]
        else:
            coordAxis = coords[:,0]
    elif backgroundDir == "Conventional":
        if(direction == 0):
            coordAxis = coords[:,0]
        elif(direction == 1):
            coordAxis = coords[:,1]
        else:
            coordAxis = coords[:,2]
        
        
        
    argMin = np.argmin(coordAxis)
    argMax = np.argmax(coordAxis)
    
    posRange = np.max(coordAxis) - np.min(coordAxis)
    
    bRange = bField[argMax] - bField[argMin]
    
    efficiency= (bRange/posRange)
    
    
    coordAxis[np.abs(coordAxis) < 0.01*posRange] = 'NaN'
    
    
    return np.nanmax((bField - efficiency*coordAxis)/(efficiency*coordAxis)),efficiency


def exportWires(contours, coilRad, direction, conjoined):
    """
    TODO: write documentation for this function and structure input/output to fit rest
    of code. Maybe move this function to the GUI method?
    """   
    wireNum = 0
    contourDict = {}
    wireLevels = contours.allsegs

    #for the X gradient the center of the smallest contour is needed for joining the wires
    if ((direction == 0) and conjoined):        
        minLength = np.inf
        for wireLevel in wireLevels:
            for wire in wireLevel:
                if(np.size(wire,0) < minLength):
                    centerHeight = np.abs(np.mean(wire[:,1])*1e3)
    
    for wireLevel in wireLevels:
        for wire in wireLevel:
            wirePath3D = np.stack((np.cos(wire[:,0])*coilRad,np.sin(wire[:,0])*coilRad,wire[:,1]*1e3),axis=1)
            if(conjoined):
                gapSize = 8 #gap in which the sections are joined
                gapAngle = gapSize/coilRad      
                centerAngle = np.mean(wire[:,0])
                
                if(direction == 0):
                    mask = (np.abs(wire[:,0] - centerAngle) > gapAngle) | (np.abs(wirePath3D[:,2]) < centerHeight)
                else:
                    mask = (np.abs(wire[:,0] - centerAngle) > gapAngle) | (wirePath3D[:,2] < 0)
                
                while mask[0]:
                    mask = np.roll(mask,1)
                    wirePath3D = np.roll(wirePath3D, 1, axis = 0)
        
                contourDict[str(wireNum)] = np.stack((wirePath3D[mask, 0],wirePath3D[mask, 1],wirePath3D[mask, 2]),axis=1)
            else:
                contourDict[str(wireNum)] = wirePath3D
            wireNum += 1
    
    if(not conjoined):
        return contourDict
    else:
        
        
        #############################################
        # Join the wires with a gap in to one array #
        #############################################
        
        numCoilSegments = 4             #Number of quadrants

        joinedContour = {}
        joinedContour[str(0)] = contourDict[str(0)]
        joinedContour[str(1)] = contourDict[str(1)]
        joinedContour[str(2)] = contourDict[str(int(2*wireNum/numCoilSegments))]
        joinedContour[str(3)] = contourDict[str(int(2*wireNum/numCoilSegments)+1)]
        
        for idx in range(1,int(wireNum/numCoilSegments)):
            joinedContour[str(0)] = np.append(joinedContour[str(0)], contourDict[str(2*idx)], axis = 0)
            joinedContour[str(1)] = np.append(joinedContour[str(1)], contourDict[str(2*idx+1)], axis = 0)
            joinedContour[str(2)] = np.append(joinedContour[str(2)], contourDict[str(int(2*wireNum/numCoilSegments) + idx*2 )], axis = 0)
            joinedContour[str(3)] = np.append(joinedContour[str(3)], contourDict[str(int(2*wireNum/numCoilSegments) + idx*2 + 1)], axis = 0)
        
        
        ############################################
        # Check for consecutive identical elements #
        ############################################
        tol = 1e-5
        for key in joinedContour:
            delta = joinedContour[key][1:,:] - joinedContour[key][:-1,:]
            delta = np.sum(np.square(delta), axis = 1)
            zeroElements = delta < tol
            joinedContour[key] = np.delete(joinedContour[key],np.nonzero(zeroElements), axis = 0)
            
        return joinedContour



def printGridding(designParameters):
    print('\t{:.3g} to {:.3g}'.format(-designParameters['gridSize']/2, 
                                       designParameters['gridSize']/2))
    print('\twith increment of {:.3g}'.format(designParameters['resolution']))

