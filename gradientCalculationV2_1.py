# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 17:02:30 2019

@author: to_reilly
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from scipy import spatial

def calcP(mm, aa, bb, kk):
    mu0 = 4 * np.pi * 10**(-7)
    dBesselk_m = special.kvp(mm, abs(kk)*aa,1)
    dBesseli_m = special.ivp(mm, abs(kk)*bb,1)
    return aa * mu0 * kk * dBesseli_m * dBesselk_m

def calcQ(mm, aa, bb, kk):
    mu0 = 4 * np.pi * 10**(-7)
    dBesselk_m = special.kvp(mm,abs(kk)*aa,1)
    Besseli_m = special.iv(mm, abs(kk) * bb)
    return mm * aa * mu0 / bb * abs(kk) / kk * Besseli_m * dBesselk_m

def calculateContour(streamF, numWires, phi, z):
    
    phi2D, z2D= np.meshgrid(phi,z)
    levels = np.linspace(np.min(streamF), np.max(streamF), numWires*2 + 4)
    levels = levels[1:-1]

    # Wire should be laid along contours between the isolines, calculate midpoint between isolines
    midPointLevels = [(levels[i]+levels[i+1])/2 for i in range(np.size(levels)-1)]
    midPointLevels = np.array(midPointLevels)[np.abs(midPointLevels) >= 1e-6] #remove zeros, account for floating point error
    
    plt.ioff()
    plt.figure()
    contWires = plt.contour(phi2D,z2D,streamF,levels = midPointLevels)
    
    return contWires

def halbachXgradient(linearLength = 140, coilRad = 135, coilLength = 350,numWires = 10, numHighOrders = 10, \
                     linearityTerm = 16, apoTerm = .05, resolution = 1e-3):
    
    gradStrength    = 1e-3
    a               = coilRad*1e-3     #gradRad
    b               = 0.001*a
    d               = linearLength*1e-3     #Linear region 
    nShape          = linearityTerm
    res             = resolution
    Zmax            = coilLength*1e-3
    N               = numHighOrders
    h               = apoTerm
    
    Nsamp           = np.int(2*Zmax/res)
    z               = np.linspace(-Zmax,Zmax,Nsamp)
    phi             = np.linspace(-1.5*np.pi,0.5*np.pi,Nsamp)
    
    
    kfov            = 1/res
    k               = np.linspace(0.00001,kfov,Nsamp).conj().T
    
    gradShape       = np.divide(z,1+(z/d)**nShape)
    gradShape_k     = np.fft.fft(gradShape.conj().T);
    
    t_k = np.exp(-2*(h*k)**2)       #apodisation term
    
    n_0    = 2*np.pi*gradShape_k*gradStrength/(calcQ(1,a,b,k)+calcP(1,a,b,k))
    streamF = 0
    
    for n in range(N+1):
        sign = (-1)**(n+1)
        scale  = 1
        for m in range(n+1):
            scale *= np.divide((calcP(2*m-1,a,b,k)-calcQ(2*m-1,a,b,k)),(calcP(2*m+1,a,b,k)+calcQ(2*m+1,a,b,k)))
        B_apo   = sign*np.fft.ifft(np.divide(np.multiply(n_0,t_k*scale),k)) 
        streamF += np.outer(B_apo,np.cos((2*n+1)*phi))
    streamF = np.real(streamF)
    

    return calculateContour(streamF, numWires, phi, z)

def halbachYgradient(linearLength = 140, coilRad = 135, coilLength = 350, numWires = 10, numHighOrders = 10, \
                     linearityTerm = 16, apoTerm = .05, resolution = 1e-3):
    
    gradStrength    = 1e-3
    a               = coilRad*1e-3     #gradRad
    b               = 0.001*a
    d               = linearLength*1e-3     #Linear region 
    Zmax            = 2*coilLength*1e-3
    N               = 0
    print("Number of higher order modes set to 0, not needed for Y and Z")
    h               = apoTerm
    
    Nsamp           = np.int(2*Zmax/resolution)
    z               = np.linspace(-Zmax,Zmax,Nsamp)
    phi             = np.linspace(-np.pi,np.pi,Nsamp)
    
    k               = np.linspace(0.00001,1/resolution,Nsamp).conj().T
    
    gradShape       = 1/(1+(z/d)**linearityTerm) - 1/(1+((z+3.5*d)/(0.5*d))**linearityTerm)-1/(1+((z-3.5*d)/(0.5*d))**linearityTerm)
    
    gradShape_k     = np.fft.fft(gradShape.conj().T);
    
    t_k             = np.exp(-2*(h*k)**2)       #apodisation term
    
    n_0             = b*2*np.pi*gradShape_k*gradStrength/(calcQ(2,a,b,k)+calcP(2,a,b,k))
    streamF         = 0
    
    for n in range(N+1):
        amp = (-1)**(n+1)
        scale = 1
        for m in range(1,n+1):
            scale *= np.divide((calcP(2*m,a,b,k)-calcQ(2*m,a,b,k)),(calcP(2*m+2,a,b,k)+calcQ(2*m+2,a,b,k)))
        B_apo   = np.fft.ifft((2/np.pi)*amp*(np.divide(1,k))*n_0*scale*t_k)
        streamF +=  np.outer(B_apo,np.sin((2*n+2)*phi))
    
    #remove sidelobes       
    d_samp=d/resolution
    streamF = streamF[int(Nsamp/2-1.5*d_samp):int(Nsamp/2+1.5*d_samp)]
    z = z[int(Nsamp/2-1.5*d_samp):int(Nsamp/2+1.5*d_samp)]
    phi2D, z2D= np.meshgrid(phi,z)

    return calculateContour(streamF.real, numWires, phi, z)

def halbachZgradient(linearLength = 140, coilRad = 135, coilLength = 350, numWires = 10, numHighOrders = 10, \
                     linearityTerm = 16, apoTerm = .05, resolution = 1e-3):
    
    gradStrength    = 1e-3
    a               = coilRad*1e-3     #gradRad
    b               = 0.001*a
    d               = linearLength*1e-3     #Linear region 
    Zmax            = 2*coilLength*1e-3
    N               = 0
    #print("Number of higher order modes set to 0, not needed for Y and Z")
    h               = apoTerm
    
    Nsamp           = np.int(2*Zmax/resolution)
    z               = np.linspace(-Zmax,Zmax,Nsamp)
    phi             = np.linspace(-np.pi,np.pi,Nsamp)+np.pi/4
    
    k               = np.linspace(0.00001,1/resolution,Nsamp).conj().T
    
    gradShape       = 1/(1+(z/d)**linearityTerm) - 1/(1+((z+3.5*d)/(0.5*d))**linearityTerm)-1/(1+((z-3.5*d)/(0.5*d))**linearityTerm)
    gradShape_k     = np.fft.fft(gradShape.conj().T);
    
    t_k             = np.exp(-2*(h*k)**2)       #apodisation term
    
    n_0             = b*2*np.pi*gradShape_k*gradStrength/(calcQ(2,a,b,k)+calcP(2,a,b,k))
    streamF         = 0
    
    for n in range(N+1):
        amp = (-1)**(n+1)
        scale = 1
        for m in range(1,n+1):
            scale *= np.divide((calcP(2*m,a,b,k)-calcQ(2*m,a,b,k)),(calcP(2*m+2,a,b,k)+calcQ(2*m+2,a,b,k)))
        B_apo   = np.fft.ifft((2/np.pi)*amp*(np.divide(1,k))*n_0*scale*t_k)
        streamF +=  np.outer(B_apo,np.cos((2*n+2)*phi))
    
    #remove sideloabs    
    d_samp=d/resolution
    streamF = streamF[int(Nsamp/2-1.5*d_samp):int(Nsamp/2+1.5*d_samp)]
    streamF[0] = 0
    streamF[-1] = 0
    z = z[int(Nsamp/2-1.5*d_samp):int(Nsamp/2+1.5*d_samp)]
    phi2D, z2D= np.meshgrid(phi,z)
    np.save("streamF_z.npy", streamF)
    np.save("phi2D.npy", streamF)
    np.save("z2D.npy", z2D)

    return calculateContour(streamF.real, numWires, phi, z)

def calculateBfield(contours, DSV, resolution, coilRad, direction):
    import numexpr as ne
    
    radius = np.float32(DSV/2)
    
    phi = np.linspace(0, 2*np.pi, int(2*np.pi*radius/resolution), dtype = np.float32)
    theta = np.linspace(0, np.pi, int(np.pi*radius/resolution), dtype = np.float32)
    phiGrid, thetaGrid = np.meshgrid(phi,theta)
    xSphere = radius*np.multiply(np.sin(thetaGrid), np.cos(phiGrid))
    ySphere = radius*np.multiply(np.sin(thetaGrid), np.sin(phiGrid))
    zSphere = radius*np.cos(thetaGrid)
    
    points = np.stack((np.ravel(xSphere),np.ravel(ySphere),np.ravel(zSphere)),axis=1)
    
    wireLevels = contours.allsegs
    gradCurrent = np.float32(contours.levels[1] - contours.levels[0])
    
    bField = np.zeros(np.shape(points)[0], dtype = np.float32)
    
    import time
    startTime = time.time()
    wireCounter = 1
    for wireLevel in wireLevels:
        for wire in wireLevel:
            print("Simulating wire %.0f of %0.f"%(wireCounter, np.size(wireLevels)))
            wire = np.array(wire, dtype = np.float32)
            wirePath3D =  np.stack((np.cos(wire[:,0])*np.float32(coilRad),np.sin(wire[:,0])*np.float32(coilRad),wire[:,1]),axis=1)
            idS = 1e-7*gradCurrent*(wirePath3D[1:,:] - wirePath3D[:-1,:])
            r = points[:,np.newaxis] - wirePath3D[:-1,:]
            r3 = np.sum(np.square(r), axis = 2)[:,:,np.newaxis]
            rNorm = r/(r3*np.sqrt(r3))
            bField += np.matmul(rNorm[:,:,2], idS[:,1]) - np.matmul(rNorm[:,:,1],idS[:,2])
            wireCounter += 1
    print("Execution time: %.2f seconds"%(time.time()-startTime))
    error = calculateError(points, bField, direction)
    return [xSphere*1e3, ySphere*1e3, zSphere*1e3], np.reshape(bField, (np.size(theta), np.size(phi))), error

def calculateError(coords, bField, direction):
    if(direction == 0):
        coordAxis = coords[:,2]
    elif(direction == 1):
        coordAxis = coords[:,1]
    else:
        coordAxis = coords[:,0]
        
    argMin = np.argmin(coordAxis)
    argMax = np.argmax(coordAxis)
    
    posRange = np.max(coordAxis) - np.min(coordAxis)
    
    bRange = bField[argMax] - bField[argMin]
    
    efficiency = bRange/posRange
    
    coordAxis[np.abs(coordAxis) < 0.01*posRange] = 'NaN'
    
    return np.nanmax((bField - efficiency*coordAxis)/(efficiency*coordAxis))
    
def exportWires(contours, coilRad, direction, conjoined):
    
    wireNum = 0
    contourDict = {}
    wireLevels = contours.allsegs
    
    if ((direction == 0) and conjoined):        #for the X gradient the center of the smallest contour is needed for joining the wires
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
            joinedContour[str(3)] = np.append(joinedContour[str(3)], contourDict[str(int(2*wireNum/numCoilSegments) + idx*2 +1)], axis = 0)
        
        
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
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    