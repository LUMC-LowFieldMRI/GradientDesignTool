# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 12:57:51 2019

@author: to_reilly
"""

import numpy as np
import tkinter as tk
from tkinter import ttk
import os
import gradientCalculationV2_1 as gradCalc

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)
from matplotlib.figure import Figure
from matplotlib.pyplot import close as closePlots
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class PopUpTextMessage(tk.Toplevel):
    def __init__(self, parent, textmsg):
        super(PopUpTextMessage, self).__init__(parent)

        self.wm_title('!')
        label = ttk.Label(self, text=textmsg)
        label.pack(side='top', fill='x', pady=10, padx=10)

        b1 = ttk.Button(self, text='Ok', command=self.cleanup)
        b1.pack(pady=10, padx=10)

    def cleanup(self):
        self.destroy()

class GDT_GUI(tk.Tk):
    
    def __init__ (self, *args, **kwargs):
        super(GDT_GUI, self).__init__(*args, **kwargs)
        # set name and icon
        # super().iconbitmap(self, default=None)
        super().wm_title('Gradient design tool')

        inputFrame = ttk.Frame(self,borderwidth=1)
        inputFrame.grid(row=0,column=0,columnspan=1, sticky=tk.W, padx=10, pady=10)
        
        inputDescriptionLabel = ttk.Label(inputFrame, text='Design parameters', font='Helvetica 12 bold')
        inputDescriptionLabel.grid(row=0,column=0, columnspan = 2, sticky=tk.W)
        
        #create variables for storing inputs
        self.linearLength       = tk.DoubleVar(inputFrame, value = 140)
        self.coilRadius         = tk.DoubleVar(inputFrame, value = 135)
        self.coilLength         = tk.DoubleVar(inputFrame, value = 400)
        self.numWires           = tk.IntVar(inputFrame, value = 14)
        self.wireDiameter       = tk.DoubleVar(inputFrame, value = 1.5)
        self.apodisationTerm    = tk.DoubleVar(inputFrame, value = 0.05)
        self.linearityOrder     = tk.IntVar(inputFrame, value = 16)
        self.higherOrderTerms   = tk.IntVar(inputFrame, value = 10)
        self.simulationSize     = tk.DoubleVar(inputFrame, value= 140)
        self.simulationRes      = tk.DoubleVar(inputFrame, value= 5)
        self.DSV                = tk.DoubleVar(inputFrame, value = 140)
        
        directionLabel = ttk.Label(inputFrame, text='Gradient direction')
        directionLabel.grid(row=1,column=0, pady=5,sticky=tk.W)
        self.directionInput = ttk.Combobox(inputFrame, values=["X", "Y", "Z"], width = 7, justify=tk.RIGHT, state = 'readonly')
        self.directionInput.current(0)
        self.directionInput.grid(row=1,column=1, padx = 20)
        
        targetRegionLabel = ttk.Label(inputFrame, text='Target linear region (mm)')
        targetRegionLabel.grid(row=2,column=0, pady=5,sticky=tk.W)
        targetRegionInput = ttk.Entry(inputFrame, textvariable=self.linearLength, width = 10, justify=tk.RIGHT)
        targetRegionInput.grid(row=2,column=1,padx = 20)
        
        coilRadiusLabel = ttk.Label(inputFrame, text='Coil radius (mm)')
        coilRadiusLabel.grid(row=3,column=0, pady=5,sticky=tk.W)
        coilRadiusInput = ttk.Entry(inputFrame, textvariable=self.coilRadius, width = 10, justify=tk.RIGHT)
        coilRadiusInput.grid(row=3,column=1,padx = 20)
        
        coilLengthLabel = ttk.Label(inputFrame, text='Coil length (mm)')
        coilLengthLabel.grid(row=4,column=0, pady=5,sticky=tk.W)
        coilLengthInput = ttk.Entry(inputFrame, textvariable=self.coilLength, width = 10, justify=tk.RIGHT)
        coilLengthInput.grid(row=4,column=1,padx = 20)
        
        numWiresLabel = ttk.Label(inputFrame, text='Wires per quadrant')
        numWiresLabel.grid(row=5,column=0, pady=5,sticky=tk.W)
        numWiresInput = ttk.Entry(inputFrame, textvariable=self.numWires, width = 10, justify=tk.RIGHT)
        numWiresInput.grid(row=5,column=1,padx = 20)
        
        wireDiameterLabel = ttk.Label(inputFrame, text='Wire diameter (mm)')
        wireDiameterLabel.grid(row=6,column=0, pady=5,sticky=tk.W)
        wireDiameterInput = ttk.Entry(inputFrame, textvariable=self.wireDiameter, width = 10, justify=tk.RIGHT)
        wireDiameterInput.grid(row=6,column=1,padx = 20)
        
        numHigherTermsLabel = ttk.Label(inputFrame, text='Number of higher order terms')
        numHigherTermsLabel.grid(row=7,column=0, pady=5,sticky=tk.W)
        numHigherTermsInput = ttk.Entry(inputFrame, textvariable=self.higherOrderTerms, width = 10, justify=tk.RIGHT)
        numHigherTermsInput.grid(row=7,column=1,padx = 20)
        
        linearityOrderLabel = ttk.Label(inputFrame, text='Linearity order')
        linearityOrderLabel.grid(row=8,column=0, pady=5,sticky=tk.W)
        linearityOrderInput = ttk.Entry(inputFrame, textvariable=self.linearityOrder, width = 10, justify=tk.RIGHT)
        linearityOrderInput.grid(row=8,column=1,padx = 20)
        
        apodisationTermLabel = ttk.Label(inputFrame, text='Apodisation term')
        apodisationTermLabel.grid(row=9,column=0, pady=5,sticky=tk.W)
        apodisationTermInput = ttk.Entry(inputFrame, textvariable=self.apodisationTerm, width = 10, justify=tk.RIGHT)
        apodisationTermInput.grid(row=9,column=1,padx=20)

        DSVLabel = ttk.Label(inputFrame, text='Simulation DSV (mm)')
        DSVLabel.grid(row=10,column=0, pady=5,sticky=tk.W)
        DSVInput = ttk.Entry(inputFrame, textvariable=self.DSV, width = 10, justify=tk.RIGHT)
        DSVInput.grid(row=10,column=1,padx = 20)

        simResolutionLabel = ttk.Label(inputFrame, text='Simulation resolution (mm)')
        simResolutionLabel.grid(row=11,column=0, pady=5,sticky=tk.W)
        simResolutionInput = ttk.Entry(inputFrame, textvariable=self.simulationRes, width = 10, justify=tk.RIGHT)
        simResolutionInput.grid(row=11,column=1,padx = 20)

        calculateWireButton = ttk.Button(inputFrame, text="Calculate wire pattern", command=self.calculateGradientWires)
        calculateWireButton.grid(row=12, column=0, pady=5, sticky=tk.W)

        self.calculateFieldButton = ttk.Button(inputFrame, text="Simulate B-field", command=self.calculateGradientField)
        self.calculateFieldButton.grid(row=12, column=1, padx=20)
        
        ''' OUTPUT FRAME'''
        outputFrame = ttk.Frame(self)
        outputFrame.grid(row=1,column=0,columnspan=1, sticky=tk.NW, padx=10, pady=20)
        
        outputDescriptionLabel = ttk.Label(outputFrame, text='Simulation output', font='Helvetica 12 bold')
        outputDescriptionLabel.grid(row=0,column=0, columnspan = 1, sticky=tk.W)

        self.gradEfficiencyString   = tk.StringVar(outputFrame, value = "-")
        self.gradErrorString        = tk.StringVar(outputFrame, value = "-")
        self.resistanceString       = tk.StringVar(outputFrame, value = "-")
        self.inductanceString       = tk.StringVar(outputFrame, value = "-")
        self.wireLengthString       = tk.StringVar(outputFrame, value = "-")
        self.exportConjoined        = tk.BooleanVar(outputFrame, value = True)
        self.zRangeString           = tk.StringVar(outputFrame, value = "-")
        
        efficiencyTextLabel = ttk.Label(outputFrame, text='Gradient efficiency:')
        efficiencyTextLabel.grid(row=1,column=0, pady=5,sticky=tk.W)
        efficiencyValueLabel = ttk.Label(outputFrame, textvariable=self.gradEfficiencyString, justify = tk.RIGHT)
        efficiencyValueLabel.grid(row=1,column=1, padx=10,sticky=tk.E)
        
        linearityTextLabel = ttk.Label(outputFrame, text='Error over %.0f mm DSV:'%(self.DSV.get()))
        linearityTextLabel.grid(row=2,column=0, pady=5,sticky=tk.W)
        linearityValueLabel = ttk.Label(outputFrame, textvariable=self.gradErrorString, justify = tk.RIGHT)
        linearityValueLabel.grid(row=2,column=1, padx=10,sticky=tk.E)

        wireLengthTextLabel = ttk.Label(outputFrame, text='Wire length:')
        wireLengthTextLabel.grid(row=3,column=0, pady=5,sticky=tk.W)
        wireLengthValueLabel = ttk.Label(outputFrame, textvariable=self.wireLengthString, justify = tk.RIGHT)
        wireLengthValueLabel.grid(row=3,column=1, padx=10,sticky=tk.E)

        resistanceTextLabel = ttk.Label(outputFrame, text='Coil resistance:')
        resistanceTextLabel.grid(row=4,column=0, pady=5,sticky=tk.W)
        resistanceValueLabel = ttk.Label(outputFrame, textvariable=self.resistanceString, justify = tk.RIGHT)
        resistanceValueLabel.grid(row=4,column=1, padx=10,sticky=tk.E)

        zRangeStringTextLabel = ttk.Label(outputFrame, text='Min/Max X: ')
        zRangeStringTextLabel.grid(row=5,column=0, pady=5,sticky=tk.W)
        zRangeStringValueLabel = ttk.Label(outputFrame, textvariable=self.zRangeString, justify = tk.RIGHT)
        zRangeStringValueLabel.grid(row=5,column=1, padx=10,sticky=tk.E)
        
        self.saveBfieldButton = ttk.Button(outputFrame, text="Export B-field",state=tk.DISABLED, command=self.exportBfield)
        self.saveBfieldButton.grid(row=6, column=0, pady=5, sticky=tk.SW)
        
        saveWireButton = ttk.Button(outputFrame, text="Export wire to CSV", command=self.exportWireCSV)
        saveWireButton.grid(row=6, column=1, pady = 5, padx=20,sticky=tk.SE)
        
        conjoinedExport = ttk.Checkbutton(outputFrame, text="Join wires", variable = self.exportConjoined)
        conjoinedExport.grid(row=7, column=1, pady = 5, padx=20,sticky=tk.SE)
        
        '''Position plots'''
        
        self.contourfig = Figure(figsize=(5,5))
        self.contourfig.gca(projection='3d')
        self.contourfig.set_tight_layout(True)
        
        self.contourCanvas = FigureCanvasTkAgg(self.contourfig, master = self)
        self.contourCanvas.draw()
        self.contourCanvas.get_tk_widget().grid(row=0,column=1)
        self.contourCanvas.mpl_connect('button_press_event', self.contourfig.gca()._button_press)
        self.contourCanvas.mpl_connect('button_release_event', self.contourfig.gca()._button_release)
        self.contourCanvas.mpl_connect('motion_notify_event', self.contourfig.gca()._on_move)
        
        self.fieldFig = Figure(figsize=(5,5))
        self.fieldFig.gca(projection='3d')
        self.fieldFig.set_tight_layout(True)
        
        self.fieldCanvas = FigureCanvasTkAgg(self.fieldFig, master = self)
        self.fieldCanvas.draw()
        self.fieldCanvas.get_tk_widget().grid(row=1,column=1)
        self.fieldCanvas.mpl_connect('button_press_event', self.fieldFig.gca()._button_press)
        self.fieldCanvas.mpl_connect('button_release_event', self.fieldFig.gca()._button_release)
        self.fieldCanvas.mpl_connect('motion_notify_event', self.fieldFig.gca()._on_move)
        self.calculateGradientWires()
        
    def calculateGradientWires(self):
        
        if (self.directionInput.current() == 0):
            self.contourData = gradCalc.halbachXgradient(linearLength = self.linearLength.get(), coilRad = self.coilRadius.get(), coilLength = self.coilLength.get(),\
                                      numWires = self.numWires.get(), numHighOrders = self.higherOrderTerms.get(), linearityTerm = self.linearityOrder.get(),\
                                      apoTerm = self.apodisationTerm.get())
            
        elif(self.directionInput.current() == 1):
            self.contourData =  gradCalc.halbachYgradient(linearLength = self.linearLength.get(), coilRad = self.coilRadius.get(), coilLength = self.coilLength.get(),\
                                      numWires = self.numWires.get(), numHighOrders = self.higherOrderTerms.get(), linearityTerm = self.linearityOrder.get(),\
                                      apoTerm = self.apodisationTerm.get())
                
        elif(self.directionInput.current() == 2):
            self.contourData = gradCalc.halbachZgradient(linearLength = self.linearLength.get(), coilRad = self.coilRadius.get(), coilLength = self.coilLength.get(),\
                                      numWires = self.numWires.get(), numHighOrders = self.higherOrderTerms.get(), linearityTerm = self.linearityOrder.get(),\
                                      apoTerm = self.apodisationTerm.get())
        wireLevels = self.contourData.allsegs

        self.gradEfficiencyString.set("%.4f mT/m/A"%(np.divide(1,self.contourData.levels[1] - self.contourData.levels[0])))
        
        self.contourfig.gca().clear()
        self.contourfig.gca().set_xlabel('Z (mm)')
        self.contourfig.gca().set_ylabel('Y (mm)')
        self.contourfig.gca().set_zlabel('X (mm)')
        np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)                 
        self.length = 0
        for idx, wireLevel in enumerate(wireLevels):
            currentLevel = self.contourData.levels[idx]
            for wire in range (np.size(wireLevel,0)):
                wirePath3D = np.stack((np.cos(wireLevel[wire][:,0])*self.coilRadius.get(), np.sin(wireLevel[wire][:,0])*self.coilRadius.get(),wireLevel[wire][:,1]*1e3),axis=1).astype(float)
                if currentLevel>0:
                    self.contourfig.gca().plot3D(wirePath3D[:,0],wirePath3D[:,1],wirePath3D[:,2],'red')
                else:
                    self.contourfig.gca().plot3D(wirePath3D[:,0],wirePath3D[:,1],wirePath3D[:,2],'blue')

                temp, indices = np.unique(wirePath3D.round(decimals=6), axis = 0, return_index = True)
                wirePath3D = wirePath3D[sorted(indices)]
                delta = wirePath3D[1:,:] - wirePath3D[:-1,:]
                segLength = np.sqrt(np.sum(np.square(delta), axis = 1))
                self.length += np.sum(segLength)# + np.sum(np.square(wirePath3D[0,:] - wirePath3D[-1,:]))

        # self.contourfig.gca().mouse_init()
        self.wireLengthString.set("%.2f meters"%(self.length*1e-3))
        self.contourCanvas.draw()
        self.contourCanvas.flush_events()
        self.wireDict = gradCalc.exportWires(self.contourData,  self.coilRadius.get(), self.directionInput.current(), self.exportConjoined.get())
        minZ = 0
        maxZ = 0
        
        for segmentKey in self.wireDict:
            try:
                segment = self.wireDict[segmentKey]
                minZSegment = np.min(segment[:,2])
                maxZSegment = np.max(segment[:,2])
                if minZSegment < minZ:
                    minZ = minZSegment
                if maxZSegment > maxZ:
                    maxZ = maxZSegment
            except:
                print("Error")
                            
        self.zRangeString.set("%.2f / %.2f mm"%(minZ, maxZ))
        self.calculateResistance()

    def calculateResistance(self):
        copperResistivity = 1.68*1e-8
        resistance = copperResistivity*self.length*1e-3/(np.pi*(self.wireDiameter.get()*1e-3 /2)**2)
        self.resistanceString.set("%.4f Ohms"%(resistance,))
        
    def calculateGradientField(self):
        self.fieldFig.gca().clear()
        self.fieldFig.gca().set_xlabel('Z (mm)')
        self.fieldFig.gca().set_ylabel('Y (mm)')
        self.fieldFig.gca().set_zlabel('X (mm)')
        self.coords, self.bField, error = gradCalc.calculateBfield(self.contourData, self.DSV.get()*1e-3, self.simulationRes.get()*1e-3, self.coilRadius.get()*1e-3, self.directionInput.current())
        bNorm = (self.bField  - np.min(self.bField))/(np.max(self.bField) - np.min(self.bField))
        self.fieldFig.gca().plot_surface(self.coords[0], self.coords[1], self.coords[2],  rstride=1, cstride=1, facecolors=cm.jet(bNorm))
        self.fieldFig.gca().mouse_init()
        self.fieldCanvas.draw()
        self.fieldCanvas.flush_events()
        self.gradErrorString.set("%.2f %%"%(error*100))
        self.saveBfieldButton['state'] = 'normal'
    
    def exportBfield(self):
        bFieldOutputFile = tk.filedialog.asksaveasfilename(defaultextension = '.csv', filetypes=(("CSV file","*.csv"), ("Text file","*.txt")))
        if bFieldOutputFile == None:
            return
        header = 'X (mm),\tY (mm),\tZ (mm),\tBz (mT)'
        delimiter  = ',\t'
        outputArray = np.zeros((np.size(self.bField),np.size(self.coords,0)+1))
        for idx in range(np.size(self.coords,0)):
            outputArray[:,idx] = np.ravel(self.coords[idx])
        outputArray[:,-1] = np.ravel(self.bField)
        np.savetxt(bFieldOutputFile.name , outputArray,delimiter  = delimiter ,header = header )
        bFieldOutputFile.close()
    
    def exportWireCSV(self):
        wireOutputFile = tk.filedialog.asksaveasfilename(defaultextension = '.csv', filetypes=(("CSV file","*.csv"), ("Text file","*.txt")))
        if wireOutputFile == '':
            return
        folder, filename = os.path.split(wireOutputFile)
        file, extension = os.path.splitext(filename)
        
        self.wireDict = gradCalc.exportWires(self.contourData,  self.coilRadius.get(), self.directionInput.current(), self.exportConjoined.get())
        
        for key in self.wireDict:
            filename = os.path.join(folder,file+key+extension)
            np.savetxt(filename,self.wireDict[key],delimiter=",", fmt='%f')
        
        
app = GDT_GUI()
app.mainloop()
closePlots('all')