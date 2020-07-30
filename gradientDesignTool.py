# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 12:57:51 2019

@author: to_reilly
"""

import numpy as np
import tkinter as tk
from tkinter import ttk, font
from tkinter.ttk import Style
import os
import gradientCalculation as gradCalc


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

    def __init__(self, *args, **kwargs):
        super(GDT_GUI, self).__init__(*args, **kwargs)
        # set name and icon
        # super().iconbitmap(self, default=None)
        #super().configure(background='gray')
        super().wm_title('Gradient design tool')

              

        '''INPUT FRAME'''
        #s = Style()
        #s.configure('My.TFrame', background='gray')
        inputFrame = ttk.Frame(self, borderwidth=1)#,style='My.TFrame')
        inputFrame.grid(row=0, column=0, columnspan=1, sticky=tk.W,
                        padx=10, pady=10)
        
        default_font = font.nametofont("TkDefaultFont")
        default_font.configure(size=12) 

        inputDescriptionLabel = ttk.Label(inputFrame, text='Design parameters',
                                          font='Helvetica 18 bold')
        inputDescriptionLabel.grid(row=0, column=0, columnspan=2, sticky=tk.W)
        

        # create variables for storing inputs
        self.linearLength = tk.DoubleVar(inputFrame, value=140)
        self.coilRadius = tk.DoubleVar(inputFrame, value=135)
        self.coilLength = tk.DoubleVar(inputFrame, value=400)
        self.numWires = tk.IntVar(inputFrame, value=15)
        self.wireDiameter = tk.DoubleVar(inputFrame, value=1.5)
        self.apodisationTerm = tk.DoubleVar(inputFrame, value=0.05)
        self.linearityOrder = tk.IntVar(inputFrame, value=16)
        self.higherOrderTerms = tk.IntVar(inputFrame, value=0)
        self.simulationSize = tk.DoubleVar(inputFrame, value=140)
        self.simulationRes = tk.DoubleVar(inputFrame, value=5)
        self.DSV = tk.DoubleVar(inputFrame, value=140)

        directionLabel = ttk.Label(inputFrame, text='Background field direction')
        directionLabel.grid(row=1, column=0, pady=5, sticky=tk.W)
        self.backgroundInput = ttk.Combobox(inputFrame, values=["Halbach", "Conventional"],
                                           width=12, justify=tk.RIGHT,
                                           state='readonly')
        self.backgroundInput.current(0)
        self.backgroundInput.grid(row=1, column=1, padx=20)

        directionLabel = ttk.Label(inputFrame, text='Gradient direction')
        directionLabel.grid(row=2, column=0, pady=5, sticky=tk.W)
        self.directionInput = ttk.Combobox(inputFrame, values=["X", "Y", "Z"],
                                           width=7, justify=tk.RIGHT,
                                           state='readonly')
        self.directionInput.current(0)
        self.directionInput.grid(row=2, column=1, padx=20)

        targetRegionLabel = ttk.Label(inputFrame, text='Target linear region (mm)')
        targetRegionLabel.grid(row=3, column=0, pady=5, sticky=tk.W)
        targetRegionInput = ttk.Entry(inputFrame, textvariable=self.linearLength, width=10, justify=tk.RIGHT)
        targetRegionInput.grid(row=3, column=1, padx=20)

        coilRadiusLabel = ttk.Label(inputFrame, text='Coil radius (mm)')
        coilRadiusLabel.grid(row=4, column=0, pady=5, sticky=tk.W)
        coilRadiusInput = ttk.Entry(inputFrame, textvariable=self.coilRadius, width=10, justify=tk.RIGHT)
        coilRadiusInput.grid(row=4, column=1, padx=20)

        coilLengthLabel = ttk.Label(inputFrame, text='Coil length (mm)')
        coilLengthLabel.grid(row=5, column=0, pady=5, sticky=tk.W)
        coilLengthInput = ttk.Entry(inputFrame, textvariable=self.coilLength, width=10, justify=tk.RIGHT)
        coilLengthInput.grid(row=5, column=1, padx=20)

        numWiresLabel = ttk.Label(inputFrame, text='Wires per quadrant')
        numWiresLabel.grid(row=6, column=0, pady=5, sticky=tk.W)
        numWiresInput = ttk.Entry(inputFrame, textvariable=self.numWires, width=10, justify=tk.RIGHT)
        numWiresInput.grid(row=6, column=1, padx=20)


        # #Todo: remove higher order terms for longidutinal background fields 
        # numHigherTermsLabel = ttk.Label(inputFrame, text='Number of higher order terms')
        # numHigherTermsLabel.grid(row=8, column=0, pady=5, sticky=tk.W)
        # numHigherTermsInput = ttk.Entry(inputFrame, textvariable=self.higherOrderTerms, width=10, justify=tk.RIGHT)
        # numHigherTermsInput.grid(row=8, column=1, padx=20)

        resolutionLabel = ttk.Label(inputFrame, text='Linearity order')
        resolutionLabel.grid(row=7, column=0, pady=5, sticky=tk.W)
        resolutionInput = ttk.Entry(inputFrame, textvariable=self.linearityOrder, width=10, justify=tk.RIGHT)
        resolutionInput.grid(row=7, column=1, padx=20)

        apodisationTermLabel = ttk.Label(inputFrame, text='Apodisation term')
        apodisationTermLabel.grid(row=8, column=0, pady=5, sticky=tk.W)
        apodisationTermInput = ttk.Entry(inputFrame, textvariable=self.apodisationTerm, width=10, justify=tk.RIGHT)
        apodisationTermInput.grid(row=8, column=1, padx=20)
        
       
        '''SIMULATION FRAME'''

        inputDescriptionLabel = ttk.Label(inputFrame, text='Simulation parameters',
                                          font='Helvetica 18 bold')
        inputDescriptionLabel.grid(row=9, column=0, columnspan=2, sticky=tk.W)
       
              
        wireDiameterLabel = ttk.Label(inputFrame, text='Wire diameter (mm)')
        wireDiameterLabel.grid(row=10, column=0, pady=5, sticky=tk.W)
        wireDiameterInput = ttk.Entry(inputFrame, textvariable=self.wireDiameter, width=10, justify=tk.RIGHT)
        wireDiameterInput.grid(row=10, column=1, padx=20)

        linearityOrderLabel = ttk.Label(inputFrame, text='Simulation volume (mm)')
        linearityOrderLabel.grid(row=11, column=0, pady=5, sticky=tk.W)
        linearityOrderInput = ttk.Entry(inputFrame, textvariable=self.simulationSize, width=10, justify=tk.RIGHT)
        linearityOrderInput.grid(row=11, column=1, padx=20)

        simDimensionsLabel = ttk.Label(inputFrame, text='Simulation resolution (mm)')
        simDimensionsLabel.grid(row=12, column=0, pady=5, sticky=tk.W)
        simDimensionsInput = ttk.Entry(inputFrame, textvariable=self.simulationRes, width=10, justify=tk.RIGHT)
        simDimensionsInput.grid(row=12, column=1, padx=20)

        DSVLabel = ttk.Label(inputFrame, text='Error DSV (mm)')
        DSVLabel.grid(row=13, column=0, pady=5, sticky=tk.W)
        DSVInput = ttk.Entry(inputFrame, textvariable=self.DSV, width=10, justify=tk.RIGHT)
        DSVInput.grid(row=13, column=1, padx=20)

        calculateWireButton = ttk.Button(inputFrame, text="Calculate wire pattern", command=self.calculateGradientWires)
        calculateWireButton.grid(row=14, column=0, pady=5, sticky=tk.W)

        self.calculateFieldButton = ttk.Button(inputFrame, text="Simulate B-field", command=self.calculateGradientField)
        self.calculateFieldButton.grid(row=14, column=1, padx=20)
        
        ''' OUTPUT FRAME'''
        outputFrame = ttk.Frame(self)
        outputFrame.grid(row=1, column=0, columnspan=1, sticky=tk.NW, padx=10, pady=10)

        outputDescriptionLabel = ttk.Label(outputFrame, text='Simulation output', font='Helvetica 18 bold')
        outputDescriptionLabel.grid(row=0, column=0, columnspan=1, sticky=tk.W)

        self.gradEfficiencyString   = tk.StringVar(outputFrame, value="-")
        self.gradErrorString        = tk.StringVar(outputFrame, value="-")
        self.resistanceString       = tk.StringVar(outputFrame, value="-")
        self.inductanceString       = tk.StringVar(outputFrame, value="-")
        self.wireLengthString       = tk.StringVar(outputFrame, value="-")
        self.DSVerrorString         = tk.StringVar(outputFrame, value="Error over {:.0f} mm DSV:".format(self.DSV.get()))
        self.exportConjoined        = tk.BooleanVar(outputFrame, value=True)

        efficiencyTextLabel = ttk.Label(outputFrame, text='Gradient efficiency:')
        efficiencyTextLabel.grid(row=1, column=0, pady=5, sticky=tk.W)
        efficiencyValueLabel = ttk.Label(outputFrame, textvariable=self.gradEfficiencyString, justify=tk.RIGHT)
        efficiencyValueLabel.grid(row=1, column=1, padx=10, sticky=tk.E)

        # TODO: This doesn't update for a change in DSV value!
        linearityTextLabel = ttk.Label(outputFrame, textvariable=self.DSVerrorString)
        linearityTextLabel.grid(row=2, column=0, pady=5, sticky=tk.W)
        linearityValueLabel = ttk.Label(outputFrame, textvariable=self.gradErrorString, justify=tk.RIGHT)
        linearityValueLabel.grid(row=2, column=1, padx=10, sticky=tk.E)

        wireLengthTextLabel = ttk.Label(outputFrame, text='Wire length:')
        wireLengthTextLabel.grid(row=3, column=0, pady=5, sticky=tk.W)
        wireLengthValueLabel = ttk.Label(outputFrame, textvariable=self.wireLengthString, justify=tk.RIGHT)
        wireLengthValueLabel.grid(row=3, column=1, padx=10, sticky=tk.E)

        resistanceTextLabel = ttk.Label(outputFrame, text='Coil resistance:')
        resistanceTextLabel.grid(row=4, column=0, pady=5, sticky=tk.W)
        resistanceValueLabel = ttk.Label(outputFrame, textvariable=self.resistanceString, justify=tk.RIGHT)
        resistanceValueLabel.grid(row=4, column=1, padx=10, sticky=tk.E)

        #inductanceTextLabel = ttk.Label(outputFrame, text='Coil inductance:')
        #inductanceTextLabel.grid(row=5, column=0, pady=5, sticky=tk.W)
        #inductanceValueLabel = ttk.Label(outputFrame, textvariable=self.inductanceString, justify=tk.RIGHT)
        #inductanceValueLabel.grid(row=5, column=1, padx=10, sticky=tk.E)

        self.saveBfieldButton = ttk.Button(outputFrame, text="Export B-field", state=tk.DISABLED, command=self.exportBfield)
        self.saveBfieldButton.grid(row=6, column=0, pady=5, sticky=tk.SW)

        saveWireButton = ttk.Button(outputFrame, text="Export wire to CSV", command=self.exportWireCSV)
        saveWireButton.grid(row=6, column=1, pady=5, padx=20, sticky=tk.SE)

        conjoinedExport = ttk.Checkbutton(outputFrame, text="Join wires", variable=self.exportConjoined)
        conjoinedExport.grid(row=7, column=1, pady=5, padx=20, sticky=tk.SE)

        '''Position plots'''

        self.contourfig = Figure(figsize=(7, 7))
        self.contourfig.gca(projection='3d')
        self.contourfig.set_tight_layout(True)

        self.contourCanvas = FigureCanvasTkAgg(self.contourfig, master=self)
        self.contourCanvas.draw()
        self.contourCanvas.get_tk_widget().grid(row=0, column=1)

        self.fieldFig = Figure(figsize=(7, 7))
        self.fieldFig.gca(projection='3d')
        self.fieldFig.set_tight_layout(True)

        self.fieldCanvas = FigureCanvasTkAgg(self.fieldFig, master=self)
        self.fieldCanvas.draw()
        self.fieldCanvas.get_tk_widget().grid(row=0, column=2)
       

        self.calculateGradientWires()

    def calculateGradientWires(self):

        designParameters = {}
        designParameters['length'] = self.linearLength.get()*1e-3
        designParameters['linearity'] = self.linearityOrder.get()

        designParameters['resolution'] = self.simulationRes.get()*1e-3
       # designParameters['gridSize'] = self.simulationSize.get()*1e-3
        designParameters['gridSize'] = self.coilLength.get()*1e-3

        designParameters['B0'] = self.backgroundInput.get()
        designParameters['coilRad'] = self.coilRadius.get()*1e-3
        # TODO: arbitrary scaling of b, why this value?? comment and win!
        designParameters['fieldRad'] = 0.001*self.coilRadius.get()*1e-3
        designParameters['apo'] = self.apodisationTerm.get()
        designParameters['nrModes'] = self.higherOrderTerms.get()
        designParameters['nrWires'] = self.numWires.get()
       
     
        # Halbach system:
        if (self.backgroundInput.current() == 0):
         
            if (self.directionInput.current() == 0):
                designParameters['type'] = "longitudinal"
                designParameters['gradDir'] = "z"

            elif(self.directionInput.current() == 2):
                designParameters['type'] = "transverse"
                designParameters['gradDir'] = "x"

            elif(self.directionInput.current() == 1):
                designParameters['type'] = "transverse"
                designParameters['gradDir'] = "y"
            
        # Conventional MRI:        
        elif (self.backgroundInput.current() == 1):
            
            if (self.directionInput.current() == 0):
                designParameters['type'] = "transverse"
                designParameters['gradDir'] = "x"
                
            elif(self.directionInput.current() == 1):
                designParameters['type'] = "transverse"
                designParameters['gradDir'] = "y"

            elif(self.directionInput.current() == 2):
                designParameters['type'] = "longitudinal"
                designParameters['gradDir'] = "z"
                


        Coil = gradCalc.generateCoil(designParameters, gradStrength=1e-3, current=2)

        self.contourData = Coil['wires']
        wireLevels = self.contourData.allsegs

        self.inductanceString.set("{:.4g} \u03BCH".format(Coil['L']*1e6))
        

        self.contourfig.gca().clear()
        
        if (self.backgroundInput.current() == 0):
            self.contourfig.gca().set_xlabel('Z/B$_{0}$ (mm)')
            self.contourfig.gca().set_ylabel('X (mm)')
            self.contourfig.gca().set_zlabel('Y (mm)')
        else:
            self.contourfig.gca().set_xlabel('X (mm)')
            self.contourfig.gca().set_ylabel('Z/B$_{0}$ (mm)')
            self.contourfig.gca().set_zlabel('Y (mm)')
            

        for idx, wireLevel in enumerate(wireLevels):
            currentLevel = self.contourData.levels[idx]
            for wire in range(np.size(wireLevel, 0)):
                wirePath3D = np.stack((np.cos(wireLevel[wire][:, 0])*self.coilRadius.get(),
                                       np.sin(wireLevel[wire][:, 0])*self.coilRadius.get(),
                                       wireLevel[wire][:, 1]*1e3), axis=1)
                if currentLevel > 0:
                    self.contourfig.gca().plot3D(wirePath3D[:, 0],
                                                 wirePath3D[:, 2],
                                                 wirePath3D[:, 1], 'red')
                else:
                    self.contourfig.gca().plot3D(wirePath3D[:, 0],
                                                 wirePath3D[:, 2],
                                                 wirePath3D[:, 1], 'blue')
                
        self.length = Coil['length']
        self.wireLengthString.set("{:.2f} metres".format(self.length))
        self.contourfig.gca().mouse_init()
        self.contourCanvas.draw()
        self.contourCanvas.flush_events()
        self.DSVerrorString.set('Error over {:.0f} mm DSV:'.format(self.DSV.get()))
        self.resistanceString.set("{:.4f} \u03A9".format(Coil['R']))

    def calculateGradientField(self):
        self.calculateFieldButton['state'] = 'disabled'
        super().update()

        if (self.backgroundInput.current() == 0):        
            self.fieldFig.gca().set_xlabel('Z/B$_{0}$ (mm)')
            self.fieldFig.gca().set_ylabel('X (mm)')
            self.fieldFig.gca().set_zlabel('Y (mm)')
        else:
            self.fieldFig.gca().set_xlabel('X (mm)')
            self.fieldFig.gca().set_ylabel('Z/B$_{0}$ (mm)')
            self.fieldFig.gca().set_zlabel('Y (mm)')

        self.coords, self.bField, self.efficiency, error= gradCalc.calculateBfield(self.contourData, self.DSV.get()*1e-3, self.simulationRes.get()*1e-3, self.coilRadius.get()*1e-3, self.directionInput.current(), self.backgroundInput.get()) 
        
        
        bNorm = (self.bField  - np.min(self.bField))/(np.max(self.bField) - np.min(self.bField))
        self.fieldFig.gca().plot_surface(self.coords[0], self.coords[2], self.coords[1],  rstride=1, cstride=1, facecolors=cm.jet(bNorm))
        
        self.fieldCanvas.draw()
        self.fieldCanvas.flush_events()
        self.gradEfficiencyString.set("{:.4f} mT/m/A".format(self.efficiency*1000))    
        self.gradErrorString.set("{:.2f} %".format(error*100))
        self.fieldFig.gca().mouse_init()

        self.calculateFieldButton['state'] = 'normal'
        self.saveBfieldButton['state'] = 'normal'
        super().update()

    def exportBfield(self):
        bFieldOutputFile = tk.filedialog.asksaveasfilename(mode='w', defaultextension='.csv', filetypes=(("CSV file", "*.csv"), ("Text file", "*.txt")))
        if bFieldOutputFile is None:
            return
        header = 'X (mm),\tY (mm),\tZ (mm),\tBz (mT)'
        delimiter = ',\t'
        outputArray = np.zeros((np.size(self.bField), np.size(self.coords, 0)+1))
        for idx in range(np.size(self.coords, 0)):
            outputArray[:, idx] = np.ravel(self.coords[idx])
        outputArray[:, -1] = np.ravel(self.bField)
        np.savetxt(bFieldOutputFile.name, outputArray, delimiter=delimiter, header=header)
        bFieldOutputFile.close()

    def exportWireCSV(self):
        wireOutputFile = tk.filedialog.asksaveasfilename(defaultextension='.csv', filetypes=(("CSV file","*.csv"), ("Text file", "*.txt")))
        if wireOutputFile == '':
            return
        folder, filename = os.path.split(wireOutputFile)
        file, extension = os.path.splitext(filename)

        self.wireDict = gradCalc.exportWires(self.contourData,
                                             self.coilRadius.get(),
                                             self.directionInput.current(),
                                             self.exportConjoined.get())

        for key in self.wireDict:
            filename = os.path.join(folder, file+key+extension)
            np.savetxt(filename, self.wireDict[key], delimiter=",", fmt='%f')


app = GDT_GUI()
app.mainloop()
closePlots('all')
