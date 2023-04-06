# -*- coding: utf-8 -*-
"""
Created on Mon Mar 06 10:50:53 2017

@author: JakobiA
"""
# import sys

# from PyQt4 import QtCore,QtGui
# # die folgende .py wurde aus BaQtDesigner_V3.ui erstellt
# # durch pyuic4 BaQtDesigner_V3.ui -o BaQtDesigner_V3.py
# # in der WIN Konsole
# # das widget layout kann dort entsprechend angepasst werden
# from BaQtDesigner_V3 import Ui_MainWindow
import numpy as np
# from PIL import Image, ImageDraw
# from matplotlib import pyplot as plt
# import os.path
# import time
# import CalculateClass_ForQA as CALCULATIONS


def NormalizeCenter100(arr, arr_x, arr_y):
    """ norm data to their center, center of array equals 100
    
    Parameters
    ---------------------
    arr: array of number values
    arr_x: x-Axis
    arr_y: y-Axis

    Returns
    ------------------
    normalized array
    """

    print("NormalizeCenter100")
    
    # Parameter für Normierung auf Zentrumsbereich (Groeße)
    # sind immer so, nicht variabel
    X_NORM_CENTRUM = 8    # 8mm x direction
    Y_NORM_CENTRUM = 8    # 8mm y direction
    
    x_axarr, y_axarr = np.meshgrid(arr_x, arr_y)
    index = (x_axarr >= -X_NORM_CENTRUM) * (x_axarr <= X_NORM_CENTRUM) *\
            (y_axarr >= -Y_NORM_CENTRUM) * (y_axarr <= Y_NORM_CENTRUM)
    
    normval = np.around(np.mean(arr[index]), decimals=1)
    normalizedArr = 100.0 * arr / normval      

    return normalizedArr, normval

#------------------------------------------------------------------------------    
def NormalizeMax100(arr):
    """ normalize the array to the maximum (max=100)
    
      Parameters
      -------------        
      arr: array of number values

      Returns
      ----------
      normalized array
     """

    print("NormalizeMax100")
    normalizedArr=100*arr/arr.max()
    
    return normalizedArr

#------------------------------------------------------------------------------        
def NormalizeSOBP(Doses, Depths):
    """ normalize the array to the average value of the middle third of the SOBP
    
      Parameters
      -------------        
      Doses: array of Dose values
      Depths: array of belonging Depths

      Returns
      ----------
      normalized dose array
      normvalue
    """
    print("NormSOBP")
    # finde SOBP Bereich anhand R95 dist und prox
    dis95 = CALCULATIONS.CalculatePos('dep',"d",95,Doses,Depths)
    prox = CALCULATIONS.CalculatePos('dep',"p",95,Doses,Depths)
    
    # Falls prox 95 nicht definiert, da Dosisplateau so hoch, nehme R98
    if (prox == 'nan'):
        print('Normierung mit prox98')            
        prox = CALCULATIONS.CalculatePos('dep',"p",98,Doses,Depths)
    
    # setze Bereich der analysiert werden soll: 1/3 des Plateaus
    dis=1.0/3.0*abs(dis95-prox)  
    proxPoint=prox+dis #  vordere normarea punkt verschoben um 1/3 des plateasu nach hinten
    distPoint=dis95-dis #  hinten normarea punkt verschoben um 1/3 des plateasu nach vorne
    
    ### normalize to average der 1/3 plateau groesse
    x1=np.where(Depths>=proxPoint)[0] # alle Positionen im array die hinter dem gew prox punktes liegen
    x2=np.where(Depths<=distPoint)[0] # alle Positionen im array die vor dem gew dist punktes liegen
    if proxPoint+0.5<distPoint: # wenn es tatsachlich ein Plateau gibt von relevanter Groesse
        normval = Doses[x1[0]:x2[-1]+1].mean() # Mittelwert des Plateaus zur Normierung
    else: # sonst normierung auf Maximum
        print("norm to maximum")
        normval = max(Doses)  

    Doses=Doses*100.0/normval #Normierung der Dosiswwerte, auf Prozent gerechnet
               
    return Doses, normval
#------------------------------------------------------------------------------        
def NormalizeLatSOBP(Doses, coords, iteration):
    
    print("NormLatSOBP")
    
    abscoords = abs(coords)
    indnum = max(np.nonzero(abscoords == min(abscoords))[0])
    
    normval=sum(Doses[indnum-20:indnum+21])/41.      # Normierung
    Doses=Doses*100./normval
       
    return Doses, normval      
