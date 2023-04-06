# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 16:52:43 2017
Taken from BaMainProgram, as extra file, usable in other files
@author: JakobiA
"""

#Bibliotheken
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

#------------------------------------------------------------------------------
def CalculatePos(curve, location, value, doses, coords):                     
    """ 
     Calculates Position:
     interpolate between:
       	the first (proximal) pair of values, which was identified 
                   (dose & according coordinate) and the pair of values 
                   one index in front
       	the last (distal) pair of values, which was identified 
                   (dose & according coordinate) and the pair of values 
                   one index after it
                   calculate the coordinates, where the dose would equals the “value””
                   
     Parameters
     ----------------
     curve: if lateral or depth dose analyses
     location:  location where the data will be examined
       	“p” – proximal
      	“d” – distal
      	“pd” – proximal-distal
        
     value: value, its exact position will be determined
     
     doses: array of the data (has to come under and over the “value”)
     
     coords: array of coordinates belong to “doses”
         identify all indices of the doses, which come over the “value”
         

    Returns
    -----------
      proximal(location=”p”) or distal(location=”d”) coordinate, 
      where doses equals “value”, or the distance between this coordinates(location=”pd”)
        """
  
    if curve == 'lat':
        threshold = 0
    else:
        threshold = max(coords[doses == max(doses)])
    
    # Calculate position for coordinates < 0
    try:
        indnum = max(np.nonzero((coords < threshold) * (doses <= value))[0])
        coord1 = coords[indnum] + \
                 (coords[indnum+1] - coords[indnum]) * \
                 ((value - doses[indnum]) / \
                 (doses[indnum+1] - doses[indnum]))

    except:
        coord1 = 'nan'
#        print "Except coord 1"
    
    # Calculate position for coordinates > 0
    try:
        indnum = min(np.nonzero((coords > threshold) * (doses <= value))[0])
        #print("Indnum hier {0}".format(indnum))
            
        coord2 = coords[indnum] + \
                 (coords[indnum] - coords[indnum-1]) * \
                 ((value - doses[indnum]) / \
                 (doses[indnum] - doses[indnum-1]))

    except:
        coord2 = 'nan'
#        print "Except coord 2"
    
    if location=="p":
        try:
            return coord1
        except:
            return coord1
    elif location=="d":
        try:
            return coord2
        except:
            return coord2
    elif location == "pd":
        try:
            return abs(coord2-coord1)
        except:
            return 'nan'

#------------------------------------------------------------------------------   
def CalculateAverage(arr, y, x, x_axarr, y_axarr, x_off=8, y_off=8, ret=False):  
    '''
    Calculates the average value around a single point with given x- and y-coordinates
    offset is fixed to 8mm (used in original source code), if not changed in function start
    
    Parameters
    -----------------------------
    arr: array(2D) of image data (float)
    y: Position shift from center
    x: Position shift from center
    x_axarr: x data of image array (x_cl as meshgrid, x central position)
    y_axarr:y data of image array (y_il as meshgrid, y central position)
    x_off: defines size of array to be averaged in x dir
    y_off: defines size of array to be averaged in y dir
    
    Returns
    --------------------------
    rounded value of the average
    '''
#    print("Calculate Average - Extern function")
    index = (x_axarr >= x - x_off) *\
            (x_axarr <= x + x_off) *\
            (y_axarr >= y - y_off) *\
            (y_axarr <= y + y_off)        
    val = np.around(np.mean(arr[index]), decimals=1)
    return [val, index]
    
#------------------------------------------------------------------------------        
# Bestimmung der Breite anhand von Wendepunkten
def CalculateWidth(dose, coord, farbe):
    ''' fuer wedge Phantom / daily QA
    Bestimmung der Reichweite anhand der Wendepunkte auf der Diagonalen Messung
    R1/R2
    '''
    # Glaettung
    smind = 4
    dose_smoothed = np.zeros_like(dose)
    
    for i in np.arange(smind, len(dose)-smind):
        dose_smoothed[i] = dose[i-smind:i+(smind+1)].mean()
        
    ableitung = np.gradient(dose_smoothed)
    
    if farbe == 'black' or farbe == 'blue':
        indmin = (coord <= 0) * (coord >= -100)
        indmax = (coord >= 0) * (coord <= 100)
    else:
        indmax = (coord <= 0) * (coord >= -100)
        indmin = (coord >= 0) * (coord <= 100)
    coord_max = coord[indmax]    
    pos_1 = coord_max[ableitung[indmax] == max(ableitung[indmax])][0]
    coord_min = coord[indmin]
    pos_2 = coord_min[ableitung[indmin] == min(ableitung[indmin])][0]
    
    width = abs(pos_1 - pos_2)     
 
    return width, pos_1, pos_2

#------------------------------------------------------------------------------  
def CalculatePenumbra(location,doses,coords):                         
    """ Calcualted Penumbra as D80-D20
    locate the coordinate at the given location, where doses equals 80
    locate the coordinate at the given location, where doses equals 20 

    Parameters
    ----------------
    location: location where the data will be examined
            “p” – proximal
           “d” – distal
    d_norm: array(1D) dose values
    coords: array of coordinates according to the doses

	

   Returns
   ---------------
       distance between the coordinates
    """
#    print("Calculate Penumbra "+ location + " - Extern function")
    
    coord_80 = CalculatePos('lat',location,80,doses,coords)
    coord_20 = CalculatePos('lat',location,20,doses,coords)     
    
    # print coord_80
    # print coord_20
    result=abs(round(float(coord_80)-float(coord_20),1))
    return result
#------------------------------------------------------------------------------
def CalculateLateralOffset(doses,coords):                             
    ''' locate the coordinate, where doses equals 50 first time (negative value)
        locate the coordinate, where doses equals 50 last time (positive value)
        summarize both coordinates and divide the result by 2 to find central axis
 
     Parameters
     -----------------------
     doses: array(1D) of lateral dose values
     coords: array of coordinates according to the doses
     
     Returns
     --------------
     rounded result of the shift from the null-location
     '''

#    print("Calculate Lateral Offset - Extern function")
    
    coord_p = CalculatePos('lat',"p",50,doses,coords)
    coord_d = CalculatePos('lat',"d",50,doses,coords)
    
    shift= round((float(coord_p) + float(coord_d))/2.0,1)
    
    return shift

#-------------------------------------------------------------------------    
def CalculateMinMax(prox,dist,d_norm,coords):                         
    """Calculates Minimum and Maximum
    	indentifies the coordinates, which are bigger than the start position and 
                  smaller than the end position, 	takes the according data, which belong to the coordinates, 
                  finds the minimum and the maximum of this particular range

    Parameters
    ------------------
    prox: start position
    dist: end position
    d_norm: array(1D) of normalized data (e.g. dose values)
    coords: array of coordinates (e.g. x-coordinates)
    
    Returns
    --------------
       value of the minimum of the particular range
       value of the maximum of the particular range
       
     """

#    print("Uniformity Region - Extern function")
    p=np.where(coords>prox)
    d=np.where(coords<dist)
    arr=d_norm[p[0][0]:d[0][-1]]
    minVal=min(arr)
    maxVal=max(arr)
    
    # Calculation Symmetry
    if len(arr) % 2 != 0:
        arr = arr[:-1]
    arr_left = arr[:int(len(arr)/2)]
    arr_right = arr[int(len(arr)/2):]
    arr_right = arr_right[::-1]
    sym = np.amax(abs((arr_left/arr_right*1.0) - 1))
    
    return minVal, maxVal, sym
    

#-------------------------------------------------------------------
def CalculateAndCorrectLateralOffset(x_cl, d_cl, y_il, d_il, arr_norm):
    """
    Calculates lateral offset of a profile at 50 % dose of the profile at central axis of image array
    corrects the image array data for it
    gets new profile doses at center of corrected image array
    check a second and third time for offset, and does the correction again, if required
    
    Parameters
    --------------------------------
     x_cl: array of coordinates for the crossline calculation 
                   (moving in x direction)
                   
       d_cl: array of doses for the crossline calculation 
                   (moving in x direction)
                   
       y_il: array of coordinates for the inline calculation 
                   (moving in y direction)
                   
       d_il: array of doses for the inline calculation 
                   (moving in y direction)
        arr_norm: normed array of doses, eg lynx measurement plane
        
    Returns
    -----------------------------------------------
    corrected x_cl, d_cl, y_il, d_il, array_norm_zentriert
    """
    for i in [0,1,2]:
        # Berechne lateralen Offset
        LateralOffset_cl = CalculateLateralOffset(d_cl, x_cl)
        LateralOffset_il = CalculateLateralOffset(d_il, y_il)   
#        print LateralOffset_cl, "latOff_cl"
#        print LateralOffset_il, "latOff_il"
        
        if str(LateralOffset_cl) == "nan":
            LateralOffset_cl = 0.0
            print("lat offset cl set to 0.0")
        if str(LateralOffset_il) == "nan":
            LateralOffset_il = 0.0
            print("lat offset il set to 0.0")
        ######## Verschiebung der Bildmatrix um lateralen offset: Runden des Offsets auf ganze Pixel, integer setzen
        latOffset_il_pixel = int(round(LateralOffset_il/abs(y_il[0]-y_il[1]),0))
        latOffset_cl_pixel = int(round(LateralOffset_cl/abs(x_cl[0]-x_cl[1]),0))    
        
#        print latOffset_cl_pixel, latOffset_il_pixel, "Pixeel"
        # Verschieben der Bildmatrix array_norm um lateralen offset, Bild dann in Mitte, x_cl, y_il dann für symmetrische Abbildung korrekte Werte
        array_norm_zentriert = arr_norm[:]*0.# anlegen 0 array
        length_il = len(y_il)
        length_cl = len(x_cl)     
        
        if LateralOffset_il >= 0:
            if LateralOffset_cl >=0:
#                print("Offset inline positiv; Offset crossline positiv")                
                array_norm_zentriert[0:length_il-latOffset_il_pixel, 0:length_cl-latOffset_cl_pixel] = arr_norm[latOffset_il_pixel:,latOffset_cl_pixel:]
            else:
#                print("Offset inline positiv; Offset crossline negativ")                
                array_norm_zentriert[0:length_il-latOffset_il_pixel, abs(latOffset_cl_pixel):] = arr_norm[latOffset_il_pixel:,0:length_cl-abs(latOffset_cl_pixel)]
        else:
            if LateralOffset_cl >=0:
#                print("Offset inline negativ; Offset crossline positiv")                
                array_norm_zentriert[abs(latOffset_il_pixel):, 0:length_cl-latOffset_cl_pixel] = arr_norm[0:length_il-abs(latOffset_il_pixel),latOffset_cl_pixel:]
            else:
#                print("Offset inline negativ; Offset crossline negativ")                
                array_norm_zentriert[abs(latOffset_il_pixel):, abs(latOffset_cl_pixel):] = arr_norm[0:length_il-abs(latOffset_il_pixel),0:length_cl-abs(latOffset_cl_pixel)]   
            
        ## berechne in folgender schleife mit neuem array
        arr_norm = array_norm_zentriert[:]*1.# Übergabe der Varuablen
        ## lese d_cl, d_il jetzt in neuem array aus --> wirklich in Bildmitte!
        ## korrigierte Daten von d_cl, d_il
        x_clarr, y_ilarr = np.meshgrid(x_cl, y_il)
        d_cl = arr_norm[y_ilarr == 0]
        d_il = arr_norm[x_clarr == 0] 
    
    return x_cl, d_cl, y_il, d_il, array_norm_zentriert

#------------------------------------------------------------------------------   
def GetSmoothedCurve(smind, profile, d, uniformity_coord_p, uniformity_coord_d, include_fringe=False):
    if include_fringe:
        first_index = np.where(profile <= uniformity_coord_p)[0][-1]
        last_index = np.where(profile >= uniformity_coord_d)[0][0]
    else:       
        first_index = np.where(profile >= uniformity_coord_p)[0][0]
        last_index = np.where(profile <= uniformity_coord_d)[0][-1]
   
    d_smoothed = np.copy(d)
    d_cutoff = np.copy(d)  
    d_cutoff[:first_index] = d[first_index]
    d_cutoff[last_index:] = d[last_index]
    for i in xrange(first_index, last_index):
        d_smoothed[i] = d_cutoff[i-smind:i+(smind+1)].mean()
    return d_smoothed, d_cutoff


#------------------------------------------------------------------------------            
def CalculateProfileData(x_cl, d_cl, y_il, d_il, typeQA, option):               
    """
     Calculate Lateral Profile Values: Field Size, Penumbra, Symm, Flatness, lateral offset, uniformity region, uniformity (Max-Min)/Max;
     used for daily, weekly, yearly QA
     in case of yearly QA: only symmetry and flatness are calculated, field size and penumbra are set to -1!
     
     Parameters
     -------------------
       x_cl: array of coordinates for the crossline calculation 
                   (moving in x direction)
                   
       d_cl: array of doses for the crossline calculation 
                   (moving in x direction)
                   
       y_il: array of coordinates for the inline calculation 
                   (moving in y direction)
                   
       d_il: array of doses for the inline calculation 
                   (moving in y direction)
       
       typeQA:  NEW!
               to deviate between daily, weekly, yearly QA
               ("Ta", "Wo", "Ja")
               
       option: For changes in the analysed size in option 8

     Returns: 
     -------
     array of fieldSize_il, y_low, y_high, flat_il, sym_il, fieldSize_cl, x_low, x_high, flat_cl, sym_cl
    """
    
#    print("Calculate Profile Data - Extern Function")
    if typeQA == "Ja":
        # jaehrliche QA: Messung lat. Profile gantrywinkelabhaengig mit Lynxhalter, d.h. ohne Apertur
        # -> Berechnung von Feldgroesse, Penumbra, Korrektur um lateral Offset nicht sinnvoll
        # -> Feldgroesse und Penumbra werden auf -1 gesetzt, um Funktion CalculateProfileData() dennoch nutzen zu koennen
        fieldSize_il = -1
        fieldSize_cl = -1
        x_low = -1
        x_high = -1
        y_low = -1
        y_high = -1
        latOffset_il = -1
        latOffset_cl = -1
        
        # uniformity region inline/crossline
        # fuer Optionen 1-7 wird klinisch max. verwendete Feldgröße durch max. Aperturdurchmesser fuer Snout180 (=16cm) begrenzt
        # mit SAD 236.8cm ergibt sich daraus eine max. Feldradius im Isozentrum von 10.14cm
        uniformity_il_coord_p = -102.0
        uniformity_il_coord_d =  102.0
        uniformity_cl_coord_p = -102.0
        uniformity_cl_coord_d =  102.0
        if option[0:3] == "DS8":
            uniformity_il_coord_p = -60.0
            uniformity_il_coord_d =  60.0
            uniformity_cl_coord_p = -60.0
            uniformity_cl_coord_d =  60.0
        
    else:
        ## Feldgroesse inline, crossline, in Lynx Zentrum!
        fieldSize_il = CalculatePos('lat',"pd", 50, d_il, y_il)
        fieldSize_cl = CalculatePos('lat',"pd", 50, d_cl, x_cl)
        # Penumbra crossline
        x_low = CalculatePenumbra("p", d_cl, x_cl)
        x_high = CalculatePenumbra("d", d_cl, x_cl)
        # Penumbra inline
        y_low=CalculatePenumbra("p",d_il,y_il)
        y_high=CalculatePenumbra("d",d_il,y_il)
        # Lateraler offset, inline crossline
        latOffset_il=CalculateLateralOffset(d_il, y_il)
        latOffset_cl=CalculateLateralOffset(d_cl, x_cl)
        
        # Korrektur der Profile um laterales Offset! --> fuer uniformity region
        y_il = y_il - latOffset_il
        x_cl = x_cl - latOffset_cl
        # uniformity region inline    
        uniformity_il_coord_p=CalculatePos('lat',"p",50,d_il,y_il)+\
                                                2*(y_low+y_high)/2
        uniformity_il_coord_d=CalculatePos('lat',"d",50,d_il,y_il)-\
                                                2*(y_low+y_high)/2
        # uniformity region crossline
        uniformity_cl_coord_p=CalculatePos('lat',"p",50,d_cl,x_cl)+\
                                               2*(float(x_low)+float(x_high))/2
        uniformity_cl_coord_d=CalculatePos('lat',"d",50,d_cl,x_cl)-\
                                               2*(float(x_low)+float(x_high))/2
                                            
        # Einschraenken der Uniformity Region fuer woechentliche und monatliche QA, DS8
        if typeQA == "Wo":
            if option[0:3] == "DS8":
                uniformity_il_coord_p = -55.0
                uniformity_il_coord_d =  55.0
                uniformity_cl_coord_p = -55.0
                uniformity_cl_coord_d =  55.0
                
        if typeQA == "Mo":
            if option[0:3] == "DS8":
                uniformity_il_coord_p = -66.1
                uniformity_il_coord_d =  66.1
                uniformity_cl_coord_p = -66.1
                uniformity_cl_coord_d =  66.1
                
        if typeQA == "TaDDLP":
            if option[0:3] == "DS8":
                uniformity_il_coord_p = -60.0
                uniformity_il_coord_d =  60.0
                uniformity_cl_coord_p = -60.0
                uniformity_cl_coord_d =  60.0
   
    # Glaettung der Dosis (inline, crossline) zur Bestimmung der Flatness und Symmetrie
    smind = 2 #Glaettungsfaktor
    
    #d_il_smoothed = np.zeros_like(d_il)
    #for i in np.arange(smind, len(d_il)-smind):
    #    d_il_smoothed[i] = d_il[i-smind:i+(smind+1)].mean()               
               
    #d_cl_smoothed = np.zeros_like(d_cl)
    #for i in np.arange(smind, len(d_cl)-smind):
    #    d_cl_smoothed[i] = d_cl[i-smind:i+(smind+1)].mean()        

    # Bessere Glättung, berücksichtigt nur Werte innerhalb des Uniformitätsbereichs
    d_cl_smoothed, _ = GetSmoothedCurve(smind, x_cl, d_cl, uniformity_cl_coord_p, uniformity_cl_coord_d)
    d_il_smoothed, _ = GetSmoothedCurve(smind, y_il, d_il, uniformity_il_coord_p, uniformity_il_coord_d)
    
    ##     Symmentrie crossline inline, sowie Daten fuer Flatness
    minVal_cl, maxVal_cl, sym_cl =CalculateMinMax(uniformity_cl_coord_p,
                                                   uniformity_cl_coord_d,
                                                   d_cl_smoothed,x_cl) 
    minVal_il, maxVal_il, sym_il =CalculateMinMax(uniformity_il_coord_p,
                                                   uniformity_il_coord_d,
                                                   d_il_smoothed,y_il)     
                                                       
    # Flatness inline crossline
    flat_cl=maxVal_cl/minVal_cl - 1.0
    flat_il = maxVal_il/minVal_il - 1.0
    
    return [fieldSize_il, y_low, y_high, latOffset_il, [round(uniformity_il_coord_p,1),round(uniformity_il_coord_d,1)], round(flat_il, 3), round(sym_il,3), fieldSize_cl, x_low, x_high, latOffset_cl, [round(uniformity_cl_coord_p,1),round(uniformity_cl_coord_d,1)], round(flat_cl, 3), round(sym_cl, 3)]





#------------------------------------------------------------------------------            
def CalculateProfileDataSmoothedCurves(x_cl, d_cl, y_il, d_il, typeQA, option):               
    """
     Calculate the smoothed curves in the same way as CalculateProfileData does and return them for plotting
     
     Parameters
     -------------------
       x_cl: array of coordinates for the crossline calculation 
                   (moving in x direction)
                   
       d_cl: array of doses for the crossline calculation 
                   (moving in x direction)
                   
       y_il: array of coordinates for the inline calculation 
                   (moving in y direction)
                   
       d_il: array of doses for the inline calculation 
                   (moving in y direction)
       
       typeQA:  NEW!
               to deviate between daily, weekly, yearly QA
               ("Ta", "Wo", "Ja")
               
       option: For changes in the analysed size in option 8

     Returns: 
     -------
     The smoothed crossline curve, smoothed inline curve, crossline data used for smoothing, inline data used for smoothing
    """
    
#    print("Calculate Profile Data - Extern Function")       
    ## Feldgroesse inline, crossline, in Lynx Zentrum!
    fieldSize_il = CalculatePos('lat',"pd", 50, d_il, y_il)
    fieldSize_cl = CalculatePos('lat',"pd", 50, d_cl, x_cl)
    # Penumbra crossline
    x_low = CalculatePenumbra("p", d_cl, x_cl)
    x_high = CalculatePenumbra("d", d_cl, x_cl)
    # Penumbra inline
    y_low=CalculatePenumbra("p",d_il,y_il)
    y_high=CalculatePenumbra("d",d_il,y_il)
    # Lateraler offset, inline crossline
    latOffset_il=CalculateLateralOffset(d_il, y_il)
    latOffset_cl=CalculateLateralOffset(d_cl, x_cl)
    
    # Korrektur der Profile um laterales Offset! --> fuer uniformity region
    y_il = y_il - latOffset_il
    x_cl = x_cl - latOffset_cl
    # uniformity region inline    
    uniformity_il_coord_p=CalculatePos('lat',"p",50,d_il,y_il)+\
                                            2*(y_low+y_high)/2
    uniformity_il_coord_d=CalculatePos('lat',"d",50,d_il,y_il)-\
                                            2*(y_low+y_high)/2
    # uniformity region crossline
    uniformity_cl_coord_p=CalculatePos('lat',"p",50,d_cl,x_cl)+\
                                           2*(float(x_low)+float(x_high))/2
    uniformity_cl_coord_d=CalculatePos('lat',"d",50,d_cl,x_cl)-\
                                           2*(float(x_low)+float(x_high))/2
                                        
    # Einschraenken der Uniformity Region fuer woechentliche und monatliche QA, DS8
    if typeQA == "Wo":
        if option[0:3] == "DS8":
            uniformity_il_coord_p = -55.0
            uniformity_il_coord_d =  55.0
            uniformity_cl_coord_p = -55.0
            uniformity_cl_coord_d =  55.0
            
    if typeQA == "Mo":
        if option[0:3] == "DS8":
            uniformity_il_coord_p = -66.1
            uniformity_il_coord_d =  66.1
            uniformity_cl_coord_p = -66.1
            uniformity_cl_coord_d =  66.1
            
    if typeQA == "TaDDLP":
        if option[0:3] == "DS8":
            uniformity_il_coord_p = -60.0
            uniformity_il_coord_d =  60.0
            uniformity_cl_coord_p = -60.0
            uniformity_cl_coord_d =  60.0
   
    # Glaettung der Dosis (inline, crossline) zur Bestimmung der Flatness und Symmetrie
    smind = 2 #Glaettungsfaktor
    
    d_cl_smoothed,cutoff_cl = GetSmoothedCurve(smind, x_cl, d_cl, uniformity_cl_coord_p, uniformity_cl_coord_d)
    d_il_smoothed,cutoff_il = GetSmoothedCurve(smind, y_il, d_il, uniformity_il_coord_p, uniformity_il_coord_d)
    return d_cl_smoothed, d_il_smoothed, cutoff_cl, cutoff_il


#------------------------------------------------------------------------------           
def Calculate_DepthDoseWaterPhantom(d_norm, coords):   
    ''' Calculates Depth Dose Parameters of Water Phantom Measurement
        prox R90, prox R95, prox R98, dist R10, dist R20, dist R50, dist R80
        dist R90, dist R95, SOBP: prox R90 – dist R90, SOBP: prox R98 – dist R90
        SOBP: prox R95 – dist R95, DFO, Uniformity: (Max-Min)/Max
    
    Parameters
    -------------------------
    d_norm: normiertes Dosisarray
    coords: array of coordinates belonging to the doses (z axis corr with WET of chamber)
    
    Returns
    -------------------------
    Array of Parameters in output order
    prox R90 / mm', 'prox R95 / mm', 'dist R90 / mm','dist R95 / mm','SOBP p-dR90 / mm','SOBP p-dR95 / mm','DFO / mm
    prox R98 / mm', 'dist R10 / mm', 'dist R20 / mm',distR50, 'dist R80 / mm', 'SOBP p98-dR90 / mm','Uniformity

    '''
    ## porximal values
    pR90=CalculatePos('dep',"p",90,d_norm,coords)                      
    pR95=CalculatePos('dep',"p",95,d_norm,coords)                      
    pR98=CalculatePos('dep',"p",98,d_norm,coords)                      
    # distal values
    dR10=CalculatePos('dep',"d",10,d_norm,coords)                      
    dR20=CalculatePos('dep',"d",20,d_norm,coords)                      
    dR50=CalculatePos('dep',"d",50,d_norm,coords)                      
    dR80=CalculatePos('dep',"d",80,d_norm,coords)                      
    dR90=CalculatePos('dep',"d",90,d_norm,coords)                      
    dR95=CalculatePos('dep',"d",95,d_norm,coords)                      
    #SOBPs
    sobp9090=CalculatePos('dep',"pd",90,d_norm,coords)                      
    sobp9890=dR90-pR98  
    sobp9595 = dR95-pR95                   
    #distal fall off
    dfo=dR20-dR80
    # uniformity
    minVal, maxVal, sym = CalculateMinMax(pR98+dfo,dR95-dfo,
                                              d_norm, coords)
    minMax=round((np.amax(d_norm)-minVal)/np.amax(d_norm),2)

    return [pR90, pR95, dR90, dR95, sobp9090,  sobp9595, dfo, pR98, dR10, dR20, dR50, dR80,sobp9890,minMax]
    
