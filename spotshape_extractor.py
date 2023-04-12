import numpy as np
import os.path, sys, time
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from PIL import Image, ImageDraw
from CalculateClass_ForQA import CalculatePos
from NormalizeClass_ForQA import NormalizeCenter100


def GrayScale(arr):     
    # normalize array to grayscale maximum (255)
    arr_greyscale = 1 * arr * 255 / arr.max()
    
    return arr_greyscale

def Gauss1D(x, *p):
    A, mu, sigma = p
    return A * np.exp(-(x - mu)**2 / (2. * sigma**2))


def FindLYNXFile(date, gantry, energy):
    pass


def LoadFromLYNXFile(inputf, GetSensor=False):
    # obtain spot intensity information from LYNX output .txt files
    try:
        data_arr = np.loadtxt(inputf, delimiter = '\t', skiprows = 33)
        data_arr = data_arr[::-1,:]
        y_il = data_arr[:,0]
        data_arr = data_arr[:,1:]
        xlen = len(data_arr[0,:])
        x_arr = np.loadtxt(inputf, delimiter = '\t', skiprows = 32, 
                            usecols=tuple(np.arange(1,xlen+1)), dtype = float)
        x_cl = x_arr[0,:]
        
        if GetSensor == True:           # read sensor type from lynx file
            lynx = np.loadtxt(inputf, delimiter = '\t', skiprows = 21, usecols=(0,), dtype = str)[0]
            lynx = lynx[9:]
        
        # normalize spotmap
        array_norm, normval = NormalizeCenter100(data_arr, x_cl, y_il)
                    
        # find center position
        x_clarr, y_ilarr = np.meshgrid(x_cl, y_il)
        d_cl = array_norm[y_ilarr == 0]
        d_il = array_norm[x_clarr == 0]                

        # richtung immer von - nach plus
        flip_x = 'n'                    # Abfragen, ob Matrix gespiegelt wurden, um sie ggf. zurückzuspiegeln
        flip_y = 'n'
        # ACHTUNG: das ist nicht IEC konform, da y (oben nach unten) eigentlich von + nach - gehen muesste        
        if x_cl[0]>x_cl[-1]:
            x_cl=x_cl[::-1]
            d_cl=d_cl[::-1]
            flip_x = 'j'
        if y_il[0]>y_il[-1]:
            y_il=y_il[::-1]
            d_il=d_il[::-1]
            flip_y = 'j'
        
    except IOError:
        print("Please correct your inputdata-path! ")

    if GetSensor == True:
        return array_norm,x_cl,y_il,d_cl,d_il, flip_x, flip_y, lynx
    else:
        return array_norm,x_cl,y_il,d_cl,d_il, flip_x, flip_y


def Calculate_Spots_GantryAngleDependent(dataarray, xarray, yarray, option='fifty'):   
    # find 9 QA spots and calculate position/fwhm

    # initialize empty dicts 
    FWHM_Dict = {}
    SpotPos_Dict = {}
    Intensity_Dict = {}
    KeyNames = []
    
    # cut out spots in 100x100px tiles from original 600x600px LYNX array
    SubDoseArray_Dict = {"SpotMitte" : dataarray[250:350, 250:350], "SpotUntenLinks" : dataarray[10:110, 10:110], "SpotUntenRechts" : dataarray[10:110, -110:-10], "SpotObenLinks" : dataarray[-110:-10, 10:110], "SpotObenRechts" : dataarray[-110:-10, -110:-10], "SpotLinks" : dataarray[250:350, 10:110], "SpotRechts" : dataarray[250:350, -110:-10], "SpotUnten" : dataarray[10:110, 250:350], "SpotOben" : dataarray[-110:-10, 250:350]} #erster Wert: oben-unten, !Lynx-Bild ist geflipt!
    LeftLowerCorner_Dict = {"SpotMitte" : [250,250], "SpotUntenLinks" : [10,10], "SpotUntenRechts" : [10,490], "SpotObenLinks" : [490,10], "SpotObenRechts" : [490,490], "SpotLinks" : [250,10], "SpotRechts" : [250,490], "SpotUnten" : [10,250], "SpotOben" : [490,250]}
    xarray_dict = {"SpotMitte" : xarray[250:350], "SpotUntenLinks" : xarray[10:110], "SpotUntenRechts" : xarray[490:590], "SpotObenLinks" : xarray[10:110], "SpotObenRechts" : xarray[490:590], "SpotLinks" : xarray[10:110], "SpotRechts" : xarray[490:590], "SpotUnten" : xarray[250:350], "SpotOben" : xarray[250:350]}
    yarray_dict = {"SpotMitte" : yarray[250:350], "SpotUntenLinks" : yarray[10:110], "SpotUntenRechts" : yarray[10:110], "SpotObenLinks" : yarray[490:590], "SpotObenRechts" : yarray[490:590], "SpotLinks" : yarray[250:350], "SpotRechts" : yarray[250:350], "SpotUnten" : yarray[10:110], "SpotOben" : yarray[490:590]}

    # iterate through QA spots
    for key in SubDoseArray_Dict:
        Spot = SubDoseArray_Dict[key]
                
        # single spots at fixed positions                        
        # normalize each spot max pixel to 100%
        Spot_norm = 100.*Spot/Spot.max()
        # find maximum position
        Spot_MaxPos = np.where(Spot_norm >99.999)
        Intensity = np.average(Spot[int(np.average(Spot_MaxPos[0]))-1:int(np.average(Spot_MaxPos[0]))+2, int(np.average(Spot_MaxPos[1]))-1:int(np.average(Spot_MaxPos[1]))+2]) # array mit 1 links von max int und 1 rechts (muss +2 sein wegen array def in python: letzter Wert der angegeben ist wird nihct genommen!)               
        
        # go on only, if the spot is really there: Intensity is not nan!                
        if str(Intensity) != "nan":  
            if option == 'fifty':
                # Normierung des Spots auf Intensitaets Wert
                Spot_norm = 100.*Spot/Intensity

                # Profiles through the 100 position (assumed center): if more than one: average position, put to integer, "closest" pixel
                xProfil = Spot_norm[int(np.average(Spot_MaxPos[0])), :]
                yProfil = Spot_norm[:, int(np.average(Spot_MaxPos[1]))]

                # FWHM calc in profile = 50-50 field size
                FWHM_Spot_x = CalculatePos("PBS","pd", 50, xProfil, xarray_dict[key])
                FWHM_Spot_y = CalculatePos("PBS","pd", 50, yProfil, yarray_dict[key])

                # FWHM at position of profiles                
                FWHM_Dict[key] = [FWHM_Spot_x, FWHM_Spot_y]
                SpotPos_Dict[key] = [np.average(yarray[LeftLowerCorner_Dict[key][1]+Spot_MaxPos[1]]), np.average(xarray[LeftLowerCorner_Dict[key][0] + Spot_MaxPos[0]])]
                Intensity_Dict[key] = Intensity
                KeyNames.append(key) 
            
            elif option == 'gauss':
                # obtain integrated 1D profiles by summation
                x_profile_integrated = Spot.sum(axis=0)
                max_index_x, max_index_y = int(np.mean(np.argmax(Spot, axis=0))), int(np.mean(np.argmax(Spot, axis=1)))
                # x_profile_integrated = Spot[:, max_index_x]  # that would mean exact Gaussian fitting through spot maximum
                x_axis = [(px + 1) * 0.5 - 0.25 for px in range(len(x_profile_integrated))]  # mm
                y_profile_integrated = Spot.sum(axis=1)
                # y_profile_integrated = Spot[:, max_index_y]
                y_axis = [(px + 1) * 0.5 - 0.25 for px in range(len(y_profile_integrated))]
                
                # fit 1D Gaussian with initial guess parameters
                popt_x, _ = curve_fit(Gauss1D, x_axis, x_profile_integrated, p0=[1., 0., 1.])  # popt contains optimization parameters A, mu, sigma
                popt_y, _ = curve_fit(Gauss1D, y_axis, y_profile_integrated, p0=[1., 0., 1.])

                # calculate FWHM from sigma
                const = 2 * np.sqrt(2 * np.log(2))
                x_fwhm = const * abs(popt_x[-1])
                y_fwhm = const * abs(popt_y[-1])
                
                # plt.plot(x_axis, x_profile_integrated, 'o', color='tab:blue', ms=2)
                # plt.plot(x_axis, Gauss1D(x_axis, *popt_x), color='tab:blue')
                # plt.plot(y_axis, y_profile_integrated, 'o', color='tab:orange', ms=2)
                # plt.plot(y_axis, Gauss1D(y_axis, *popt_y), color='tab:orange')
                # plt.show()
                # return None

                SpotPos_Dict[key] = [np.average(yarray[LeftLowerCorner_Dict[key][1]+Spot_MaxPos[1]]), np.average(xarray[LeftLowerCorner_Dict[key][0] + Spot_MaxPos[0]])]
                Intensity_Dict[key] = Intensity
                FWHM_Dict[key] = [x_fwhm, y_fwhm]
                KeyNames.append(key)
    
    return FWHM_Dict, SpotPos_Dict, Intensity_Dict, KeyNames


def AnalyzeSpotsGantryAngleDependent(path):
    # analyze spots extracted from LYNX output file    
    array_norm,x_cl,y_il, d_cl, d_il, flip_x, flip_y= LoadFromLYNXFile(path)

    # wenn geflippt wurde (die richtung x, y immer von - nach plus), muss auch bildMatrix noch geflippt werden
    if flip_x == "j":        
        arr_norm=array_norm[:,::-1]
    if flip_y == "j":
        arr_norm=array_norm[::-1,:]
    
    FWHM_Dict, SpotPos_Dict, Intensity_Dict, KeyNames = Calculate_Spots_GantryAngleDependent(arr_norm, x_cl, y_il, option='gauss')
    
    # KeyNames in same order as in tables:
    KeyNames = ("SpotMitte", "SpotUntenLinks", "SpotUntenRechts", "SpotObenLinks", "SpotObenRechts", "SpotLinks", "SpotRechts", "SpotUnten", "SpotOben")
    for i in range(len(KeyNames)):
        # spot position
        x_pos, y_pos = SpotPos_Dict[KeyNames[i]][0], SpotPos_Dict[KeyNames[i]][1]
        # spot fwhm
        x_fwhm, y_fwhm = FWHM_Dict[KeyNames[i]][0], FWHM_Dict[KeyNames[i]][1]
        # spot intensity
        intensity = Intensity_Dict[KeyNames[i]]

        # ... append data to array/dataframe ...                    
    
    print('Exit code 0')

if __name__ == '__main__':
    testfile = r'Q:\PHYSIK\PROTONEN\QA\QA_Halbjaehrlich\Spotmessungen_winkelabhaengig\2023-03\Spots_G000_100MeV_2023-03-14.txt'
    AnalyzeSpotsGantryAngleDependent(testfile)