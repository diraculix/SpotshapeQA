import sys
sys.path.insert(1, r' N:\fs4-HPRT\HPRT-Docs\Lukas\BA_Anna_Leimbach\Programm_V3.1.6_Aeneis_Kopie')
import spotshape_extractor

lynx_file = spotshape_extractor.FindLYNXFile('20230314', 90., 100.)
x_fwhm, y_fwhm = spotshape_extractor.AnalyzeSpotsGantryAngleDependent(lynx_file, option='gauss')  # other option 'fifty'
print(x_fwhm, y_fwhm)