#!/usr/bin/env python
# -*- coding: latin-1 -*-
from pcraster import *



def assertWithinRange(map, Lower, Upper):
   """ Checks the range of maps
   map Input PCRaster map"""
   Minimum = cellvalue(mapminimum(map), 0, 0)
   Maximum = cellvalue(mapmaximum(map), 0, 0)
   assert Minimum[0] >= Lower and Maximum[0] <= Upper

def writeVar(varName, varContent):
   """Writes the value of a variable to a file
   varName string name of variable
   varContent value of variable"""
   
   varContent = str(varContent)
   checkFile.write(varName+" = "+varContent+"\n")
   return

def writeLoc(varName, map):
   """Writes values of a variable at a specific location to a file
   varName string name of variable
   map name of map"""
   
   findLoc = cellvalue(map, rowy, colx)
   map = findLoc[0]
   writeVar(varName, map)
   return

def median(map):
   """Function to calculate median of a map

   map Input PCRaster map"""

   OrderMap = order(map)
   Mid = roundoff(mean(OrderMap))
   MidMap = ifthenelse(OrderMap == Mid, map, 0)
   Median = cellvalue(mapmaximum(MidMap), 0, 0)
   assert Median[0] > 0.0
   return Median[0]

def mean(map):
   """Calculates the mean value of a PCRaster map

   map Input PCRaster map"""
   
   Total = cellvalue(maptotal(map), 0, 0)
   NumCells = cellvalue(maparea(map) / cellarea(), 0, 0)
   assert NumCells[0] != 0
   mean = Total[0] / NumCells[0]
   return mean

# SEBS functions
def Rswd(DEM, Lat, Trans, DOY, Time):
   """ Potential Radiation Equator model
   (c) O. van Dam, UU, Tropenbos-Guyana
   Version 5, June 2000
   NOTE: Copyright: This program is free to use provided·
         you refer to the manualfor citation.
         Do not distribute without prior approval of the author.
         Manual and additional info: O.vanDam@geog.uu.nl

   -----------------------------------------------------
                   Model for calculation
               incoming potential light energy
   -----------------------------------------------------
   
   DEM Input Digital Elevation Model (spatial)
   Lat Latitude in decimal degrees (non-spatia)
   Trans Transmissivity tau (Gates, 1980) (non-spatial)
   DOY Day of Year (non-spatial)
   Time Time in hours (non-spatial)"""
   
   # constants
   pi       = 3.1415          # pi
   Sc       = 1367.0          # Solar constant (Gates, 1980) [W/m2]

   SlopMap = scalar(atan(slope(DEM)))
   AspMap  = scalar(aspect(DEM)) # aspect [deg]
   AtmPcor = ((288.0-0.0065*DEM)/288.0)**5.256 # atmospheric pressure correction [-]

   # Solar geometry
   # ----------------------------
   # SolDec  :declination sun per day  between +23 & -23 [deg]
   # HourAng :hour angle [-] of sun during day
   # SolAlt  :solar altitude [deg], height of sun above horizon
   SolDec  = -23.4*cos(360.0*(DOY+10.0)/365.0)
   HourAng = 15.0*(Time-12.01)
   SolAlt  = scalar(asin(scalar(sin(Lat)*sin(SolDec)+cos(Lat)*cos(SolDec)*cos(HourAng))))

   # Solar azimuth
   # ----------------------------
   # SolAzi  :angle solar beams to N-S axes earth [deg]
   SolAzi = scalar(acos((sin(SolDec)*cos(Lat)-cos(SolDec)*sin(Lat)*cos(HourAng))/cos(SolAlt)))
   SolAzi = ifthenelse(Time <= 12.0, SolAzi, 360.0 - SolAzi)
   # Additonal extra correction by R.Sluiter, Aug '99
   SolAzi = ifthenelse(SolAzi > 89.994 and SolAzi < 90.0, 90.0, SolAzi)
   SolAzi = ifthenelse(SolAzi > 269.994 and SolAzi < 270.0, 270.0, SolAzi)

   # Surface azimuth
   # ----------------------------
   # cosIncident :cosine of angle of incident; angle solar beams to angle surface
   cosIncident = sin(SolAlt)*cos(SlopMap)+cos(SolAlt)*sin(SlopMap)*cos(SolAzi-AspMap)

   # Critical angle sun
   # ----------------------------
   # HoriAng  :tan maximum angle over DEM in direction sun, 0 if negÂ·
   # CritSun  :tan of maximum angle in direction solar beams
   # Shade    :cell in sun 1, in shade 0
   HoriAng   = horizontan(DEM,directional(SolAzi))
   HoriAng   = ifthenelse(HoriAng < 0.0, scalar(0.0), HoriAng)
   CritSun   = ifthenelse(SolAlt > 90.0, scalar(0.0), scalar(atan(HoriAng)))
   Shade   = ifthenelse(SolAlt > CritSun, scalar(1), scalar(0))
   
   # Radiation outer atmosphere
   # ----------------------------
   OpCorr = Trans**((sqrt(1229.0+(614.0*sin(SolAlt))**2.0)-614.0*sin(SolAlt))*AtmPcor)     # correction for air masses [-]Â·
   Sout   = Sc*(1.0+0.034*cos(360.0*DOY/365.0)) # radiation outer atmosphere [W/m2]
   Snor   = Sout*OpCorr                    # rad on surface normal to the beam [W/m2]

   # Radiation at DEM
   # ----------------------------
   # Sdir   :direct sunlight on a horizontal surface [W/m2] if no shade
   # Sdiff  :diffuse light [W/m2] for shade and no shade
   # Stot   :total incomming light Sdir+Sdiff [W/m2] at Hour
   # PotRad :avg of Stot(Hour) and Stot(Hour-HourStep)
   Sdir   = ifthenelse(Snor*cosIncident*Shade < 0.0, 0.0, Snor*cosIncident*Shade)
   Sdiff  = ifthenelse(Sout*(0.271-0.294*OpCorr)*sin(SolAlt) < 0.0, 0.0, Sout*(0.271-0.294*OpCorr)*sin(SolAlt))
   #Rswd  = Sdir + Sdiff                                         # Rad [W/m2]
   Rswd = Snor
   return Rswd

def LAINDVI(NDVI):
   """ Calculates initial Leaf Area Index from NDVI (Su, 1996). Output is non-spatial
   
   NDVI Input Normalized Difference Vegetation Index Map (scalar, ratio between 0 and 1)"""

   nd_max = cellvalue(mapmaximum(NDVI), 0, 0)
   nd_min = cellvalue(mapminimum(NDVI), 0, 0)
   nd_mid = median(NDVI)
   nd_df = nd_max[0] - nd_min[0]
   if nd_df == 0.0:
      nd_df == 1.0

   LAI = sqrt(nd_mid * (1.0 + nd_mid) / (1.0 - nd_mid + 1.0E-6))
   if LAI > 6.0:
      LAI == 6.0
   return LAI, nd_max[0], nd_min[0], nd_mid, nd_df

def u_pbl(NDVI):
   """Calculates Planetary Boundary Layer wind speed [m s-1] from NDVI
   
   NDVI Input PCRaster NDVI map (scalar, ratio between 0 and 1)"""
   
   z0m = 0.005 + 0.5 * (nd_mid/nd_max) ** 2.5
   assert z0m >= 0.0
   fc = ((nd_mid - nd_min) / nd_df) ** 2.0    # fractional vegetation cover == Wfol (-)
   assert fc >= 0.0
   h = z0m / 0.136                            # total height of vegetation (m)
   d = 2.0/3.0 * h			      # zero plane displacement (m)
   u_c = ln((z_pbl - d) / z0m) / ln((z_ms - d) / z0m)
   u_pbl = u_s * u_c
   return u_pbl, z0m, d, fc, h

# FUNCTIONS FOR DETERMINATION OF ROUGHNESS LENGTH FOR HEAT TRANSFER
def FKB_1(u_zref, zref, h, LAI, Wfol, Ta, pa):
   """Initial determination of roughness length for heat transfer (non-spatial)
   KB-1 function according to Massman, 1999
   Convention of variable names:
   f_z = f(z)
   d2h = d/h
   
   u_zref Input wind speed at reference height [m s-1]
   zref Input reference height [m]
   h Input canopy height [m]
   LAI Input canopy total Leaf Area Index [-]
   Wfol Input Fractional canopy cover [-]
   Ta Input ambient temperature [degrees Celsius]
   pa Input ambient air pressure [Pa]"""

   # Constants
   C_d = 0.2   # foliage drag coefficient
   C_t = 0.01  # heat transfer coefficient
   k = 0.41     # Von Karman constant
   Pr = 0.7    # Prandtl number
   hs = 0.009  # height of soil roughness obstacles (0.009-0.024)
   
   # Calculations
   Wsoil = 1.0 - Wfol
   if Wfol == 0.0: # for bare soil take soil roughness
      h = hs
   assert Wfol >= 0.0 and Wfol <= 1.0 and Wsoil >= 0.0 and Wsoil <= 1.0
   z0 = 0.136 * h   # Brutsaert (1982)
   u_h0 = u_zref * ln(2.446) / ln ((zref - 0.667 * h) / z0) # wind speed at canopy height
   u_h0 = cellvalue(u_h0, 0, 0)
   u_h0 = u_h0[0]
   assert u_h0 >= 0.0
   ust2u_h = 0.32 - 0.264/exp(15.1 * C_d * LAI)
   ustarh = ust2u_h * u_h0
   nu0 = 1.327E-5 * (101325.0 / pa) * (Ta / 273.15 + 1.0) ** 1.81 # kinematic viscosity
   n_h = C_d * LAI / (2.0 * ust2u_h ** 2.0)
   # First term
   if n_h != 0.0:
      F1st = k * C_d / (4.0 * C_t * ust2u_h * (1.0 - exp(pcrumin(n_h)/2.0))) * Wfol ** 2.0
   else:
      F1st = 0.0
   # Second term
   S2nd = k * ust2u_h * 0.136 * Pr ** (2.0/3.0) * sqrt(ustarh * h / nu0) * Wfol ** 2.0 * Wsoil ** 2.0
   # Third term
   T3rd = (2.46 * (u_zref * k / ln(zref/hs) * hs / nu0) ** 0.25 - ln(7.4)) * Wsoil ** 2.0
   
   return F1st + S2nd + T3rd

def z0h(KB_1, z0m):
   """Calculates the scalar roughness height for heat transfer (z0h)
   KB_1 Input KB_1 values
   z0m Input scalar roughness height for momentum"""
   
   z0h = z0m / exp(KB_1)
   return z0h

def GKB_1(u_zref, zref, h, LAI, Wfol, Ta, pa):
   """Same as FKB_1, but then for spatial in- and output"""
   
   # Constants
   C_d = 0.2   # foliage drag coefficient
   C_t = 0.05  # heat transfer coefficient
   k = 0.41     # Von Karman constant
   Pr = 0.7    # Prandtl number
   hs = 0.009  # height of soil roughness obstacles (0.009-0.024)
   
   # Calculations
   Wsoil = 1.0 - Wfol
   h = ifthenelse(Wfol == 0.0, hs, h)
   
   z0 = 0.136 * h   # Brutsaert (1982)
   u_h0 = u_zref * ln(2.446) / ln ((zref - 0.667 * h)/z0) # wind speed at canopy height
   ust2u_h = 0.32 - 0.264 / exp(15.1 * C_d * LAI)
   ustarh = ust2u_h * u_h0
   nu0 = 1.327E-5 * (101325.0/pa) * (Ta / 273.15 + 1.0) ** 1.81 # kinematic viscosity
   n_h = C_d * LAI / (2.0 * ust2u_h ** 2.0)
   # First term
   F1st = ifthenelse(pcrne(n_h, 0.0), k * C_d / (4.0 * C_t * ust2u_h * (1.0 - exp(pcrumin(n_h)/2.0))) * Wfol ** 2.0, 0.0)
   # Second term
   S2nd = k * ust2u_h * 0.136 * Pr ** (2.0/3.0) * sqrt(ustarh * h / nu0) * Wfol ** 2.0 * Wsoil ** 2.0
   # Third term
   T3rd = (2.46 * (u_zref * k / ln(zref/hs) * hs / nu0) ** 0.25 - ln(7.4)) * Wsoil ** 2.0
   
   return F1st + S2nd + T3rd

def Rn(Alfa, Rswd, Eair, t_pbl, ems, T):
#def Rn(Alfa, Rswd, Eair, t_pbl, Eground):
   """ Calculation of surface net radiation [W m-2]
   Alfa Input albedo map [-]
   Rswd Input downward solar radiation [W m-2], PCRaster map from POTRAD
   Eair Input emissivity air [-]
   Eground Input PCRaster emissivity map [-]
   t_pbl Input PBL temperature map [K]
   T Surface Kinetic Temperature [K]"""
   
   print("Calculating Net Radiation map...")
   # constants
   sigma = 5.678E-8   #Stefan Boltzmann's constant (W m-2 K-4)
   
   # calculations
   Rn = (1.0 - Alfa) * Rswd + 5.678 * ems * (Eair * (t_pbl/100.0)**4.0 - (T/100.0)**4.0)
   return Rn

def G0(Rn, cover):
   """Calculates Soil Heat Flux [W m-2]
   Rn Input Surface Net Radiation [W m-2]
   cover Input fractional canopy cover [-]"""
   
   print("Calculating soil heat flux map...")
   # constants:
   Gamma_c = 0.05   # ratio of G0 to Rn for full vegetation canopy (Monteith, 1973)
   Gamma_s = 0.315  # ratio of G0 to Rn for bare soil (Kustas & Daughtry, 1989)
   
   # calculation
   G0 = Rn * (Gamma_c + (1.0 - cover) * (Gamma_s - Gamma_c))
   return G0

def FRUstar(z_pbl,hst):
   """Iteration to calculate RUstar
   z_pbl Input PBL depth [m]
   hst Input height of the ASL [m]"""
   
   print("Starting iterations to derive stability parameters..." )
   
   RUstar = ku / zdm
   RH = CH * RUstar / zdh
   RH0 = RH
   Reps = 10.0
   Isteps = 0
   RHA = RH
   RHB = RH
   RH0A = RH0
   RH0B = RH0
   RUstarA = RUstar
   RUstarB = RUstar
   IstepsA = Isteps
   IstepsB = Isteps
   RepsA = Reps
   RepsB = Reps
   itNr = 100.0
   itThreshold = 0.01
     
   while RepsA > itThreshold and IstepsA < itNr:
      RLA = CL * RUstarA ** 3.0 / RHA
      tempBw = Bw(z_pbl, RLA, z0m)
      RUstarA = ku / (zdm - tempBw)
      tempCw = Cw(z_pbl, RLA, z0m, z0h)
      RHA = CH * RUstarA / (zdh - tempCw)
      RepsA = mapmaximum(abs(RH0A - RHA))
      difa = abs(RH0A - RHA)
      min = mapminimum(difa)
      meandif = mean(difa)
      RH0A = RHA
      IstepsA = IstepsA + 1
      percentage = (IstepsA/itNr)*100
      print("Iteration A:", int(percentage), "% completed\r",)
   print()
   
   while RepsB > itThreshold and IstepsB < itNr:
      RLB = CL * RUstarB ** 3.0 / RHB
      tempPSIm_y1 = PSIm_y(zd0/ RLB)
      tempPSIm_y2 = PSIm_y(z0m / RLB)
      RUstarB = ku / (zdm - tempPSIm_y1 + tempPSIm_y2)
      tempPSIh_y1 = PSIh_y(zd0 / RLB)
      tempPSIh_y2 = PSIh_y(z0h / RLB)
      RHB = CH * RUstarB / (zdh - tempPSIh_y1 + tempPSIh_y2)
      RepsB = mapmaximum(abs(RH0B - RHB))
      difb = abs(RH0B - RHB)
      meandif = mean(difb)
      min = mapminimum(difb)                                                    
      RH0B = RHB
      IstepsB =  IstepsB + 1
      percentage = (IstepsB/itNr)*100
      print("Iteration B:", int(percentage), "% completed\r",)
   print() 
   RUstar = ifthenelse(z_pbl >= hst, RUstarA, RUstarB)
   RL = ifthenelse(z_pbl >= hst, RLA, RLB)
   dif = ifthenelse(z_pbl >= hst, difa, difb)
   return RUstar, RL

# MOS STABILITY CORRECTION FUNCTIONS
def PSIma(f, g):
   a = 0.33
   b = 0.41
   pi = 3.141592654
   tangens = scalar(atan((2.0 * g - 1.0) / sqrt(3.0))) * pi /180
   tangens = ifthenelse(tangens > pi/2.0, tangens - 2.0 * pi, tangens) 
   PSIma = ln(a + f) - 3.0 * b * f ** (1.0 / 3.0) + b * a ** (1.0 / 3.0) / 2.0 * ln((1 + g) ** 2.0 / (1.0 - g + sqr(g))) + sqrt(3.0) * b * a ** (1.0 / 3.0) * tangens
   return PSIma
   
def PSIm_y(Y):
   # Integrated stability correction function for momentum
   # Inputs
   # Y = -z/L, where z is the height, L the Obukhov length
   # test values
   
   # Constants (Brutsaert, 1999)
   a = 0.33
   b = 0.41
   m = 1.0
   pi= 3.141592654
   
   # Calculation
   #//HK 040902 
   Y = abs(Y) #abs(Y)
   x = (Y/a) ** (1.0/3.0)
   PSI0 = pcrumin(ln(a)) + sqrt(3.0) * b * a ** (1.0 / 3.0) * pi / 6.0
   b_3 = b ** -3.0
   PSIm_y = ifthenelse(Y <= b_3, PSIma(Y, x) + PSI0, PSIma(b_3, ((b_3/a)**(1.0/3.0))) + PSI0)
   #PSIm_y = ifthenelse(Y <= b_3, PSIma(Y, x) + PSI0, (1.0 / (PSIma(b_3, ((b_3/a)**(1.0/3.0))))) + PSI0)
   return PSIm_y
   
def PSIh_y(Y):
   # Integrated stability correction function for heat
   # Inputs
   # Y = -z/L, z is the height, L the Obukhov length
   # constants (Brutsaert, 1999)
   c = 0.33
   d = 0.057
   n = 0.78
   # Calculation
   Y =  abs(Y)
   PSIh_y = (1.0 - d) / n * ln((c + Y ** n) / c)
   return PSIh_y

# BAS STABILITY CORRECTION FUNCTIONS
def Bw(hi, L, z0):
   # constants (Brutsaert, 1999)
   alfa = 0.12
   beta = 125.0
   
   # calculations
   B0 = (alfa / beta) * hi
   B1 = -1.0 *z0 / L
   B11 = -alfa * hi / L
   B21 = hi / (beta * z0)
   B22 = -beta * z0 / L
   tempB11 = PSIm_y(B11)
   tempB1 = PSIm_y(B1)
   B = ifthenelse(z0 < B0, -1.0 * ln(alfa) + PSIm_y(B11) - PSIm_y(B1), ln(B21) + PSIm_y(B22) - PSIm_y(B1))
   Bw = ifthenelse(B < 0.0, 0.0, B) # This results from unfortunate parameter combination!
   return Bw

def Cw(hi, L, z0, z0h):
   alfa = 0.12
   beta = 125.0
   C0 = (alfa / beta) * hi
   C1 = pcrumin(z0h) / L
   C11 = -alfa * hi / L
   C21 = hi / (beta * z0)
   C22 = -beta * z0 / L
   C = ifthenelse(z0 < C0, pcrumin(ln(alfa)) + PSIh_y(C11) - PSIh_y(C1), ln(C21) + PSIh_y(C22) - PSIh_y(C1))
   Cw = ifthenelse(C < 0.0, 0.0, C) # This results from unfortunate parameter combination!
   return Cw

def esat(t):
   """Calculation of saturated vapour pressure [Pa]
   
   t Input temperature in degrees Celsius"""
   # constants
   e0 = 610.7   # saturated water vapour pressure at 273.15K
   A = 7.5
   B = 237.3
   
   # Calculation
   esat = e0 * 10.0 ** ((A * t) / (B + t))
   return esat


#----------------------------------------------------------------------------
# INPUT

# Validation pixel
#rowy = 0				# row number of validation pixel
#colx = 0				# column number of validation pixel
#checkFile = file("check.txt", "w")	# name of validation textfile

# Define inputs
# maps

DEM = readmap('./example/dem90.map')	# Digital Elevation Model [m]
nd = readmap('./example/ndvi90.map') # NDVI map [-]
T = readmap('./example/tkin90.map') # Surface temperature [Kelvin]
albedo = readmap('./example/albedo90.map') # Albedo map [-]
ems = readmap('./example/emissivity90.map') # emissivity [-]

# parameters

Trans = 0.788606 #transmissivity [0-1]
Lat = 33.9932 #Latitude [dd]
DOY = 294.0 #Day of year
Time = 11.217 #Time of overpass [decimal hours]
z_pbl = 1000.0#PBL height [m]
alt_ms = 2.5 #Measurement height [m]
u_s = 4.313 #Wind speed [m/s]
t_s = 27.35 #Air temperature [Celcius]
p_s = 100000.0 #Air pressure [Pa]
hr_s = 0.5055 #Relative humidity [0-1]
z_ms = alt_ms

# Define output files
lemap = './example/le.map'
hmap = './example/h.map'
gmap = './example/g0.map'
rnmap = './example/rn.map'
evaprmap = './example/evapr.map'
evapfrmap = './example/evapfr.map'
etmap = './example/et.map'

print("Initializing SEBS.",)
# Initialize model starttime for calculation runtime
starttime = time()
# Check input data
nd = ifthenelse(pcror(pcrlt(nd,0.0),pcrgt(nd,1.0)), 1.0, nd) # Convert waterbodies to 1.0 --> soilflux is minimal
assertWithinRange(nd, 0.0, 1.0)
minimumDEM = cellvalue(mapminimum(DEM), 0, 0)
assert minimumDEM[0] >= 0.0
minimumT = cellvalue(mapminimum(T), 0, 0)
assert minimumT[0] >= 0.0
assert DOY >= 0.0 and DOY <= 366
assert Time >= 0.0 and Time <= 24.0
assert alt_ms >= 0.0
assert u_s >= 0.0
assert hr_s >= 0.0 and hr_s <= 1.0
assert z_pbl >= 0.0
albedo = ifthenelse(pcror(pcrlt(albedo,0.0),pcrgt(albedo,1.0)), 0.0, albedo)
assertWithinRange(albedo, 0.0, 1.0)
ems = ifthenelse(pcror(pcrlt(ems,0.0), pcrgt(ems,1.0)), 0.0, ems)
assertWithinRange(ems, 0.0, 1.0)
T = ifthen(T >= 273.15, T)
print("\b.",)   

# INITIALIZE MODEL
# Calculating initial LAI
LAINDVI = LAINDVI(nd)
LAI = LAINDVI[0]
assert(LAI >= 0.0 and LAI <= 6.0)
print("\b.",)
nd_max = LAINDVI[1]
nd_min = LAINDVI[2]
nd_mid = LAINDVI[3]
nd_df = LAINDVI[4]
print("\b.",)
# Calculate initial PBL parameters
Fu_pbl = u_pbl(nd)
u_pbl = Fu_pbl[0]
u_pbl = cellvalue(u_pbl, 0, 0)
u_pbl = u_pbl[0]
assert u_pbl >= 0.0
print("\b.",)
z0m = Fu_pbl[1]
d = Fu_pbl[2]
fc = Fu_pbl[3]
h = Fu_pbl[4]

# Calculating initial KB-1 and z0h
KB_1 = FKB_1(u_pbl, z_pbl, h, LAI, fc, t_s, p_s)
KB_1 = cellvalue(KB_1, 0, 0)
KB_1 = KB_1[0]
print("\b.",)
z0h = cellvalue(z0h(KB_1, z0m), 0, 0)
z0h = z0h[0]
print("\b.")

# Calculating initial temperatures and pressures"
t_c = ln((z_pbl - d) / z0h) / ln((alt_ms - d) / z0h)
t_s = t_s + 273.15
t_pbl_A = T * (1.0 - t_c) + t_s * t_c
p_s_A = p_s * ((44331.0 - DEM) / (44331.0 - alt_ms)) ** (1.0 / 0.1903)   # surface pressure
z_pbl_A = z_pbl
p_pbl_A = p_s * ((44331.0 - (DEM + z_pbl_A)) / (44331.0 - alt_ms)) ** (1.0 / 0.1903)
helpvar1= DEM/44331.0
helpvar2 = 1.0 - helpvar1
T0 = T / helpvar2 ** 1.5029
t_pbl_A = t_pbl_A / (1.0 - DEM / 44331.0) ** 1.5029
T_0pbl = 0.5 * (T0 + t_pbl_A)   # mean potential temperature
Tcn = T_0pbl - 273.15    # mean PBL temperature converted to degrees Celcius
esat = 611.0 * exp(17.502 * Tcn / (Tcn + 240.97))   # Pa
hr_pbl = hr_s
eact = hr_pbl * mean(esat) # actual vapour pressure
q_pbl_A = 5.0 / 8.0 * eact / p_pbl_A
z_pbl = z_pbl_A
ps = p_s_A
Ta = T_0pbl - 273.15
t_pbl = t_pbl_A
LAI = sqrt(nd * (1.0 + nd)/ (1.0 + 1.0E-6 - nd))
LAI = ifthenelse(LAI > 6.0, 6.0, LAI)
#assert LAI >= 0.0
fc = ((nd - nd_min) / nd_df) ** 2.0
assertWithinRange(fc, 0.0, 1.0)
p_pbl = p_pbl_A
q_pbl = q_pbl_A
z0m = 0.005 + 0.5 * (nd / nd_max) ** 2.5
d = z0m * 4.9
h = z0m / 0.136
KB_1 = GKB_1(u_pbl, z_pbl, h, LAI, fc, Ta, p_pbl)
z0h = z0m / exp(KB_1)
Tsk = T					# potential surface temperature
Theta_s = T0
Theta_v = Tsk * (1.0 + 0.61 * q_pbl)	# surface virtual temperature
#Theta_a = t_pbl				# potential air temperature at reference height (K)
Theta_a = t_pbl * (101325/p_pbl) ** 0.286
T0ta = Theta_s - Theta_a
Rv = 461.05				# specific gas constant water vapour (J kg-1 K-1)
Rd = 287.04				# specific gas constant dry air (J kg-1 K-1)
Cp = 1005.0				# specific heat (J kg-1 K-1)
eact = p_pbl * q_pbl * (Rv / Rd)	# actual vapour pressure
rhoa = ps / (Rd * Theta_v)		# surface air density (kg m-3)
rhoam = (ps / (Rd * Tsk)) * (1.0 - 0.378 * eact / ps) # moist air density (kg m-3)
rhoacp = rhoa * Cp			# specific air heat capacity (J K-1 m­3)
alfa = 0.12
beta = 125.0
g = 9.81
k = 0.4
hst = max((alfa * z_pbl), (beta * z0m)) # height of ASL (m)
zd0 = z_pbl - d
ku = k * u_pbl
zdm = ln(zd0/z0m)
zdh = ln(zd0/z0h)
CH = T0ta * k * rhoacp
CL = pcrumin(rhoam) * Cp * Theta_v / (k * g)


# Calculate energy balance
print("Calculating Energy Balance terms...")
Rswd = Rswd(DEM, Lat, Trans, DOY, Time)
Eair = 9.2 * (t_pbl/1000.0) ** 2.0
Rn = Rn(albedo, Rswd, Eair, t_pbl, ems, T)
report(Rn, rnmap)
G0 = G0(Rn, fc)
report(G0,gmap)
R_G = Rn - G0
# Dry-limit heat flux
print("Calculating Dry Limit...")
H_d = R_G
FRUstar = FRUstar(z_pbl,  hst)
RUstar = FRUstar[0]
RL = FRUstar[1]
print("Calculating Wet Limit...")
# For completely wet areas
# Wet-limit stability length
L_e = 2.430E+06   # Latent heat of vapourization (J kg-1) (Brutsaert, 1982)
L_w = pcrumin(RUstar ** 3.0) * rhoam / (0.61 * k * g * R_G / L_e)
C_wet = ifthenelse(z_pbl >= hst, Cw(z_pbl, L_w, z0m, z0h), PSIh_y(pcrumin(z_pbl/L_w)))
# Wet-limit external resistance
re_w = (zdh - C_wet) / (k * RUstar)
re_w = ifthenelse(re_w <= 0.0, zdh / (k * RUstar), re_w)
# Wet-limit heat flux
slopef = 17.502 * 240.97 * esat / (Ta + 240.97) ** 2.0
gamma = 67.0   # psychrometric constant (Pa K-1)
H_w = (R_G - (rhoacp / re_w) * ((esat - eact) / gamma)) / (1.0 + slopef / gamma)
LEwet = Rn - G0 - H_w
# Sensible Heat flux
print("Calculating sensible heat flux...")
C_i = ifthenelse(z_pbl >= hst, Cw(z_pbl, RL, z0m, z0h), PSIh_y(pcrumin(z_pbl)/RL))
# external resistance
re_i = (zdh - C_i) / (k * RUstar)
H_i = rhoacp * T0ta / re_i
H_i = ifthenelse(H_i > H_d, H_d, H_i)
H_i = ifthenelse(H_i < H_w, H_w, H_i)
report(H_i, hmap)

# Calculate evaporation variables
print("Calculating relative evaporation and evaporative fraction...")
# Calculate relative evaporation
ev_r = 1.0 - (H_i - H_w) / (H_d - H_w) # set water and wet surfaces to 1.0
report(ev_r, evaprmap)
# Calculate evaporative fraction
Evapfr = ev_r * (1.0 - H_w / H_d)
report(Evapfr, evapfrmap)

# Calculate latent energy flux
print("Calculating Latent Energy Flux...")
labdaE = Evapfr * (Rn - G0)
labdaE2 = Rn - G0 - H_i
#assert(labdaE == labdaE2) # Check on closure of energy balance components!
report(labdaE, lemap)

# Calculate evapotranspiration flux
print("Calculating Evapotranspiration Flux...")
rhow = 998.0 # density of water [kg m-3]
E = labdaE / (L_e * rhow) #[m/s]
report(E, etmap)
#Ehour = E * 3600.0 * 100.0 # cm/h for data assimilation with sm model
#report(Ehour, "ehour_aster90.map")
#checkFile.close()
endtime = time()
deltaTime = endtime - starttime
print()
print("=============================================")
print("The model has been running for %5.2f seconds." % deltaTime)
print
print("Credits:")
print("Bob Su (ITC)")
print("Ambro Gieske (ITC)")
print("Wim Timmermans (ITC)")
print("Victor Jetten (UU)")
print("Steven de Jong (UU)")
print("Li Jia (WUR)")
print("Kor de Jong (UU)")
print("Derek Karssenberg (UU)")
print("=============================================")


