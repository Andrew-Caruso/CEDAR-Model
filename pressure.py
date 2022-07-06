#task 3 monte carlo simulation of successful rings for different pressures assuming the simplest cone model 
import random as rnd 
import numpy as np
import matplotlib.pyplot as plt 
import collections as col 

#global constants all in CGS 
radiator = "N2" 
print("radiator:",radiator)
#n depends on pressure, temperature, and wavelength 
#all naught variables are together
#all naught variables are FIXED
n0 = 1.00029910 #better value
#from website https://refractiveindex.info/?shelf=main&book=N2&page=Borzsonyi
#assuming 0.42 micrometer wavelength of cherenkov rad
T0 = 273.15 #K 0 deg celsius 
P0 = 1E+6 #Ba 1 bar 
T = 293 #K 
#let T be fixed for all particles and pressure (see index fun) 

c = 29979245800 #cm/s 
hbar = 1.05457266E-27 #g cm^2 / s 
#1 eV = 1.60217733E-12 g cm^2 / s^2
factor = 1.60217733E-12 #erg / eV 
mpi = 139.57039E+6 #eV 140MeV/c mass of +pi
mka = 493.677E+6 #eV 494MeV/c mass of +kaon
mpr = 938.27208816E+6 #eV 938MeV/c mass of proton 
mel = 0.51099895000E+6 #eV 0.51 MeV/c mass of electron 
#convert masses from eV/c to g 
mpi = mpi*factor/c**2 
mka = mka*factor/c**2 
mpr = mpr*factor/c**2
mel = mel*factor/c**2 
l = 300 #cm or 3m distance of cone 
pi = np.pi 
theta_min = 0
theta_max = pi*2
P_max = 2E+6 #Ba 2 bar 
P_min = 1.6E+6 #Ba 1.6 bar
print("P_min: ",P_min,"Ba")
print("P_max: ",P_max,"Ba")
#true P min and max: vacuum and 5 bar respectively 
P_step = 1
p = (75*(1E+9))*factor/c #75 GeV/c to g cm/s momentum  
#momentum is FIXED 

'''
DISCARD: 
R_min = 0.003 #cm 0.03 mm, minimum inner radius of ring (aperture lower bound)
R_max = 2 #cm 20 mm, maximum inner radius of ring (aperture upper bound)
D_max = 10 #cm 100 m #maximum outer radius of the ring (diaphram radius)
#because the aperture size can be changed then 
#radius of aperture is selected such that R_min < r < R_max 
'''
#notes on pressure P 
#1 barye (Ba) = 1 g /(cm s^2) = 0.1 Pa  
#1 Pa = 10 Ba 
#1 bar = 1E+5 Pa 
#1 bar = 1E+6 Ba 
#pressure goes from vacuum to 5 bar
#p vs P 
#keep momentum p constant at 75 GeV 

#radius function 
def radius(m,n):
	beta = 1 / (np.sqrt(((m*c/p)**2)+1))
	#print(beta*n)
	if (beta*n < 1):
		stat = 0
		r = 0
	else:
		theta = np.arccos(1/(n*beta))
		r = l*np.tan(theta)
	#print(r)
	return r 

#check the ring's radius 
def check(r,P,m,n,r_min,r_max,P_min,P_max):
	boolean = 0
	beta = 1 / (np.sqrt(((m*c/p)**2)+1))
	if ((beta*n >= 1) and (r_min < r) and (r_max > r) and (P_min < P) and (P_max > P)): 
		boolean = 1
	else: 
		boolean = 0
	return boolean 

def findPressure(r,m):
	beta = 1 / (np.sqrt(((m*c/p)**2)+1))
	Pout = (1/(beta*np.cos(np.arctan(r/l))) -1 )* P0 / (n0 -1) 
	return Pout 

#generate 18 photons on a valid ring using 2D and polar coor 
#valid ring meaning r_min < r < r_min, fits in the aperature 
#then determine which sector each photon is in and output 
#wherein cherenkov angle (theta) = polar angle off of the beam axis
#phi = azimuthal angle off of x axis, in the transverse plane in respect to beam 

def photonGeneration(r): 
	numG = 1 #number of photons to generate on the ring 
	for i in range(0,numG,1): 
		phi = rnd.uniform(0,2*pi)
		x = r*np.cos(phi)
		y = r*np.sin(phi)
		print("hypotenuse:",np.sqrt((x**2)+(y**2)))
		print("r:",r) 
		print("x:",x)
		print("y:",y) 
		print("phi:",phi) 
		if (phi > 0 and phi < pi/4): 
			print("sector 2") 
		elif (phi > pi/4 and phi < pi/2):
			print("sector 1")
		elif (phi > pi/2 and phi < 3*pi/4): 
			print("sector 8")
		elif phi > 3*pi/4 and phi < pi:
			print("sector 7")
		elif phi > pi and phi < 5*pi/4:
			print("sector 6") 
		elif phi > 5*pi/4 and phi < 3*pi/2:
			print("sector 5") 
		elif phi > 3*pi/2 and phi < 7*pi/4:
			print("sector 4")
		elif phi > 7*pi/4 and phi < 2*pi: 
			print("sector 3") 
		else:
			print("no sector") 

#index of refraction 
def index(P):
	#n = (n0-1)*(P/P0) + 1 
	n = (n0-1)*(T0/P0)*(P/T)+1
	return n

#dimensions of ring-like diaphram and aperture 
#use the radius of the kaon ring, i.e. center the diaphram on the kaon ring 
#with deviations of 0.9 mm from the radius, thus a thickness of 1.8 mm
#with deviations of 0.09 cm from the radius, thus a thickness of 0.18 cm
#NOTE: use the kaon ring at 1.75 bar and adjust it by 0.9 mm to create the 
#aperature of the ring-like diaphram 

#N2 gas at 1.75 bar, find n and r
P_N2 = 1.75E+6 #Ba or 1.75 bar
n = index(P_N2)  
r_kaon = radius(mka,n)
print("r_kaon:",r_kaon,"cm")
dev = 0.09 #cm deviation from radius 
r_min = r_kaon - dev #minimum radius of the aperture 
r_max= r_kaon + dev #maximum radius of aperture 
thick = 0.18 #cm thickness of diaphram ring  
#print(r_max-r_min) #precision issue?? 
print("r_min: ", r_min,"cm")
print("r_max: ", r_max,"cm")

#all radii of rings that are not within r_min and r_max are blocked and out of focus, do not pass the diaphram's aperture, are denoted by 0
#all radii of rings that are within r_min and r_max are in focus on the diaphram, thus pass and are detected by PMs, denoted by 1  
#vary the pressure with temperature, momentum, and mass fixed 
#change pressure -> change refractive index n  -> change theta -> change r   

#set defaults for variables 
n = 1 #index of refraction will be changed by changing pressure via index fun 
P = 1E+6 #Ba 1 bar, as default, but will be changed 
r = 0 #radius of ring
stat = 0 #boolean that is either 1 or 0 depending if the ring is in the boundary or not respectively
#number of radii to generate  
#100 points is good enough
N = int(1E+3) #CHANGE ME

#array of P list
#2 rows N columns 
stat_array1 = np.array([]) #empty 
stat_array2 = np.array([]) #empty 
stat_array3 = np.array([]) #empty 

#monte carlo simulation for kaon 
print("for kaon") 
for i in range(0,N,1):
	P = rnd.uniform(P_min,P_max)
	n = index(P)
	r = radius(mka,n) 
	stat = check(r,P,mka,n,r_min,r_max,P_min,P_max)
	if stat == 1: 
		stat_array1 = np.append(stat_array1, P)
		photonGeneration(r) 
		'''
		print("iteration: ",i)
		print("pressure: ",P,"Ba")
		print("refractive index: ",n) 
		print("radius: ",r,"cm")
		print("stat: ",stat,"\n")
		'''
#monte carlo simulation for proton  
print("for proton") 
for i in range(0,N,1):
	P = rnd.uniform(P_min,P_max)
	n = index(P)
	r = radius(mpr,n) 
	stat = check(r,P,mpr,n,r_min,r_max,P_min,P_max)
	if stat == 1: 
		stat_array2 = np.append(stat_array2, P)
		photonGeneration(r) 
		'''
		print("iteration: ",i)
		print("pressure: ",P,"Ba")
		print("refractive index: ",n) 
		print("radius: ",r,"cm")
		print("stat: ",stat,"\n")
		'''
#monte carlo simulation for pion  
print("for pion:") 
for i in range(0,N,1):
	P = rnd.uniform(P_min,P_max)
	n = index(P)
	r = radius(mpi,n) 
	stat = check(r,P,mpi,n,r_min,r_max,P_min,P_max)
	if stat == 1: 
		stat_array3 = np.append(stat_array3, P)
		photonGeneration(r) 
		'''
		print("iteration: ",i)
		print("pressure: ",P,"Ba")
		print("refractive index: ",n) 
		print("radius: ",r,"cm")
		print("stat: ",stat,"\n")
		'''

#create the histogram
#print(np.shape(stat_array1))
#number of bins, usually 200 is good 
num = 200 #CHANGE ME  
y1,x1, _ = plt.hist(stat_array1,bins=num,range=[P_min,P_max],label="kaon")
y2,x2, _ = plt.hist(stat_array2,bins=num,range=[P_min,P_max],label="proton")
y3,x3, _ = plt.hist(stat_array3,bins=num,range=[P_min,P_max],label="pion")
#y are bins, x are values 
#find most common value and the count 
'''
print("y1:",y1)
print("x1:",x1)
print("y2:",y2)
print("x2:",x2)
print("y3:",y3)
print("x3:",x3)
'''
indexarray1 = np.where(y1 == y1.max())  
indexarray2 = np.where(y2 == y2.max()) 
indexarray3 = np.where(y3 == y3.max()) 
index1 = indexarray1[0] 
index2 = indexarray2[0] 
index3 = indexarray3[0] 
'''
print(type(index1)) 
print("kaon greatest frequency:",y1.max(),"at index:",index1,"pressure:",x1[index1],"Ba")
print("proton greatest frequency:",y2.max(),"at index:",index2,"pressure:",x2[index2],"Ba")
print("pion greatest frequency:",y3.max(),"at index:",index3,"pressure:",x3[index3],"Ba")
'''
plt.xlabel("pressure (Ba)")
plt.ylabel("frequency") 
plt.legend(loc="best") 
plt.show()

