#task 3 monte carlo simulation of successful rings for different pressures assuming the simplest cone model 
import random as rnd 
import time
import matplotlib.pyplot as plt 
from numpy import sqrt as sqrt 
from numpy import array as array 
from numpy import append as append 
from numpy import zeros as zeros 
from numpy import cos as cos 
from numpy import sin as sin 
from numpy import arctan as arctan 
from numpy import arccos as arccos
from numpy import tan as tan 
from numpy import shape as shape 
from numpy import where as where 
from numpy import pi as pi 
from numpy import random as rand 

#global constants all in CGS 
start_time = time.time() #start the timer
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
#notes on pressure P 
#1 barye (Ba) = 1 g /(cm s^2) = 0.1 Pa  
#1 Pa = 10 Ba 
#1 bar = 1E+5 Pa 
#1 bar = 1E+6 Ba 
#pressure goes from vacuum to 5 bar
#p vs P 
#keep momentum p constant at 75 GeV 
#off set the aperture (via off set photons) 
#see photonGeneration fun 
x_off = 0.1 #cm 
y_off = 0.1 #cm 

#create all the functions
#radius function 
def radius(m,n):
	beta = 1 / (sqrt(((m*c/p)**2)+1))
	#print(beta*n)
	if (beta*n < 1):
		stat = 0
		r = 0
	else:
		theta = arccos(1/(n*beta))
		r = l*tan(theta)
	#print(r)
	return r 

#check the ring's radius 
def check(r,P,m,n,r_min,r_max,P_min,P_max):
	boolean = 0
	beta = 1 / (sqrt(((m*c/p)**2)+1))
	if ((beta*n >= 1) and (r_min < r) and (r_max > r) and (P_min < P) and (P_max > P)): 
		boolean = 1
	else: 
		boolean = 0
	return boolean 

def findPressure(r,m):
	beta = 1 / (sqrt(((m*c/p)**2)+1))
	Pout = (1/(beta*cos(arctan(r/l))) -1 )* P0 / (n0 -1) 
	return Pout 

#generate 18 photons on a valid ring using 2D and polar coor 
#valid ring meaning r_min < r < r_min, fits in the aperature 
#then determine which sector each photon is in and output 
#wherein cherenkov angle (theta) = polar angle off of the beam axis
#phi = azimuthal angle off of x axis, in the transverse plane in respect to beam 
#need coincidences to be 5,6,7 or 8 to graph as valid detection by photomultipliers PMs
#8 sectors = 8 PMs 
#presumption: need 1 or more photons to count as a detection by 1 sector (i.e. PM) 

#photon smearing: instead of having each photon exactly on the ring, given them their own radii that is not exactly on the radius r 
#this is done by chosing a radius from a normal/gaussian distribution centered on r, for each photon 

#shifting the aperture in the tranverse xy plane. 
#Instead of adjusting the actual poisition of the aperture, adjust the photon x and y positions, thus the radius of each photon 

def photonGeneration(r): 
	numG = 18 #number of photons to generate on the ring 
	chart = zeros(8) #index = sector-1  
	mu = r #mean for normal/gaussian distribution 
	sigma = 0.01 #cm std, smearing of the photons  
	#loop through each photon 
	for i in range(0,numG,1): 
		phi = rnd.uniform(0,2*pi) #generate random angle
		r0 = rand.normal(mu,sigma) #generate photon's radius from center
		x = r0*cos(phi) + x_off 
		y = r0*sin(phi) + y_off 
		r1 = sqrt(x**2 + y**2) 
		#if the photon radius is no longer in the ring, throw the iteration out and move on to the next photon  
		if r1 < r_min or r1 > r_max:
			continue
		'''
		print("hypotenuse:",sqrt((x**2)+(y**2)))
		print("r:",r) 
		print("x:",x)
		print("y:",y) 
		print("phi:",phi) 
		'''
		if (phi > 0 and phi < pi/4): 
			#print("sector 2") 
			if chart[1] != 1:
				chart[1] = 1 
		elif (phi > pi/4 and phi < pi/2):
			#print("sector 1")
			if chart[0] != 1:
				chart[0] = 1	
		elif (phi > pi/2 and phi < 3*pi/4): 
			#print("sector 8")
			if chart[7] != 1:
				chart[7] = 1	
		elif phi > 3*pi/4 and phi < pi:
			#print("sector 7")
			if chart[6] != 1:
				chart[6] = 1	
		elif phi > pi and phi < 5*pi/4:
			#print("sector 6") 
			if chart[5] != 1:
				chart[5] = 1	
		elif phi > 5*pi/4 and phi < 3*pi/2:
			#print("sector 5") 
			if chart[4] != 1:
				chart[4] = 1	
		elif phi > 3*pi/2 and phi < 7*pi/4:
			#print("sector 4")
			if chart[3] != 1:
				chart[3] = 1	
		elif phi > 7*pi/4 and phi < 2*pi: 
			#print("sector 3") 
			if chart[2] != 1:
				chart[2] = 1	
		else:
			print("no sector") 
	#print("chart:",chart)
	totCount = sum(chart) #total number of sectors that made a detection (1 or more photons)
	#total count is the number of coincidences Nc, number of PMs that coincide with the same detection (recieve photons from same source/particle 
	return totCount 

#get color based on the Nc
def getColor(Nc):
	color = "white"
	if Nc == 5:
		color = "magenta"
	elif Nc == 6:
		color = "darkred"
	elif Nc == 7:
		color = "lime"
	elif Nc == 8:
		color = "cyan"
	return color 

#index of refraction 
def index(P):
	#n = (n0-1)*(P/P0) + 1 
	n = (n0-1)*(T0/P0)*(P/T)+1
	return n

#generate binaries for 4 diff hists 
#for number of coincidences being 5+ (i.e. >=5), 6+, 7+, 8
def valuesHist5p(Nc):
	boolean = 0 
	if Nc >= 5:
		boolean = 1
	else:
		boolean = 0
	return boolean 

def valuesHist6p(Nc):
	boolean = 0 
	if Nc >= 6:
		boolean = 1
	else:
		boolean = 0
	return boolean 
		 
def valuesHist7p(Nc):
	boolean = 0 
	if Nc >= 7:
		boolean = 1
	else:
		boolean = 0
	return boolean 

def valuesHist8p(Nc):
	boolean = 0 
	if Nc == 8:
		boolean = 1
	else:
		boolean = 0
	return boolean 
'''
#filter out the zeros in the bins array b and corresponding values array v 
#get the midpoint values from an array 
def filterZerosGetMid(b,v):
	output = array([[],[]])
	#output[0] reps b output[1] reps v  
	indexb = where(b != 0) #get all indices where b is not zero 
	tempv = v[indexb] #temporary values  
	print("shape of output[0]",shape(output[0]))
	print("shape of b[indexb]",shape(b[indexb]))
	print("shape of output[1]",shape(output[1]))
	print("shape of tempv",shape(tempv))
	for i in range(0,len(b[indexb])):
		output[0] = append(output[0], b[indexb][i]) #new b values omitting the zeros  
	#get midpoints in tempv 
	#for i in range(0,len(tempv)):
		#output[1] = append(output[1], ((tempv[:-1] + tempv[1:]) / 2) )
	return output  
'''

#dimensions of ring-like diaphram and aperture 
#use the radius of the kaon ring, i.e. center the diaphram on the kaon ring 
#with deviations of 0.9 mm from the radius, thus a thickness of 1.8 mm
#with deviations of 0.09 cm from the radius, thus a thickness of 0.18 cm
#NOTE: use the kaon ring at 1.75 bar and adjust it by 0.9 mm to create the 
#aperature of the ring-like diaphram 
#N2 gas at 1.75 bar, find n and r
P_N2 = 1.75E+6 #Ba or 1.75 bar of N2 gas 
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

#create empty arrays for the histograms 
#array of P list
stat_array1 = array([]) #empty 
stat_array2 = array([]) #empty 
stat_array3 = array([]) #empty 
#array of Nc 5+, 6+, 7+, 8 
stat5p_array1 = array([]) 
stat6p_array1 = array([]) 
stat7p_array1 = array([]) 
stat8p_array1 = array([]) 
stat5p_array2 = array([]) 
stat6p_array2 = array([]) 
stat7p_array2 = array([]) 
stat8p_array2 = array([]) 
stat5p_array3 = array([]) 
stat6p_array3 = array([]) 
stat7p_array3 = array([]) 
stat8p_array3 = array([]) 

#monte carlo simulation for kaon 
#use spectrum of momentum NOT random generated momentum, to decrease fluctations in the histograms 
#for each momentum (i.e. ring's radius)  generate 1E+2 PM detections of 18 photons  
print("for kaon") 
#for i in range(0,N,1):
for P in range(int(P_min),int(P_max),2000):
	#P = rnd.uniform(P_min,P_max)
	n = index(P)
	r = radius(mka,n) 
	stat = check(r,P,mka,n,r_min,r_max,P_min,P_max)
	stat5p1 = 0 #for Nc >= 5, stat 5+ 
	stat6p1 = 0 #for Nc >= 6
	stat7p1 = 0 #for Nc >= 7
	stat8p1 = 0 #for Nc = 7 or 7+ 
	if stat == 1: 
		stat_array1 = append(stat_array1, P)
	for i in range(0,int(1E+2)): 
		Nc1 = photonGeneration(r) #number of coincidences (see photonGen fun)  
		#print("Nc1:",Nc1)
		color1 = getColor(Nc1)
		stat5p1 = valuesHist5p(Nc1) 
		stat6p1 = valuesHist6p(Nc1) 
		stat7p1 = valuesHist7p(Nc1) 
		stat8p1 = valuesHist8p(Nc1) 
		if stat5p1 == 1:
			stat5p_array1 = append(stat5p_array1, P)
		if stat6p1 == 1:
			stat6p_array1 = append(stat6p_array1, P)
		if stat7p1 == 1:
			stat7p_array1 = append(stat7p_array1, P)
		if stat8p1 == 1:
			stat8p_array1 = append(stat8p_array1, P)
		'''
		print("iteration: ",i)
		print("pressure: ",P,"Ba")
		print("refractive index: ",n) 
		print("radius: ",r,"cm")
		print("Nc: ",Nc1)
		print("stat: ",stat,"\n")
		'''

#monte carlo simulation for proton  
print("for proton") 
#for i in range(0,N,1):
for P in range(int(P_min),int(P_max),2000):
	#P = rnd.uniform(P_min,P_max)
	n = index(P)
	r = radius(mpr,n) 
	stat = check(r,P,mpr,n,r_min,r_max,P_min,P_max)
	stat5p2 = 0 #for Nc >= 5, stat 5+ 
	stat6p2 = 0 #for Nc >= 6
	stat7p2 = 0 #for Nc >= 7
	stat8p2 = 0 #for Nc = 8 
	if stat == 1: 
		stat_array2 = append(stat_array2, P)
		for i in range(0,int(1E+2)):
			Nc2 = photonGeneration(r) #number of coincidences (see photonGen fun)  
			stat5p2 = valuesHist5p(Nc2) 
			stat6p2 = valuesHist6p(Nc2) 
			stat7p2 = valuesHist7p(Nc2) 
			stat8p2 = valuesHist8p(Nc2) 
			if stat5p2 == 1:
				stat5p_array2 = append(stat5p_array2, P)
			if stat6p2 == 1:
				stat6p_array2 = append(stat6p_array2, P)
			if stat7p2 == 1:
				stat7p_array2 = append(stat7p_array2, P)
			if stat8p2 == 1:
				stat8p_array2 = append(stat8p_array2, P)
			'''
			print("iteration: ",i)
			print("pressure: ",P,"Ba")
			print("refractive index: ",n) 
			print("radius: ",r,"cm")
			print("Nc: ",Nc2)
			print("stat: ",stat,"\n")
			'''

#monte carlo simulation for pion  
print("for pion") 
#for i in range(0,N,1):
for P in range(int(P_min),int(P_max),2000):
	#P = rnd.uniform(P_min,P_max)
	n = index(P)
	r = radius(mpi,n) 
	stat5p3 = 0 #for Nc >= 5, stat 5+ 
	stat6p3 = 0 #for Nc >= 6
	stat7p3 = 0 #for Nc >= 7
	stat8p3 = 0 #for Nc = 8 
	stat = check(r,P,mpi,n,r_min,r_max,P_min,P_max)
	if stat == 1: 
		stat_array3 = append(stat_array3, P)
		for i in range(0,int(1E+2)):
			Nc3 = photonGeneration(r) #number of coincidences (see photonGen fun)  
			stat5p3 = valuesHist5p(Nc3) 
			stat6p3 = valuesHist6p(Nc3) 
			stat7p3 = valuesHist7p(Nc3) 
			stat8p3 = valuesHist8p(Nc3) 
			if stat5p3 == 1:
				stat5p_array3 = append(stat5p_array3, P)
			if stat6p3 == 1:
				stat6p_array3 = append(stat6p_array3, P)
			if stat7p3 == 1:
				stat7p_array3 = append(stat7p_array3, P)
			if stat8p3 == 1:
				stat8p_array3 = append(stat8p_array3, P)
			'''
			print("iteration: ",i)
			print("pressure: ",P,"Ba")
			print("refractive index: ",n) 
			print("radius: ",r,"cm")
			print("Nc: ",Nc3)
			print("stat: ",stat,"\n")
			'''

#create histograms for each particle of
#valid pressures only (whose generated radius falls within the ring of the diaphram)  
#number of bins, usually 200 is good 
num = 200 #CHANGE ME  
#create figure/canvas and axes/pads 
fig1, ax1 = plt.subplots(1,1)
'''
print("stat_array1:",stat_array1)
print("stat_array2:",stat_array2)
print("stat_array3:",stat_array3)
'''
print("shape of stat_array1:",shape(stat_array1))
#alternatively equivalent: fig, (ax1, ax2)
b1,v1, _ = ax1.hist(stat_array1,bins=num,range=[P_min,P_max],label="kaon",color='tab:blue')
b2,v2, _ = ax1.hist(stat_array2,bins=num,range=[P_min,P_max],label="proton",color='tab:orange')
b3,v3, _ = ax1.hist(stat_array3,bins=num,range=[P_min,P_max],label="pion",color='tab:green')
#b are bins, v are values of momentum marking the start and end of each bin, and _ are patches (ignored) 
ax1.set_xlabel("pressure (Ba)")
ax1.set_ylabel("frequency")
ax1.legend(loc="best") 
print("b1 shape:",shape(b1))
print("v1 shape:",shape(v1))
'''
print("b1",b1)
print("v1",v1)
'''

#create histograms for Nc 5+,6+,7+,8+ for each particle 
fig3,ax3 = plt.subplots(1,1) 
b5p1,v5p1, _ = ax3.hist(stat5p_array1,bins=num,range=[P_min,P_max],label="5+",color='magenta')
b6p1,v6p1, _ = ax3.hist(stat6p_array1,bins=num,range=[P_min,P_max],label="6+",color='darkred')
b7p1,v7p1, _ = ax3.hist(stat7p_array1,bins=num,range=[P_min,P_max],label="7+",color='lime')
b8p1,v8p1, _ = ax3.hist(stat8p_array1,bins=num,range=[P_min,P_max],label="8",color='cyan')

b5p2,v5p2, _ = ax3.hist(stat5p_array2,bins=num,range=[P_min,P_max],color='magenta')
b6p2,v6p2, _ = ax3.hist(stat6p_array2,bins=num,range=[P_min,P_max],color='darkred')
b7p2,v7p2, _ = ax3.hist(stat7p_array2,bins=num,range=[P_min,P_max],color='lime')
b8p2,v8p2, _ = ax3.hist(stat8p_array2,bins=num,range=[P_min,P_max],color='cyan')

b5p3,v5p3, _ = ax3.hist(stat5p_array3,bins=num,range=[P_min,P_max],color='magenta')
b6p3,v6p3, _ = ax3.hist(stat6p_array3,bins=num,range=[P_min,P_max],color='darkred')
b7p3,v7p3, _ = ax3.hist(stat7p_array3,bins=num,range=[P_min,P_max],color='lime')
b8p3,v8p3, _ = ax3.hist(stat8p_array3,bins=num,range=[P_min,P_max],color='cyan')
ax3.set_xlabel("pressure (Ba)")
ax3.set_ylabel("frequency")
ax3.legend(loc="best") 

#get midpoints and plot "errorbars" in figure 2 
#create figure 2
fig2,ax2 = plt.subplots(1,1) 
mid1 = (v1[:-1] + v1[1:]) / 2 
mid2 = (v2[:-1] + v2[1:]) / 2 
mid3 = (v3[:-1] + v3[1:]) / 2 
print("mids1 shape:",shape(mid1))
ax2.errorbar(x=mid1,y=b1,yerr=0,fmt='o',capsize=2)
ax2.errorbar(x=mid2,y=b2,yerr=0,fmt='o',capsize=2)
ax2.errorbar(x=mid3,y=b3,yerr=0,fmt='o',capsize=2)

#get midpoints for Nc 5+,6+,7+,8+ 
mid5p1 = (v5p1[:-1] + v5p1[1:]) / 2 
mid6p1 = (v6p1[:-1] + v6p1[1:]) / 2 
mid7p1 = (v7p1[:-1] + v7p1[1:]) / 2 
mid8p1 = (v8p1[:-1] + v8p1[1:]) / 2 
ax2.errorbar(x=mid5p1,y=b5p1,yerr=0,fmt='o',capsize=2,color="magenta")
ax2.errorbar(x=mid6p1,y=b6p1,yerr=0,fmt='o',capsize=2,color="darkred")
ax2.errorbar(x=mid7p1,y=b7p1,yerr=0,fmt='o',capsize=2,color="lime")
ax2.errorbar(x=mid8p1,y=b8p1,yerr=0,fmt='o',capsize=2,color="cyan")

mid5p2 = (v5p2[:-1] + v5p2[1:]) / 2 
mid6p2 = (v6p2[:-1] + v6p2[1:]) / 2 
mid7p2 = (v7p2[:-1] + v7p2[1:]) / 2 
mid8p2 = (v8p2[:-1] + v8p2[1:]) / 2 
ax2.errorbar(x=mid5p2,y=b5p2,yerr=0,fmt='o',capsize=2,color="magenta")
ax2.errorbar(x=mid6p2,y=b6p2,yerr=0,fmt='o',capsize=2,color="darkred")
ax2.errorbar(x=mid7p2,y=b7p2,yerr=0,fmt='o',capsize=2,color="lime")
ax2.errorbar(x=mid8p2,y=b8p2,yerr=0,fmt='o',capsize=2,color="cyan")

mid5p3 = (v5p3[:-1] + v5p3[1:]) / 2 
mid6p3 = (v6p3[:-1] + v6p3[1:]) / 2 
mid7p3 = (v7p3[:-1] + v7p3[1:]) / 2 
mid8p3 = (v8p3[:-1] + v8p3[1:]) / 2 
ax2.errorbar(x=mid5p3,y=b5p3,yerr=0,fmt='o',capsize=2,color="magenta")
ax2.errorbar(x=mid6p3,y=b6p3,yerr=0,fmt='o',capsize=2,color="darkred")
ax2.errorbar(x=mid7p3,y=b7p3,yerr=0,fmt='o',capsize=2,color="lime")
ax2.errorbar(x=mid8p3,y=b8p3,yerr=0,fmt='o',capsize=2,color="cyan")

#find most common value and the count 
'''
print("b1:",b1)
print("v1:",v1)
print("b2:",b2)
print("v2:",v2)
print("b3:",b3)
print("v3:",v3)
'''
indexarray1 = where(b1 == b1.max())  
indexarray2 = where(b2 == b2.max()) 
indexarray3 = where(b3 == b3.max()) 
index1 = indexarray1[0] 
index2 = indexarray2[0] 
index3 = indexarray3[0] 
#print(type(index1)) 
v1max = ((v1[index1]+v1[index1+1])/2)
v2max = ((v2[index2]+v2[index2+1])/2)
v3max = ((v3[index3]+v3[index3+1])/2)
'''
print("kaon greatest frequency:",b1.max(),"at index:",index1,"pressure:",v1max,"Ba")
print("proton greatest frequency:",b2.max(),"at index:",index2,"pressure:",v2max,"Ba")
print("pion greatest frequency:",b3.max(),"at index:",index3,"pressure:",v3max,"Ba") 
'''
#3 peaks should roughly be in the ratio pion:proton:kaon of 70:24:6
print("expected pion:proton:kaon = 70:24:6") 

#end the timer
print("%s seconds" % (time.time()-start_time))
#display the 3 figures
plt.show()
