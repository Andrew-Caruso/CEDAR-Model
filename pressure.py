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
from numpy import isnan as isnan 
from numpy import nan as nan 
from numpy import reshape as reshape
from matplotlib.patches import Circle as circle 
from matplotlib.patches import Wedge as wedge 
from matplotlib.lines import Line2D as line


#global constants all in CGS 
#off set the aperture (via off set photons) 
#see photonGeneration fun 
x_off = 0.4 #cm 
y_off = 0.0 #cm 
#NOTE: the x_off set in the photons = - x_off set of the actual aperature
#if photons are shifted left the aperature essentially was shifted right 
P_select = 1.8E+6 #baryes printed later in the code right before kaon rings are generated
#smear the photons (slight deviation from the actual radius)
#see photonGeneration fun 
sigma = 0.05 #cm std, smearing of the photons  
#originally 0.1 cm for both offs and 0.01 cm std
print("x_off:",x_off,"cm")
print("y_off:",y_off,"cm")
print("sigma:",sigma,"cm") 


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

#valid ring meaning r_min < r < r_min, fits in the aperture 
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

#sectors 1 and 8 = up
#sectors 4 and 5 = down 
#sectors 2 and 3 = right
#sectors 6 and 7 = left 
def photonGeneration(r): 
	numG = 18 #number of photons to generate on the ring 
	chart = zeros(8) #index = sector-1, used as binary 0 or 1 
	#for a PM of the sector detecting 1 or more photons 
	#from sector 1 (index 0) to sector 8 (index 7)
	count = zeros(8) #used as number of photons in each sector 
	mu = r #mean for normal/gaussian distribution 
	left = 0
	right = 0
	up = 0
	down = 0
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
			count[1] += 1 #increase count of photons in sector 2
		elif (phi > pi/4 and phi < pi/2):
			#print("sector 1")
			if chart[0] != 1:
				chart[0] = 1	
			count[0] += 1 #increase count of photons in sector 1
		elif (phi > pi/2 and phi < 3*pi/4): 
			#print("sector 8")
			if chart[7] != 1:
				chart[7] = 1	
			count[7] += 1 #increase count of photons in sector 8
		elif phi > 3*pi/4 and phi < pi:
			#print("sector 7")
			if chart[6] != 1:
				chart[6] = 1	
			count[6] += 1 #increase count of photons in sector 7
		elif phi > pi and phi < 5*pi/4:
			#print("sector 6") 
			if chart[5] != 1:
				chart[5] = 1	
			count[5] += 1 #increase count of photons in sector 6
		elif phi > 5*pi/4 and phi < 3*pi/2:
			#print("sector 5") 
			if chart[4] != 1:
				chart[4] = 1	
			count[4] += 1 #increase count of photons in sector 5
		elif phi > 3*pi/2 and phi < 7*pi/4:
			#print("sector 4")
			if chart[3] != 1:
				chart[3] = 1	
			count[3] += 1 #increase count of photons in sector 4
		elif phi > 7*pi/4 and phi < 2*pi: 
			#print("sector 3") 
			if chart[2] != 1:
				chart[2] = 1	
			count[2] += 1 #increase count of photons in sector 3
		else:
			print("no sector") 
	#print("chart:",chart)
	totCount = sum(chart) #total number of sectors that made a detection (1 or more photons)
	#total count is the number of coincidences Nc, number of PMs that coincide with the same detection (recieve photons from same source/particle 
	up = count[0] + count[7] #sector 1 and 8
	down = count[3] + count[4] #sector 4 and 5
	right = count[1] + count[2] #sector 2 and 3
	left = count[5] + count[6] #sector 6 ad 7 
	#asymmetires
	return totCount,up,down,left,right 

#get the asymmetries in number of photons in up,down,left,right zones
#asymmetries max and min values are 1 and -1 
#technical asymmetry denotes asymmetry index
#where asymmetry index for up/down is:
#asymmetry index = count up - count down / (up count + down count)
#if up >> down then index = 1
#if down >> up then index = -1
def getAsymmetries(up,down,left,right):
	#only print out if the value is not nan 
	if left+right != 0:
		A_LR = (left-right) / (left+right)
	if up+down !=0: 
		A_UD = (up-down) / (up+down)
	else: 
		A_LR = 0
		A_UD = 0
	print("A_LR:",A_LR)
	print("A_UD:",A_UD,"\n")
	return 0 

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

#dimensions of ring-like diaphram and aperture 
#use the radius of the kaon ring, i.e. center the diaphram on the kaon ring 
#with deviations of 0.9 mm from the radius, thus a thickness of 1.8 mm
#with deviations of 0.09 cm from the radius, thus a thickness of 0.18 cm
#NOTE: use the kaon ring at 1.75 bar and adjust it by 0.9 mm to create the 
#aperture of the ring-like diaphram 
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
#empty list the number of photons per iteration for each zone (up,down,left,right) 
listup1 = []
listdown1 = []
listleft1 = []
listright1 = []
listup2 = []
listdown2 = []
listleft2 = []
listright2 = []
listup3 = []
listdown3 = []
listleft3 = []
listright3 = []

#select a single pressure instead of scanning all pressures, do:
#was 1.8E+6 or 1.75E+6 baryes
#see earliest lines  
print("P_select:",P_select,"Ba\n") 


#monte carlo simulation for kaon 
#use spectrum of momentum NOT random generated momentum, to decrease fluctations in the histograms 
#for each momentum (i.e. ring's radius)  generate 1E+2 PM detections of 18 photons  
#removed the other for loop and the random pressure in order to do a pressure scan instead
print("for kaon") 
#for i in range(0,N,1):
for P in range(int(P_min),int(P_max),2000):
	#P = rnd.uniform(P_min,P_max)
	if P != P_select:   
		continue 
	n = index(P)
	r = radius(mka,n) 
	print("selected pressure radius:",r) 
	stat = check(r,P,mka,n,r_min,r_max,P_min,P_max)
	stat5p1 = 0 #for Nc >= 5, stat 5+ 
	stat6p1 = 0 #for Nc >= 6
	stat7p1 = 0 #for Nc >= 7
	stat8p1 = 0 #for Nc = 7 or 7+ 
	up1 = 0 #number of photons in sectors 1 & 8
	down1 = 0 #number of photons in sectors 4 & 5
	left1 = 0 #number of photons in sectors 2 & 3 
	right1 =0 #number of photons in sectors 6 & 7
	sumup1 = 0 #total number of photons detected in up zone 
	sumdown1 = 0
	sumleft1 = 0 
	sumright1 = 0 
	totalcount1 = 0
	ratioup1 = 0
	ratiodown1 = 0
	ratioleft1 = 0
	ratioright1 = 0
	if stat == 1: 
		stat_array1 = append(stat_array1, P)
	#for each pressure, generate more rings
	#originally 1E+2 
	#generate rings
	for i in range(0,int(1E+4)): 
		#number of coincidences (see photonGen fun)  
		Nc1,up1,down1,left1,right1= photonGeneration(r) 		
		#add to respective lists 
		listup1.append(up1)
		listdown1.append(down1)		
		listleft1.append(left1) 
		listright1.append(right1)
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
#get the total sums of photons per zone 
sumup1 = sum(i for i in listup1 if i != nan) 
sumdown1 = sum(i for i in listdown1 if i != nan) 
sumleft1 = sum(i for i in listleft1 if i != nan) 
sumright1 = sum(i for i in listright1 if i != nan) 
totalcount1 = sumup1 + sumdown1 + sumleft1 + sumright1 
if totalcount1 == 0:
	totalcount1 = 1
	#to handel division by zero
ratioup1 = sumup1 / totalcount1 
ratiodown1 = sumdown1 / totalcount1 
ratioleft1 = sumleft1 / totalcount1 
ratioright1 = sumright1 / totalcount1 
print("up",sumup1) 
print("down",sumdown1)
print("left",sumleft1)
print("right",sumright1)
print("up ratio:",ratioup1)
print("down ratio:",ratiodown1)
print("left ratio:",ratioleft1)
print("right ratio:",ratioright1)
getAsymmetries(sumup1,sumdown1,sumleft1,sumright1) 
kaonring = circle((0,0),r-5,fill=False,color="tab:blue",label="kaon",alpha=1)



#monte carlo simulation for proton  
print("for proton") 
#for i in range(0,N,1):
for P in range(int(P_min),int(P_max),2000):
	#P = rnd.uniform(P_min,P_max)
	if P != P_select:   
		continue 
	n = index(P)
	r = radius(mpr,n) 
	print("selected pressure radius:",r)
	stat = check(r,P,mpr,n,r_min,r_max,P_min,P_max)
	stat5p2 = 0 #for Nc >= 5, stat 5+ 
	stat6p2 = 0 #for Nc >= 6
	stat7p2 = 0 #for Nc >= 7
	stat8p2 = 0 #for Nc = 8 
	up2 = 0 #number of photons in sectors 1 & 8
	down2 = 0 #number of photons in sectors 4 & 5
	left2 = 0 #number of photons in sectors 2 & 3 
	right2 =0 #number of photons in sectors 6 & 7
	sumup2 = 0 #total number of photons detected in up zone 
	sumdown2 = 0
	sumleft2 = 0 
	sumright2 = 0 
	totalcount2 = 0
	ratioup2 = 0
	ratiodown2 = 0
	ratioleft2 = 0
	ratioright2 = 0
	if stat == 1: 
		stat_array2 = append(stat_array2, P)
	for i in range(0,int(1E+2)):
		#number of coincidences (see photonGen fun)  
		Nc2,up2,down2,left2,right2= photonGeneration(r) 		
		#add to respective lists 
		listup2.append(up2)
		listdown2.append(down2)		
		listleft2.append(left2) 
		listright2.append(right2)
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
#get the total sums of photons per zone 
sumup2 = sum(i for i in listup2 if i != nan) 
sumdown2 = sum(i for i in listdown2 if i != nan) 
sumleft2 = sum(i for i in listleft2 if i != nan) 
sumright2 = sum(i for i in listright2 if i != nan) 
totalcount2 = sumup2 + sumdown2 + sumleft2 + sumright2 
if totalcount2 == 0:
	totalcount2 = 1
	#to handel division by zero
ratioup2 = sumup2 / totalcount2 
ratiodown2 = sumdown2 / totalcount2 
ratioleft2 = sumleft2 / totalcount2 
ratioright2 = sumright2 / totalcount2 
print("up",sumup2) 
print("down",sumdown2)
print("left",sumleft2)
print("right",sumright2)
print("up ratio:",ratioup2)
print("down ratio:",ratiodown2)
print("left ratio:",ratioleft2)
print("right ratio:",ratioright2)
getAsymmetries(sumup2,sumdown2,sumleft2,sumright2) 
protonring = circle((0,0),r-5,fill=False,color="tab:orange",label="proton",alpha=1)



#monte carlo simulation for pion  
print("for pion") 
#for i in range(0,N,1):
for P in range(int(P_min),int(P_max),2000):
	#P = rnd.uniform(P_min,P_max)
	if P != P_select:
		continue 
	n = index(P)
	r = radius(mpi,n) 
	print("selected pressure radius:",r)
	stat5p3 = 0 #for Nc >= 5, stat 5+ 
	stat6p3 = 0 #for Nc >= 6
	stat7p3 = 0 #for Nc >= 7
	stat8p3 = 0 #for Nc = 8 
	stat = check(r,P,mpi,n,r_min,r_max,P_min,P_max)
	up3 = 0 #number of photons in sectors 1 & 8
	down3 = 0 #number of photons in sectors 4 & 5
	left3 = 0 #number of photons in sectors 2 & 3 
	right3 =0 #number of photons in sectors 6 & 7
	sumup3 = 0 #total number of photons detected in up zone 
	sumdown3 = 0
	sumleft3 = 0 
	sumright3 = 0 
	totalcount3 = 0
	ratioup3 = 0
	ratiodown3 = 0
	ratioleft3 = 0
	ratioright3 = 0
	if stat == 1: 
		stat_array3 = append(stat_array3, P)
	for i in range(0,int(1E+2)):
		#number of coincidences (see photonGen fun)  
		Nc3,up3,down3,left3,right3= photonGeneration(r) 		
		#add to respective lists 
		listup3.append(up3)
		listdown3.append(down3)		
		listleft3.append(left3) 
		listright3.append(right3)
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
#get the total sums of photons per zone 
sumup3 = sum(i for i in listup3 if i != nan) 
sumdown3 = sum(i for i in listdown3 if i != nan) 
sumleft3 = sum(i for i in listleft3 if i != nan) 
sumright3 = sum(i for i in listright3 if i != nan) 
totalcount3 = sumup3 + sumdown3 + sumleft3 + sumright3 
if totalcount3 == 0:
	totalcount3 = 1
	#to handel division by zero
ratioup3 = sumup3 / totalcount3 
ratiodown3 = sumdown3 / totalcount3 
ratioleft3 = sumleft3 / totalcount3 
ratioright3 = sumright3 / totalcount3 
print("up",sumup3) 
print("down",sumdown3)
print("left",sumleft3)
print("right",sumright3)
print("up ratio:",ratioup3)
print("down ratio:",ratiodown3)
print("left ratio:",ratioleft3)
print("right ratio:",ratioright3)
getAsymmetries(sumup3,sumdown3,sumleft3,sumright3) 
pionring = circle((0,0),r-5,fill=False,color="tab:green",label="pion",alpha=1)


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
ax1.set_title("Pressure Scan of Particles")
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
ax3.set_title("x_off: %fcm y_off: %fcm std: %fcm" % (x_off,y_off,sigma))
ax3.legend(loc="best") 

#get midpoints and plot "errorbars" in figure 2 
#create figure 2
fig2,ax2 = plt.subplots(1,1) 
mid1 = (v1[:-1] + v1[1:]) / 2 
mid2 = (v2[:-1] + v2[1:]) / 2 
mid3 = (v3[:-1] + v3[1:]) / 2 
#print("mids1 shape:",shape(mid1))
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
#for figure 1
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
#3 peaks should roughly be in the ratio pion:proton:kaon of 70:24:6
print("expected pion:proton:kaon = 70:24:6") 
'''

#find max value for figure 3 
indexVec5p1 = where(b5p1 == b5p1.max())
indexVec6p1 = where(b6p1 == b6p1.max()) 
indexVec7p1 = where(b7p1 == b7p1.max()) 
indexVec8p1 = where(b8p1 == b8p1.max()) 

indexVec5p2 = where(b5p2 == b5p2.max())
indexVec6p2 = where(b6p2 == b6p2.max()) 
indexVec7p2 = where(b7p2 == b7p2.max()) 
indexVec8p2 = where(b8p2 == b8p2.max()) 

indexVec5p3 = where(b5p3 == b5p3.max())
indexVec6p3 = where(b6p3 == b6p3.max()) 
indexVec7p3 = where(b7p3 == b7p3.max()) 
indexVec8p3 = where(b8p3 == b8p3.max()) 

index5p1 = indexVec5p1[0] 
index6p1 = indexVec6p1[0]
index7p1 = indexVec7p1[0] 
index8p1 = indexVec8p1[0] 

index5p2 = indexVec5p2[0] 
index6p2 = indexVec6p2[0]
index7p2 = indexVec7p2[0] 
index8p2 = indexVec8p2[0] 

index5p3 = indexVec5p3[0] 
index6p3 = indexVec6p3[0]
index7p3 = indexVec7p3[0] 
index8p3 = indexVec8p3[0] 

vmax5p1 = ( ((v5p1[index5p1]+v5p1[index5p1+1]))/2 )
vmax6p1 = ( ((v6p1[index6p1]+v6p1[index6p1+1]))/2 )
vmax7p1 = ( ((v7p1[index7p1]+v7p1[index7p1+1]))/2 )
vmax8p1 = ( ((v8p1[index8p1]+v8p1[index8p1+1]))/2 )

vmax5p2 = ( ((v5p2[index5p2]+v5p2[index5p2+1]))/2 )
vmax6p2 = ( ((v6p2[index6p2]+v6p2[index6p2+1]))/2 )
vmax7p2 = ( ((v7p2[index7p2]+v7p2[index7p2+1]))/2 )
vmax8p2 = ( ((v8p2[index8p2]+v8p2[index8p2+1]))/2 )

vmax5p3 = ( ((v5p3[index5p3]+v5p3[index5p3+1]))/2 )
vmax6p3 = ( ((v6p3[index6p3]+v6p3[index6p3+1]))/2 )
vmax7p3 = ( ((v7p3[index7p3]+v7p3[index7p3+1]))/2 )
vmax8p3 = ( ((v8p3[index8p3]+v8p3[index8p3+1]))/2 )

#create wedges of the zones for rings graph  
#the ring-like (annulus) aperture is off set by negative offset x and y
#because the offset moves the photons' rings instead of moving the actual aperture position 
#so essentially the aperture is shifted by -off_x and -off_y 
maxring = circle((-x_off,-y_off),r_max-5,fill=False,color="black",label="ringlike aperture",alpha=1)
minring = circle((-x_off,-y_off),r_min-5,fill=False,color="black",alpha=1)

#NOTE: all radii for fig4 have 5 cm removed to make scaling better for visual
#start of fig4 (for single pressure, NOT for pressure scan or monte carlo) 
#create a figure for graphing the rings 
fig4,ax4 = plt.subplots(1,1) 
originPoint= circle((0,0),0.05,fill=True,color="black",alpha=1)
#plot borders to segregate the sectors 
length = int(r_max)-3
xpos = [-length,length]
ypos = [-length,length]
xneg = [length,-length]
yneg = [-length,length] 
xhor = [-length,length]
yhor = [0,0]
xver = [0,0]
yver = [-length,length]
ax4.plot(xpos,ypos,linestyle="dashed",color="black") 
ax4.plot(xneg,yneg,linestyle="dashed",color="black") 
ax4.plot(xhor,yhor,linestyle="dashed",color="black") 
ax4.plot(xver,yver,linestyle="dashed",color="black") 
#create wedges (sectors of left right up and down)
rightwedge = wedge((0,0),length/4,-45,45,width=None,color="red",label="right",alpha=0.2)
leftwedge = wedge((0,0),length/4,135,225,width=None,color="green",label="left",alpha=0.2)
upwedge = wedge((0,0),length/4,45,135,width=None,color="blue",label="up",alpha=0.2)
downwedge = wedge((0,0),length/4,225,315,width=None,color="orange",label="down",alpha=0.2)
#label each sector
#by default text is written from the upper right of each text position point 
#horizontal alignment (ha) = center, right, left
#center on the point, right of the point, or left of the point, horizontally  
#vertical alignment = center, top, bottom, baseline, center_baseline
#center on the point, top of point, bottom of point vertically 
#baseline means the bottom of text is on where the center point is vertically
#center)baseline means slightly lower than center point 
#ax4.set_anchor("N") 
ax4.text(length/3,length/length,"sector 2",va="center",ha="left",rotation=0)
ax4.text(length/3,-length/length,"sector 3",va="center",ha="left",rotation=0)
ax4.text(-length/3,length/length,"sector 7",va="center",ha="right",rotation=0)
ax4.text(-length/3,-length/length,"sector 6",va="center",ha="right",rotation=0)
ax4.text(length/4,length-length*0.55,"sector 1",va="center",ha="right",rotation=0)
ax4.text(-length/4,length-length*0.55,"sector 8",va="center",ha="left",rotation=0)
ax4.text(-length/4,-length+length*0.55,"sector 5",va="center",ha="left",rotation=0)
ax4.text(length/4,-length+length*0.55,"sector 4",va="center",ha="right",rotation=0)
#add all the patches (or shape, i.e. polygon, circle, wedge, etc) 
ax4.add_patch(originPoint)
ax4.add_patch(kaonring)
ax4.add_patch(protonring)
ax4.add_patch(pionring)
ax4.add_patch(rightwedge)
ax4.add_patch(leftwedge)
ax4.add_patch(upwedge)
ax4.add_patch(downwedge)
ax4.add_patch(maxring)
ax4.add_patch(minring)
#edit the axis
ax4.set_xscale("linear") 
ax4.set_yscale("linear")
ax4.set_ylabel("y in cm")
ax4.set_xlabel("x in cm")
ax4.set_ylim(ymin=-length,ymax=length)
ax4.set_xlim(xmin=-length,xmax=length)
ax4.legend(loc="best", prop={'size':7}) 
#end of fig4

#end the timer
print("%s seconds" % (time.time()-start_time))
#display the 4 figures
plt.show()
