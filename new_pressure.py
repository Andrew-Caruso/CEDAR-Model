#import necessary modules 
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
from numpy import abs as abs 
from matplotlib.patches import Circle as circle 
from matplotlib.patches import Wedge as wedge 
from matplotlib.lines import Line2D as line
#layout of code in sections:
#declaration of functions
#initalization of parameters
#computation 



	#declare all functions

#create radius 
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

#determine the pressure 
def findPressure(r,m):
	beta = 1 / (sqrt(((m*c/p)**2)+1))
	Pout = (1/(beta*cos(arctan(r/l))) -1 )* P0 / (n0 -1) 
	return Pout 

#generate photons 
def photonGeneration(r,m): 
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
	totCount = sum(chart) #total number of sectors that made a detection (1 or more photons)
	#total count is the number of coincidences Nc, number of PMs that coincide with the same detection (of photons from same source/particle 
	up = count[0] + count[7] #sector 1 and 8
	down = count[3] + count[4] #sector 4 and 5
	right = count[1] + count[2] #sector 2 and 3
	left = count[5] + count[6] #sector 6 ad 7 
	return totCount,up,down,left,right 

#calculate the asymmetries 
def getAsymmetries(up,down,left,right,totalcount):
	#only print out if the value is not nan 
	if left+right !=0:
		A_LR = (left-right) / (left+right)
	if up+down !=0: 
		A_UD = (up-down) / (up+down)
	else: 
		A_LR = 0
		A_UD = 0
	print("A_LR:",A_LR)
	print("A_UD:",A_UD)
	#to handle division by zero 
	if totalcount == 0:
		totalcount = 1
	ratioup = up / totalcount
	ratiodown = down / totalcount
	ratioleft = left / totalcount
	ratioright = right / totalcount
	print("up ratio:",ratioup)
	print("down ratio:",ratiodown)
	print("left ratio:",ratioleft)
	print("right ratio:",ratioright,"\n")
	return 0 

#calculate the refractive index 
def index(P):
	n = (n0-1)*(T0/P0)*(P/T)+1
	return n

#generate binaries for 4 different histograms
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

#to print kaon, proton, or pion 
def printing(number):
	if number == 1:
		print("for kaon")
	elif number == 2:
		print("for proton")
	else:
		print("for pion")

#pressure scan 
def pressureScan(determiner,stat_array,mass,numRings,stat5p_array,stat6p_array,stat7p_array,stat8p_array):
	printing(determiner)
	for P in range(P_min,P_max,NumPress):
		#initialization per pressure 
		stat5p = 0 #for Nc >= 5 stat 5+ 
		stat6p = 0 #for Nc >= 6 stat 6+
		stat7p = 0 #for Nc >= 7 stat 7+
		stat8p = 0 #for Nc = 8 stat 8 
		n = index(P)
		r = radius(mka,n) 
		#check if radius is within constraints
		stat = check(r,P,mka,n,r_min,r_max,P_min,P_max)
		#store the pressures if true  
		if stat == 1: 
			stat_array = append(stat_array, P)
		#for every pressure generate photons (more rings)
		#all radii are used because off set of aperture and smearing of particle rings
		#allows for certain invalid radii to still be inside the aperture to some degree 
		for i in range(0,numRings):
			Nc,up,down,left,right = photonGeneration(r,mass)
			stat5p = valuesHist5p(Nc) 
			stat6p = valuesHist6p(Nc) 
			stat7p = valuesHist7p(Nc) 
			stat8p = valuesHist8p(Nc) 
			if stat5p == 1:
				stat5p_array = append(stat5p_array, P)
			if stat6p == 1:
				stat6p_array = append(stat6p_array, P)
			if stat7p == 1:
				stat7p_array = append(stat7p_array, P)
			if stat8p == 1:
				stat8p_array = append(stat8p_array, P)
	print(stat_array)
	return stat_array,stat5p_array,stat6p_array,stat7p_array,stat8p_array 

#pressure selection for alignment (rings graph)
def pressureSelection(determiner,mass,numRings,hue,name): 
	printing(determiner)
	#sections: up, down, left, right
	#sectors 1 and 8 = up
	#sectors 4 and 5 = down 
	#sectors 2 and 3 = right
	#sectors 6 and 7 = left 
	#initialization 
	sumup = 0 
	sumdown = 0
	sumleft = 0 
	sumright = 0 
	for P in range(P_min,P_max,NumPress):
		if P != P_select:
			continue
		up = 0 
		down = 0
		left = 0
		right = 0
		totalcount = 0
		print(n)
		r = radius(mass,n)
		for i in range(0,numRings):
			Nc,up,down,left,right = photonGeneration(r,mass) 		
			sumup += up
			sumdown += down
			sumleft += left
			sumright += right
	totalcount = sumup + sumdown + sumleft + sumright
	print("up:",sumup)
	print("down:",sumdown)
	print("left:",sumleft)
	print("right:",sumright)
	print("total count:",totalcount)
	getAsymmetries(sumup,sumdown,sumleft,sumright,totalcount) 
	ring = circle((0,0),r-5,fill=False,color=hue,label=name,alpha=1)
	return ring 



#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#initialization of parameters (unfixed, fixed, and default)




	#start the timer
start_time = time.time()

	#declare unfixed parameters
x_off = 0.0 #cm shift in x of the photons (opposite of shifting the aperture)
y_off = 0.0 #cm shift in y of the photons (opposite of shifting the aperture)
P_select = 1.8E+6 #baryes the selected pressure 
sigma = 0.01 #cm smearing of photons per particle ring
T = 293 #K 
l = 300 #cm or 3m distance of cone 
NumPress = int(2E+3) #number of pressures to generate for the pressure scan
doScan = 0 #1 means do pressure scan while 0 means do pressure select
useN2 = 1 #1 means use N2 gas while 0 means use H2 gas  

	#declare fixed parameters 
P0 = 1E+6 #baryes 1E+5 Pa 1 bar  
T0 = 273.15 #K 0 deg celsius 
c = 29979245800 #cm/s 
factor = 1.60217733E-12 #erg / eV 
mpi = 139.57039E+6 #eV 140MeV/c mass of +pi
mka = 493.677E+6 #eV 494MeV/c mass of +kaon
mpr = 938.27208816E+6 #eV 938MeV/c mass of proton 
#convert masses from eV/c to g 
mpi = mpi*factor/c**2 
mka = mka*factor/c**2 
mpr = mpr*factor/c**2
theta_min = 0
theta_max = pi*2
#momentum is FIXED 
p = (75*(1E+9))*factor/c #g cm/s or 75 GeV/c 
#70% pion, 24% proton, and 6% kaon rings generated per pressure 
num1 = int(6E+2) #number of kaon rings to generate per pressure in pressure scan
num2 = int(24E+2) #number of proton rings to generate per pressure in pressure scan
num3 = int(70E+2) #number of pion rings to generate per pressure in pressure scan
#pressures for each gas to set the size of the ring-like aperture 
P_N2 = 1.75E+6 #baryes 1.75 bar of N2 gas 
P_H2 = 0 #baryes 0 bar of H2 gas 
dev = 0.09 #cm deviation from kaon radius to calculate the size of ring-like aperture 
name1 = "kaon"
name2 = "proton"
name3 = "pion"
color1 = "tab:blue"
color2 = "tab:orange"
color3 = "tab:green"

	#set defaults for parameters that will be changed 
n = 1 #index of refraction will be changed by changing pressure via refractive index 
P = 1E+6 #baryes 1 bar
r = 0 #radius of ring
#for 1,2,3 or kaon, proton, pion respectively 
stat1 = 0 #boolean that is either 1 or 0 depending if the ring is in the boundary or not respectively
stat2 = 0 
stat3 = 0
#create empty arrays for the histograms 
#array of P list
vec1 = array([]) #empty 
vec2 = array([]) #empty 
vec3 = array([]) #empty 
#array of Nc 5+, 6+, 7+, 8 
#Nc means number of coincidences 
vec5p1 = array([]) 
vec6p1 = array([]) 
vec7p1 = array([]) 
vec8p1 = array([]) 
vec5p2 = array([]) 
vec6p2 = array([]) 
vec7p2 = array([]) 
vec8p2 = array([]) 
vec5p3 = array([]) 
vec6p3 = array([]) 
vec7p3 = array([]) 
vec8p3 = array([]) 

    #print the unfixed parameters
print("\nParameters:")
print("temperature:",T,"K")
print("length:",l,"cm")
print("x_off:",x_off,"cm")
print("y_off:",y_off,"cm")
print("sigma:",sigma,"cm") 
print("number of pressures for each scan:",NumPress)
print("selected pressure:",P_select/1E+6,"bar")

	#set parameters according to chosen gas 
if useN2 == 1:
	#for N2 gas 
	print("radiator: N2")
	#set the n0 refractive index at P0 and T0 according to the website
	#https://refractiveindex.info/?shelf=main&book=N2&page=Borzsonyi
	#assuming 0.42 micrometer wavelength of cherenkov radiation 
	n0 = 1.00029910 
	#calculate the ring-like aperature radii (inner and outer) thus the size
	#using the kaon ring at a particular pressure as a base
	n = index(P_N2)
	r_kaon = radius(mka,n)
	r_min = r_kaon - dev #minimum radius of the aperture 
	r_max= r_kaon + dev #maximum radius of aperture 
	apertureSize = r_max - r_min 
	#set the max and min pressures for the scan 
	P_max = int(2E+6) #baryes 2 bar 
	P_min = int(1.6E+6) #baryes 1.6 bar

elif useN2 == 0:
	#for H2 gas 
	print("radiator: H2")
	#set the n0 refractive index at P0 and T0 according to the website 
	#https://refractiveindex.info/?shelf=main&book=H2&page=Koch 
	#assuming 0.42 micrometer wavelength of cherenkov radiation 
	n0 = 1.00014231 
	#calculate the ring-like aperature radii (inner and outer) thus the size
	#using the kaon ring at a particular pressure as a base
	n = index(P_H2)
	r_kaon = radius(mka,n)
	r_min = r_kaon - dev #minimum radius of the aperture 
	r_max= r_kaon + dev #maximum radius of aperture 
	apertureSize = r_max - r_min 
    #set the max and min pressures for the scan 
	P_max = 4.2E+6 #baryes 4.2 bar 
	P_min = 3.55E+6 #baryes 3.55 bar

	#print the new parameters 
print("r_min:",r_min,"cm")
print("r_max:",r_max,"cm")
print("P_min: ",P_min/1E+6,"bar")
print("P_max: ",P_max/1E+6,"bar")
print("aperture size:",apertureSize,"cm\n")


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#computation 


	#perform the pressure scan
if doScan == 1:
	print("Performing Pressure Scan:\n")
	vec1,vec5p1,vec6p1,vec7p1,vec8p1= pressureScan(1,vec1,mka,num1,vec5p1,vec6p1,vec7p1,vec8p1)
	vec2,vec5p2,vec6p2,vec7p2,vec8p2= pressureScan(2,vec2,mpr,num2,vec5p2,vec6p2,vec7p2,vec8p2)
	vec3,vec5p3,vec6p3,vec7p3,vec8p3= pressureScan(3,vec3,mpr,num3,vec5p3,vec6p3,vec7p3,vec8p3)

	#perform the selected pressure alignment (rings graph)
elif doScan == 0:
	print("Performing Pressure Selection:\n")
	kaonring = pressureSelection(1,mka,num1,color1,name1)
	protonring = pressureSelection(2,mpr,num2,color1,name2)
	pionring = pressureSelection(3,mpi,num3,color1,name3)

#end the timer
print("%s seconds" % (time.time()-start_time))