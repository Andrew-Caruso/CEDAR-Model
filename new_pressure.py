#import necessary modules 
from os import POSIX_SPAWN_OPEN
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
from numpy import max as max 
from numpy import pi as pi 
from numpy import random as rand 
from numpy import isnan as isnan 
from numpy import nan as nan 
from numpy import reshape as reshape
from numpy import abs as abs 
from numpy import arange as arange 
from numpy import histogram2d as hist2d 
from matplotlib.patches import Circle as circle 
from matplotlib.patches import Wedge as wedge 
from matplotlib.lines import Line2D as line
#layout of code in sections:
#start timer 
#declaration of functions
#initalization of parameters
#computation 
#stop timer 

	#start timer
start_time = time.time()


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
def photonGeneration(r,m,a_min,a_max,numG): 
		#for the photon scan 
	#numG is the number of photons to generate on the ring 
	#numG may be set to 18 e.g. but use poission distribution to pick a value near or at it
	numG = rand.poisson(numG) #number of photons for the time resolution from a poission distribution located a numG instead of using rand numbers 
	tres1 = 280 #picoseconds, time resolution for 1 photon, 280 ps  
	mu1 = 0 #mean for normal/gaussian distribution   
	tres = tres1 / sqrt(n) #time resolution for any number of photons 
	t = rand.normal(mu1,tres) #mean of zero, std of the time resolution (given an numG)

		#for the diaphragm scan, pressure scan, and the pressure selection  
	chart = zeros(8) #index = sector-1, used as binary 0 or 1 
	#for a PM of the sector detecting 1 or more photons 
	#from sector 1 (index 0) to sector 8 (index 7)
	count = zeros(8) #used as number of photons in each sector 
	mu2 = r #mean for normal/gaussian distribution 
	left = 0
	right = 0
	up = 0
	down = 0
	#loop through each photon 
	for i in range(0,numG,1): 
		phi = rnd.uniform(0,2*pi) #generate random angle
		r0 = rand.normal(mu2,sigma) #generate photon's radius from center
		x = r0*cos(phi) + x_off 
		y = r0*sin(phi) + y_off 
		r1 = sqrt(x**2 + y**2) 
		#if the photon radius is no longer in the ring, throw the iteration out and move on to the next photon  
		if r1 < a_min or r1 > a_max:
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
	return totCount,up,down,left,right,t 

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
	print("right ratio:",ratioright)
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
		print("\tfor kaon")
	elif number == 2:
		print("\tfor proton")
	else:
		print("\tfor pion")

#pressure scan 
def pressureScan(determiner,statvec,mass,numRings,stat5pvec,stat6pvec,stat7pvec,stat8pvec):
	printing(determiner)
	for P in bins1:  
		#initialization per pressure 
		stat5p = 0 #for Nc >= 5 stat 5+ 
		stat6p = 0 #for Nc >= 6 stat 6+
		stat7p = 0 #for Nc >= 7 stat 7+
		stat8p = 0 #for Nc = 8 stat 8 
		n = index(P)
		r = radius(mass,n) 
		#check if radius is within constraints
		stat = check(r,P,mass,n,r_min,r_max,P_min,P_max)
		#store frequencies of the pressures for the histogram values array
		statindex = int((int(P) - P_min) / step)
		if stat == 1: 
			statvec[statindex] += 1
		#for every pressure generate photons (more rings)
		#all radii are used because off set of aperture and smearing of particle rings
		#allows for certain invalid radii to still be inside the aperture to some degree 
		for i in range(0,numRings):
			Nc,up,down,left,right,_ = photonGeneration(r,mass,r_min,r_max,18)
			stat5p = valuesHist5p(Nc) 
			stat6p = valuesHist6p(Nc) 
			stat7p = valuesHist7p(Nc) 
			stat8p = valuesHist8p(Nc) 
			#store the frequency in each array (for each coincidence)
			#each array is of length numPress, number of pressures to generate, so each slot or pressure of the array 
			#corresponds to the frequency of 5+,6+,7+ or 8 coincidence at said slot or pressure  
			if stat5p == 1:
				stat5pvec[statindex] += 1
			if stat6p == 1:
				stat6pvec[statindex] += 1
			if stat7p == 1:
				stat7pvec[statindex] += 1
			if stat8p == 1:
				stat8pvec[statindex] += 1
	return statvec,stat5pvec,stat6pvec,stat7pvec,stat8pvec

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
		r = radius(mass,n)
		print("radius:",r,"cm")
		for i in range(0,numRings):
			Nc,up,down,left,right,_ = photonGeneration(r,mass,r_min,r_max,18) 		
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
	return ring, sumup, sumdown,sumleft,sumright 
	
#vary the width of aperture, thus a diaphragm scan at a fixed pressure
def diaphragmScan(w_min,w_max,mass):
	stepval = 0.005 
	i = int(0)
	widths = arange(w_min,w_max,stepval)
	coincidences = zeros(len(widths))
	for w in widths:
		#these are in photonGeneration 
		a_min = r_kaon - w/2.0
		a_max = r_kaon + w/2.0
		for x in range(0,1000):
			num,_,_,_,_,_ = photonGeneration(r_kaon,mass,a_min,a_max,18)
			if num >= 6:
				coincidences[i] += 1 
			#number of coincidences is total number of photons detected in any of the sectors 
		i += int(1) 
	return widths,coincidences
 

#photon scan 
def photonScan(numG,numR,bool,numG_array):
	#numG is the number of photons 
	#numR is the number of rings 
	#numC is the number of coincidences (total number of photons detected in any sector), rather the number of sectors that made a detection
	#bool is either 1 to make single plots or 0 to not make single plots
	#photonScan has two modes bool=1 to make single histogram for only 1 numG and bool=0 to use outputs from multiple undrawn plots to make a 2D histogram
	radius = r_kaon
	mass = mka
	b_min = r_min
	b_max = r_max
	num_it = int(numG*numR) #number of iterations
	print("number of interations:",num_it)
	coin_array= zeros(num_it)  
	time_array = zeros(num_it)	
	fig6,ax6 = plt.subplots(1,2) #create a figure with 1 row 2 columns, i.e. 2 horizontal axes 
	#mode 1: create single histogram for single numG
	if bool == 1:
		#generate 2 histograms for single numG
		print("Selected photon scan mode 1 (1D histogram)")
		num = numG 
		for i in range(0,numR): 
			numC,_,_,_,_,time= photonGeneration(radius,mass,b_min,b_max,numG)
			coin_array[i] = numC  
			time_array[i] = time
		bins_c,val_c,_ = ax6[0].hist(coin_array,bins_Gscan)
		x_max1 = max(val_c)
		y_max1 = max(bins_c)
		print("sector max:",x_max1,",",y_max1)
		ax6[0].annotate("numG: %i" % num,(x_max1*0.9,y_max1*0.9),ha='center',va='center')
		ax6[0].set_xlabel("Number of sectors")
		ax6[0].set_ylabel("Frequency")
		#create the time histogram 
		bins_t,val_t,_ = ax6[1].hist(time_array,bins_Gscan)
		x_max2 = max(val_t)
		y_max2 = max(bins_t)
		print("time max:",x_max2,",",y_max2)
		ax6[1].annotate("numG: %i" % num,(x_max2*0.8,y_max2*0.9),ha='center',va='center')
		ax6[1].set_xlabel("time (ps)")

	#mode 2: use outputs from many histograms and many numG to create 2D histogram
	else: 	
		#generate 2D histograms for many numG 
		print("Selected photon scan mode 2 (2D histogram)")
		coin_array= zeros(numG)  
		time_array = zeros(numG)	
		print("coin shape:",shape(coin_array))
		print("time shape:",shape(time_array))
		print("numG shape:",shape(numG_array))
		for num in range(1,numG+1):
			for i in range(0,numR): 
				numC,_,_,_,_,time = photonGeneration(radius,mass,b_min,b_max,num)
				coin_array[num-1] = numC  
				time_array[num-1] = time
		#reshape the coin and time arrays 
		coin_array = reshape(coin_array,(numG,-1))
		time_array = reshape(time_array,(numG,-1))
		print("coin shape:",shape(coin_array))
		print("time shape:",shape(time_array))
		print("numG shape:",shape(numG_array))
		'''
		ax6[1].hist2d(numG_array,coin_array)
		ax6[0].hist2d(numG_array,time_array) 
		'''


	return coin_array,time_array 

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#initialization of parameters (unfixed, fixed, and default)



	#declare unfixed parameters
x_off = 0.0 #cm shift in x of the photons (opposite of shifting the aperture)
y_off = 0.0 #cm shift in y of the photons (opposite of shifting the aperture)
sigma = 0.02 #cm smearing of photons per particle ring
T = 293 #K 
l = 300 #cm or 3m distance of cone 
NumPress = int(2E+2) #number of pressures to generate for the pressure scan
useN2 = 1 #1 means use N2 gas while 0 means use H2 gas  
doScan = 3
#0 for pressure scan
#1 for pressure selection (alignment)
#2 for diaphragm scan
#3 for photon scan 

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
P_H2 = 3.8E+6 #baryes 3.8 bar of H2 gas 
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
#create zeros numpy arrays for the histograms 
#array of frequencies of pressures for each particle  
vec1 = zeros(NumPress) 
vec2 = zeros(NumPress) 
vec3 = zeros(NumPress) 
#array of Nc 5+, 6+, 7+, 8 for each particle  
#Nc denotes number of coincidences 
vec5p1 = zeros(NumPress) 
vec6p1 = zeros(NumPress) 
vec7p1 = zeros(NumPress) 
vec8p1 = zeros(NumPress) 
vec5p2 = zeros(NumPress) 
vec6p2 = zeros(NumPress) 
vec7p2 = zeros(NumPress) 
vec8p2 = zeros(NumPress) 
vec5p3 = zeros(NumPress) 
vec6p3 = zeros(NumPress) 
vec7p3 = zeros(NumPress) 
vec8p3 = zeros(NumPress) 

    #print the unfixed parameters
print("\n\tParameters:")
print("temperature:",T,"K")
print("length:",l,"cm")
print("x_off:",x_off,"cm")
print("y_off:",y_off,"cm")
print("sigma:",sigma,"cm") 
print("number of steps:",NumPress)

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
	width0 = r_max - r_min #aperture size i.e. width of diaphragm  
	#set the max and min pressures for the scan 
	P_max = int(2E+6) #baryes 2 bar 
	P_min = int(1.6E+6) #baryes 1.6 bar
	step = (P_max - P_min) / NumPress 
	P_select = 1.8E+6 #baryes or 1.8 bar the selected pressure 

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
	width0 = r_max - r_min 
    #set the max and min pressures for the scan 
	P_max = int(4.35E+6) #baryes 4.4 bar 
	P_min = int(3.55E+6) #baryes 3.55 bar
	step = (P_max - P_min) / NumPress 
	P_select = 3.7E+6 #baryes or 1.8 bar the selected pressure 

#create bins array for all histograms
bins1 = arange(P_min,P_max,step) 

	#print the new parameters 
print("step size:",step)
print("r_min:",r_min,"cm")
print("r_max:",r_max,"cm")
print("aperture size:",width0,"cm")
print("P_min: ",P_min/1E+6,"bar")
print("P_max: ",P_max/1E+6,"bar")
print("selected pressure:",P_select/1E+6,"bar\n")


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#computation 

	#perform the pressure scan
if doScan == 0:
	print("\tPerforming Pressure Scan:")
	vec1,vec5p1,vec6p1,vec7p1,vec8p1= pressureScan(1,vec1,mka,num1,vec5p1,vec6p1,vec7p1,vec8p1)
	vec2,vec5p2,vec6p2,vec7p2,vec8p2= pressureScan(2,vec2,mpr,num2,vec5p2,vec6p2,vec7p2,vec8p2)
	vec3,vec5p3,vec6p3,vec7p3,vec8p3= pressureScan(3,vec3,mpi,num3,vec5p3,vec6p3,vec7p3,vec8p3)
	#create the first histogram 
	#manually make the histograms using numpy.bar instead of numpy.hist
	fig1, ax1 = plt.subplots(1,1)
	ax1.bar(bins1,vec1,width=step,label=name1,color=color1,align='edge')
	ax1.bar(bins1,vec2,width=step,label=name2,color=color2,align='edge')
	ax1.bar(bins1,vec3,width=step,label=name3,color=color3,align='edge')
	ax1.set_xlabel("pressure (Ba)")
	ax1.set_ylabel("frequency")
	ax1.legend(loc="best") 
	ax1.set_title("Pressure Scan of Particles")
	#create the second histogram 
	#manually make the histograms using numpy.bar instead of numpy.hist
	#create histograms for Nc 5+,6+,7+,8+ for each particle 
	fig3,ax3 = plt.subplots(1,1) 
	ax3.bar(bins1,vec5p1,width=step,label="5+",color='magenta')
	ax3.bar(bins1,vec6p1,width=step,label="6+",color='darkred')
	ax3.bar(bins1,vec7p1,width=step,label="7+",color='lime')
	ax3.bar(bins1,vec8p1,width=step,label="8",color='cyan')
	ax3.bar(bins1,vec5p2,width=step,color='magenta')
	ax3.bar(bins1,vec6p2,width=step,color='darkred')
	ax3.bar(bins1,vec7p2,width=step,color='lime')
	ax3.bar(bins1,vec8p2,width=step,color='cyan')
	ax3.bar(bins1,vec5p3,width=step,color='magenta')
	ax3.bar(bins1,vec6p3,width=step,color='darkred')
	ax3.bar(bins1,vec7p3,width=step,color='lime')
	ax3.bar(bins1,vec8p3,width=step,color='cyan')
	ax3.set_xlabel("pressure (Ba)")
	ax3.set_ylabel("frequency")
	ax3.set_title("x_off: %fcm y_off: %fcm std: %fcm" % (x_off,y_off,sigma))
	ax3.legend(loc="best") 

	#perform the selected pressure alignment (rings graph)
elif doScan == 1:
	print("\tPerforming Pressure Selection:")
	kaonring,up1,down1,left1,right1 = pressureSelection(1,mka,num1,color1,name1)
	protonring,up2,down2,left2,right2= pressureSelection(2,mpr,num2,color2,name2)
	pionring,up3,down3,left3,right3= pressureSelection(3,mpi,num3,color3,name3)
	#get the total asymmetries include all three particles
	totup = up1+up2+up3
	totdown = down1+down2+down3
	totleft = left1+left2+left3
	totright = right1+right2+right3
	totnum = totup+totdown+totleft+totright
	print("\tfor total")
	print("up: ",totup)
	print("down: ",totdown)
	print("left: ",totleft)
	print("right: ",totright)
	getAsymmetries(totup,totdown,totleft,totright,totnum)
	#create the rings graph
	#NOTE: all radii for fig4 have 5 cm removed to make scaling better for visual
	#start of fig4 (for single pressure, NOT for pressure scan or monte carlo) 
	fig4,ax4 = plt.subplots(1,1) 
	originPoint= circle((0,0),0.05,fill=True,color="black",alpha=1)
	maxring = circle((-x_off,-y_off),r_max-5,fill=False,color="black",label="ringlike aperture",alpha=1)
	minring = circle((-x_off,-y_off),r_min-5,fill=False,color="black",alpha=1)
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
	ax4.set_title("x_off: %fcm y_off: %fcm std: %fcm" % (x_off,y_off,sigma))	
	ax4.legend(loc="best", prop={'size':7}) 
	#end of fig4

	#perform the diaphragm scan
	#want to replicate Fig 6.5 page 156 from "Leptonic Decays and Kaon Identification at the NA62 Experiment at CERN"
	#by Angela Romano, thesis submitted to the Univerisyt of Birmingham for the degree of Doctor of Philosophy 
	#as width increased from 0.4 mm to 2 mm 
	#P = P_max of the pion peak 
if doScan == 2:
	print("\tPerforming Diaphragm Scan:")
	'''
	#perform pressure scan of pion 
	_,_,vec6p3,_,_= pressureScan(3,vec3,mpi,num3,vec5p3,vec6p3,vec7p3,vec8p3)
	#get the frequency at the peak 
	#get the pressure at the peak
	peak3 = where(vec6p3 == max(vec6p3)) #returns a list from the array 
	peak3_index = peak3[0][0] #get the list then return the elemenet of the list  
	print("peak index:",peak3_index)
	'''
	P_select = 1.75E+6 #baryes 1.75 bar
	num6p3 = max(vec6p3)
	width_min = 0.0 #cm  
	width_max = 0.2 #cm  
	#print the set parameters 
	print("peak frequency: ",num6p3)
	print("peak pressure:",P_select/1E+6,"bar")
	print("min width: ",width_min)
	print("max width: ",width_max)
	#perform the diaphragm scan at the selected pressure (pressure at the peak)
	widthVal, numCoin = diaphragmScan(width_min,width_max,mka)
	'''
	 #create the single histogram for pion 6+ coincidence 
	#to confirm the numbers 
	fig2, ax2 = plt.subplots(1,1)
	ax2.bar(bins,vec6p3,width=step,label="6+",color='darkred')
	ax2.set_xlabel("pressure (Ba)")
	ax2.set_ylabel("frequency")
	ax2.set_title("x_off: %fcm y_off: %fcm std: %fcm" % (x_off,y_off,sigma))
	ax2.legend(loc="best") 
	'''
	#create the effifiency vs width plot
	fig5, ax5 = plt.subplots(1,1)
	ax5.plot(widthVal,numCoin)
	ax5.set_xlabel("diaphragm width (cm)")
	ax5.set_ylabel("number of coincidences")
	ax5.set_title("x_off: %fcm y_off: %fcm std: %fcm" % (x_off,y_off,sigma))	
	plt.show()

	#perform photon scan
if doScan == 3:
	numP = int(18) #number of photons
	numRings = 1000 #number of rings to generate per photon 
	bins_Gscan = 100 
	gen_single = 0 #to generate histogram for single numG else, generate histogram for many numG  
	print("\tPerforming Photon Scan")
	print("Number of photons:",numP)
	print("Number of rings per photon:",numRings)
	print("bins:",bins_Gscan)
	numPVec = arange(1,numP+1) #from 1 to numP
	photonScan(numP,numRings,gen_single,numPVec)


	#stop timer
print("\n\tRun time:")
print("%s seconds" % (time.time()-start_time))
print("%s minutes" % ((time.time()-start_time)/60))
#plt.show()
