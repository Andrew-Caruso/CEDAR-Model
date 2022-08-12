#find the kaon ring assuming all other parameters are ideal
#neglect the mirror and circular-ring aperture 
#simplest model: use the cherenkov angle to find the radius of the cone produced (neglect the mirror and diaphram). Simplest model. 
import numpy as np
import matplotlib.pyplot as plt 
#units of choice: CGS  

#need to convert the following (since they are all in natural units where c=1 & hbar=1)  
#by using 1eV = 1.60218E-12erg, 1erg=gcm^2/s^2, E=mc^2, p=mv, hbar=1.05erg*s=1, c=3E10cm/s 
p_fixed = 75E+9 #eV 75GeV/c momentum for +pi, +kaon, and proton  
mpi = 139.57039E+6 #eV 140Mev/c mass of +pi
mka = 493.677E+6 #eV 494MeV/c mass of +kaon
mpr = 938.27208816E+6 #eV 938MeV/c mass of proton 
mel = 0.51099895000E+6 #eV 0.51 MeV/c mass of positron 
mmu = 105.6583755E+6 #eV 106 MeV/c mass of muon 
c = 29979245800 #cm/s 
hbar = 1.05457266E-27 #g cm^2 / s 
#1 eV = 1.60217733E-12 g cm^2 / s^2
factor = 1.60217733E-12 #erg / eV 
p_fixed = p_fixed*factor/c 
mpi = mpi*factor/(c**2) 
mka = mka*factor/(c**2) 
mpr  = mpr*factor/(c**2) 
mel = mel*factor/(c**2)
mmu = mmu*factor/(c**2)
l = 300 #cm or 3m distance of cone 
P = 1 #bar big P pressure for the particular index of refraction in "CERN-82-13.pdf" page 16 


#set fixed parameters
radiator = "N2" 
print("radiator:",radiator)
n = 1.0004775549188156 #refractive index property of the gas 
#from page 16, from n/(n-1) = 2095 at pressure of 1 bar 
#depends on pressure, temperature, and wavelength 
D = 0.18 #cm or 1.8mm diameter of diaphram aperture from 0.03 to 20 mm OR 0.003 cm to 2 cm
print("refractive index:",n)
print("momentum p_fixed:",p_fixed,"g cm/s")
print("pressure:",P,"bar\n")

#define radius function for brevity 
def radius(m,p):
	beta = 1 / (np.sqrt(((m*c/p)**2)+1))
	theta = np.arccos(1/(n*beta))
	r = l*np.tan(theta)
	#print(r) 
	return r  

#for kaon 
beta1 = 1 / (np.sqrt(((mka*c/p_fixed)**2)+1))
theta1 = np.arccos(1/(n*beta1))
r1 = l*np.tan(theta1)
print("kaon")
print("mka:",mka,"g")
print("speed of kaon:",beta1*c/c*100,"% c")
print("theta:",theta1,"radians")
print("beta:",beta1) 
print("radius from function for kaon ring:",radius(mka,p_fixed),"cm") 
print("radius of kaon ring:",r1,"cm\n") 

#for pion 
beta2 = 1 / (np.sqrt(((mpi*c/p_fixed)**2)+1))
theta2 = np.arccos(1/(n*beta2))
r2 = l*np.tan(theta2)
print("pion")
print("mpi:",mpi,"g")
print("speed of pion:",beta2*c/c*100,"% c")
print("theta:",theta2,"radians")
print("beta:",beta2) 
print("radius of pion ring:",r2,"cm\n") 

#for proton  
beta3 = 1 / (np.sqrt(((mpr*c/p_fixed)**2)+1))
theta3 = np.arccos(1/(n*beta3))
r3 = l*np.tan(theta3)
print("proton")
print("mpr:",mpr,"g") 
print("speed of proton:",beta3*c/c*100,"% c")
print("theta:",theta3,"radians")
print("beta:",beta3) 
print("radius from function for proton ring:",radius(mpr,p_fixed),"cm") 
print("radius of proton ring:",r3,"cm\n") 

#for positron 
beta4 = 1 / (np.sqrt(((mel*c/p_fixed)**2)+1))
theta4 = np.arccos(1/(n*beta4))
r4 = l*np.tan(theta4)
print("positron")
print("mel:",mel,"g") 
print("speed of positron:",beta4*c/c*100,"% c")
print("theta:",theta4,"radians")
print("beta:",beta4) 
print("radius of positron ring:",r4,"cm\n") 

#for muon 
'''
beta5 = 1 / (np.sqrt(((mmu*c/p_fixed)**2)+1))
theta5 = np.arccos(1/(n*beta5))
r5 = l*np.tan(theta5)
print("muon") 
print("mmu:",mmu,"g")
print("speed of kaon:",beta5*c/c*100,"% c")
print("theta:",theta5,"radians")
print("beta:",beta5)
print("radius of muon ring:",r5,"cm\n") 
'''

#comparison
print("Comparison of different from smallest mass to greatest mass and corresponding radii")
#masses
#express also in natural units where MeV/c -> MeV
print("mel:",mel,"g"," OR 0.5 MeV")
#print("mmu:",mmu,"g"," OR 106 MeV")
print("mpi:",mpi,"g"," OR 140 MeV")
print("mka:",mka,"g"," OR 494 Mev")
print("mpr:",mpr,"g"," OR 938 MeV\n")
#ring radii
print("radius of positron ring:",r4,"cm") 
#print("radius of muon ring:",r5,"cm") 
print("radius of pion ring:",r2,"cm") 
print("radius of kaon ring:",r1,"cm") 
print("radius of proton ring:",r3,"cm") 
print("mass of particle is inversely proportional to ring radius")
#omitted the muon from the code because muons do NOT appear in hadronic beam
#muons only appear in the beam halo, around the beam fro pion and kaon decays 
#see helicity suppreision 



#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------




#second task is to graph radius vs momentum 
#first need to find the max and min momentum values

#solve for the bounds based on the angle restrictions
#cos(theta) = 1/(beta*n), so the lowest defined value of p = pmin
#is when beta*n = 1
#supervisor suggested using mass of proton to find the upper bound and mass of pion to find the lower bound 

N = 3 #to be multplied to p_min2 for p_max 
cgsp_min = mpi*c / np.sqrt(n**2 - 1) #pion lowest momentum, use as low bound for all 
cgsp_min2= mpr*c / np.sqrt(n**2 - 1) #proton lowest momentum, do NOT use as max, will omit proton 
cgsp_max = N*cgsp_min2
#instead let p_max be a mulitple of p_min2 
print("in CGS")
print("\ncgs p_min:",cgsp_min,"g cm/s")
print("cgs p_min2:",cgsp_min2,"g cm/s")
print("cgs p_max:",cgsp_max,"g cm/s")
#reassign p_max and p_min in natural units 
p_min = cgsp_min*c/(factor*1E+9) 
p_min2 = cgsp_min2*c/(factor*1E+9)
p_max = N*p_min2
print("in natural units")
print(p_min,"GeV") 
print(p_min2,"GeV") 
print(p_max,"GeV") 

#create x values 
p_scale = np.linspace(p_min,p_max,100) #for the graph NOT for radius fun
p = np.linspace(cgsp_min,cgsp_max,100) #for radius fun NOT for graph 
#create y values 
#use the radius function 
#1 cm = 10 mm
#NOTE: the radius function is in CGS units NOT natural units
rka = radius(mka,p)*10 
rpi = radius(mpi,p)*10 
rpr = radius(mpr,p)*10 
rmu = radius(mmu,p)*10 
rel = radius(mel,p)*10 
#plot these functions 
plt.plot(p_scale,rka,label="kaon")
plt.plot(p_scale,rpi,label="pion")
plt.plot(p_scale,rpr,label="proton")
#plt.plot(p_scale,rmu,label="muon")
plt.plot(p_scale,rel,label="positron")
plt.xlabel("momentum (GeV)") 
plt.ylabel("radius (cm)") 
plt.title("radius of ring versus momentum") 
plt.legend(loc="best")
plt.xticks(np.arange(0,int(p_max)+1,10))
plt.yticks(np.arange(0,int(p_max)+1,10))
plt.show()


