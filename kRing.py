#find the kaon ring assuming all other parameters are ideal
#neglect the mirror and circular-ring aperture 
import numpy as np
#units of choice: CGS  

#need to convert the following (since they are all in natural units where c=1 & hbar=1)  
p = 75E+9 #eV 75GeV/c momentum for +pi, +kaon, and proton  
mpi = 139.57039E+6 #eV 140Mev/c mass of +pi
mka = 493.677E+6 #eV 494MeV/c mass of +kaon
mpr = 938.27208816E+6 #eV 938MeV/c mass of proton 
c = 29979245800 #cm/s 
hbar = 1.05457266E-27 #g cm^2 / s 
#1 eV = 1.60217733E-12 g cm^2 / s^2
factor = 1.60217733E-12 #g cm^2 / s^2  i.e. ergs  
p = p*factor/c 
mpi = mpi*factor/(c**2) 
mka = mka*factor/(c**2) 
mpr  = mpr*factor/(c**2) 
l = 300 #cm or 3m distance of cone 
P = 1 #bar big P pressure for the particular index of refraction in "CERN-82-13.pdf" page 16 

print("p:",p,"g cm/s")
print("mpi:",mpi,"g")
print("mka:",mka,"g")
print("mpr:",mpr,"g\n") 

#set fixed parameters
radiator = "N2" 
print("radiator:",radiator)
n = 1.0004775549188156 #refractive index property of the gas 
#from page 16, from n/(n-1) = 2095 at pressure of 1 bar 
#depends on pressure, temperature, and wavelength 
D = 0.18 #cm or 1.8mm diameter of diaphram aperture from 0.03 to 20 mm OR 0.003 cm to 2 cm
m0 = mka

#set unfixed parameters 
beta = 1 / (np.sqrt(((m0*c/p)**2)+1))
theta = np.arccos(1/(n*beta))
r = l*np.tan(theta)
print("speed of Kaon:",beta*c/c*100,"% c")
print("theta:",theta,"radians")
print("radius of kaon ring:",r,"cm") 

