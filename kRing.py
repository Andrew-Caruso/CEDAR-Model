#find the kaon ring assuming all other parameters are ideal
import numpy as np


#set fixed parameters
radiator = "N2" 
c = 299792458 #m/s 
v = 1 #velocity of particle 
n = 1 #refractive index property of gas 
p = 75E+9 #eV 75GeV/c momentum for +pi, +kaon, and proton  
mpi = 139.57039E+6 #eV 140Mev/c mass of +pi
mka = 493.677E+6 #eV 494MeV/c mass of +kaon
mpr = 938.27208816E+6 #eV 938MeV/c mass of proton 

#set unfixed parameters 
theta = np.arccos(1/(n*beta))
beta = v / c 
#p=gamma*m

