"""
This code was created to size a helium tank in an appendix of a project 
Written by: Hugo Bruins
"""
import numpy as np
import matplotlib.pyplot as plt

# a = 0.0346 barL^2/mol = 0.00346 Pa m^3/mol
a = 0.00346
# b = 0.0238 L/mol = 0.0238e-3 m^3/mol
b = 0.0238e-3

p = np.arange(100*10**5, 450*10**5, 1000)
P = 20e5
V = 0.721 # 0.362 + 0.359, total tank volume
T = 290
R = 8.31446261815324
ntanks = 1

#ti-6al-4v-grade-5
# Yield strength is reduced by a safety factor of 1.25
sigma = 880e6
density = 4430

v = (V*P*(b*p+R*T))/(R*T*(p-P))
# amount of moles of helium
n = (P*(v+V))/(P*b + R*T) 

v /= ntanks
n /= ntanks

# Radius of the pressure tank
r = np.cbrt((3*v)/(4*np.pi))
# Thickness of the pressure tank
t = (p*r)/(2*sigma)
# Volume of the material of the pressure tank
tankvolume = 4*np.pi*r**2*t
# Mass of the pressure tank
tankmass = tankvolume * density
wettankmass = tankmass + n * 0.004
# To compare the error that comes with assuming a = 0
actualpressure = (n*R*T)/(v-n*b) - a*n**2/v**2

def plot(x,y, plotnr, xlabel = "", ylabel ="", plotlabel= ""):
    plt.subplot(plotnr)
    plt.plot(x,y, label=plotlabel)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid()

plt.show()
plot(p*1e-5, (actualpressure-p)*1e-5, 221, "pressure [bar]", "Pressure error [bar]")
plt.ylim(-8,0)
plot(p*1e-5, v*1000, 222, "pressure [bar]", "Required tank volume [L]")
plt.ylim(0, 200)
plot(p*1e-5, t*1000, 223, "pressure [bar]", "Wall thickness [mm]")
plt.ylim(0,6)
plot(p*1e-5, wettankmass, 224,  "pressure [bar]", "Wet Mass [kg]")

for i in range(len(p)):
    if p[i] == 175e5 or p[i] == 225e5 or p[i] == 300e5 or p[i] == 400e5:
        print("tank starting pressure [bar]:    ",round(p[i]*1e-5))
        print("tank volume [L]:                 ",np.ceil(v[i]*10000)/10)
        print("Amount of helium [mol]:          ",np.ceil(n[i]*10)/10)
        print("Tank radius [mm]:                ",np.ceil(r[i]*10000)/10)
        print("Tank thickness [mm]:             ",np.ceil(t[i]*100000)/100)
        print("Tank dry mass [kg]:              ",np.ceil(tankmass[i]*100)/100)
        print("Helium mass [kg]:                ",np.ceil(n[i]*0.004*10)/10)
        print("Tank wet mass [kg]:              ",np.ceil(wettankmass[i]*100)/100)
        print("")
