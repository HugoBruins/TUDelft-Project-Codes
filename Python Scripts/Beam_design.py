import numpy as np
import math

#Characteristics
L=3 #Length of the cylinder (height of spacecraft)
g_axial=  6          #Acceleration in Axial direction
g_transverse= 2      #Acceleration in Transverse Direction
m_fuel_tank= 500       #mass of the fuel tank
m_ox_tank= 500        #mass of the oxidizer tank
R=0.3             #Radius of Cylinder

#Calculation of distance between wall and first tank CoM

a=1

#Axial Loading as Function of z

z = np.linspace(0, L, 100)
P_1=9.81*g_axial*m_fuel_tank     #Load introduced by Fuel Tank
P_2=9.81*g_transverse*m_ox_tank  #Load introduced by Oxidizer Tank
Load_axial=np.zeros(len(z))

for i in range(len(z)):
    if z[i]<a:
        Load_axial[i]=P_1*(1-a/L)+P_2*(1-(a+2*R)/L)
    elif z[i]>=a and z[i]<(a+2*R):
        Load_axial[i] =P_1*(-a/L)+P_2*(1-(a+2*R)/L)
    elif z[i]>= (a + 2*R) and z[i]<=L:
        Load_axial[i] = -P_1*(a/L)-P_2*((a+2*R)/L)



#Lateral Loading as Function of z



#Maximum Stress State as Function of Thickness



#Minimum Required Thickness



#Buckling Check



#Mass Calculation
