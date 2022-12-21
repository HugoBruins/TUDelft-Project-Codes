import numpy as np
import matplotlib.pyplot as plt
import math

resolution = 1000
thickness_amount=10    #amount of thickness used in calculations

# Axial and lateral acceleration
a_axial = 6*9.81
a_lateral = 2*9.81

total_length_tanks = 2.327  # total length of the tanks next to each other
r = 0.352                   # radius fuel tank
m1 = 500.3 + 9.05/2         # mass first fuel tank (N2O4)
m2 = 303.2 + 9.95/2         # mass second fuel tank (MMH)
L = 3                       # length cilinder
a0 = 0                      # distance of the first tank to the floor

length_single_tank = total_length_tanks/2
a = a0 + length_single_tank/2
b = a + length_single_tank

# Forces caused by acceleration
P1z = a_axial*m1
P2z = a_axial*m2
P1y = a_lateral*m1
P2y = a_lateral*m2


z = np.linspace(0, L, resolution)

# Getting the first force and moment due to compatibility
A = np.array([[L**2/2, L],
              [L**3/6, L**2/2]])
A_inverse = np.linalg.inv(A)
junk = np.array([P1y/2 * (L-a)**2 + P2y/2 * (L-b)**2, 
                 P1y/6 * (L-a)**3 + P2y/6 * (L-b)**3])
V = A_inverse.dot(junk)
Ay = V[0]
MA = V[1]

# Calculating the internal normal force as a function of z
def calculate_normal_force(z):
    force =     P1z*(1-a/L) + P2z*(1-b/L) \
                - P1z*(z-a)**0*np.heaviside(z-a, 0) \
                - P2z*(z-b)**0*np.heaviside(z-b, 0)
    return force

# Calculating the internal shear force as a function of z
def calculate_shear_force(z):
    force =     Ay \
                - P1y*(z-a)**0*np.heaviside(z-a, 0) \
                - P2y*(z-b)**0*np.heaviside(z-b, 0)
    return force

# Calculating the internal moment about x axis as a function of z
def calculate_internal_moment(z):
    moment =    MA + Ay*z \
                - P1y*(z-a)**1*np.heaviside(z-a, 0) \
                - P2y*(z-b)**1*np.heaviside(z-b, 0)
    return moment

# Calculating the theta as a function of z, to get actual theta multiply this
# by stiffness E*I
def calculate_thetaEI(z):
    theta = MA*z + 0.5*Ay*z**2 \
            - 0.5*P1y*(z-a)**2*np.heaviside(z-a, 0) \
            - 0.5*P2y*(z-b)**2*np.heaviside(z-b, 0)
    return -theta

# Calculating the deflection as a function of z, to get actual deflection
# multiply this by stiffness E*I
def calculate_vEI(z):
    v = 0.5*MA*z**2 + 1/6*Ay*z**3 \
        - 1/6*P1y*(z-a)**3*np.heaviside(z-a, 0) \
        - 1/6*P2y*(z-b)**3*np.heaviside(z-b, 0)
    return -v

N = calculate_normal_force(z)
Vy = calculate_shear_force(z)
Mx = calculate_internal_moment(z)
theta = calculate_thetaEI(z)
v = calculate_vEI(z)

# plt.subplot(511)
# plt.plot(z,N, label="Internal normal force")
# plt.legend()
# plt.grid()
# plt.subplot(512)
# plt.plot(z,Vy, label="Internal Vy")
# plt.legend()
# plt.grid()
# plt.subplot(513)
# plt.plot(z,Mx, label="Internal Mx")
# plt.legend()
# plt.grid()
# plt.subplot(514)
# plt.plot(z,theta, label="Theta*E*Ixx")
# plt.legend()
# plt.grid()
# plt.subplot(515)
# plt.plot(z,v, label="v*E*Ixx")
# plt.legend()
# plt.grid()
# plt.show()


#Internal Stress Calculations:

sigma_yield=241E6   #Normal Yield strength of material
tau_yield=207E6     #Shear Strength of material
E=68E9          #Young's Modulus of material
rho=2720           #Density of material

t=0
dt=0.001
sigma_z_max=sigma_yield+1
tau_z_max=tau_yield+1
while (abs(sigma_z_max)>=sigma_yield \
       or abs(tau_z_max)>=tau_yield \
       or abs(sigma_z_max)>=sigma_cr) and t<0.5:
    t+=dt
    
    #normal stress check
    sigma_z_normal=N/(2*math.pi*r*t)
    sigma_z_bending=64*Mx/(math.pi*t*r**2)
    sigma_z=sigma_z_normal+sigma_z_bending
    sigma_z_max=max(sigma_z)
    
    #shear stress check
    tau_z=t*Vy*64/(math.pi*r)
    tau_z_max=max(tau_z)
    
    #buckling check
    Ixx=math.pi*(2*r)**3*t/8
    r_gyr=math.sqrt(Ixx/(2*math.pi*r*t))
    sigma_cr=(math.pi)**2*E/((0.5*3/r_gyr)**2)
    
if t>=0.5:
    print("material incompatible with design, thickness becomes unrealistic")
else:
    print("thickness required is:",t)
    
mass=rho*3*2*math.pi*r*t
print("with this thickness and material, the mass of the cylinder becomes:",mass)
