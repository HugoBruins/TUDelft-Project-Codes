import numpy as np
import matplotlib.pyplot as plt

resolution = 1000

# Axial and lateral acceleration
a_axial = 6*9.81
a_lateral = 2*9.81

m1 = 700 # mass first fuel tank
m2 = 500 # mass second fuel tank
L = 3    # length cilinder
r = 0.3  # radius cilinder
a = 0.3  # distance from the first tank to the wall

# Forces caused by acceleration
P1z = a_axial*m1
P2z = a_axial*m2
P1y = a_lateral*m1
P2y = a_lateral*m2


z = np.linspace(0, L, resolution)

# Getting the first force and moment due to compatibility
A = np.array([[L**2/2, L],[L**3/6, L**2/2]])
A_inverse = np.linalg.inv(A)
junk = np.array([P1y/2 * (L-a)**2 + P2y/2 * (L-a-2*r)**2, P1y/6 * (L-a)**3 + P2y/6 * (L-a-2*r)**3])
V = A_inverse.dot(junk)
Ay = V[0]
MA = V[1]

# Calculating the internal normal force as a function of z
def calculate_normal_force(z):
    force = P1z*(1-a/L) + P2z*(1-(a+2*r)/L) - P1z*(z-a)**0*np.heaviside(z-a, 0) - P2z*(z-a-2*r)**0*np.heaviside(z-a-2*r, 0)
    return force

# Calculating the internal shear force as a function of z
def calculate_shear_force(z):
    force = Ay - P1y*(z-a)**0*np.heaviside(z-a, 0) - P2y*(z-a-2*r)**0*np.heaviside(z-a-2*r, 0)
    return force

# Calculating the internal moment about x axis as a function of z
def calculate_internal_moment(z):
    moment = MA + Ay*z - P1y*(z-a)**1*np.heaviside(z-a, 0) - P2y*(z-a-2*r)**1*np.heaviside(z-a-2*r, 0)
    return moment

# Calculating the theta as a function of z, to get actual theta multiply this
# by stiffness E*I
def calculate_thetaEI(z):
    theta = MA*z + 0.5*Ay*z**2 - 0.5*P1y*(z-a)**2*np.heaviside(z-a, 0) - 0.5*P2y*(z-a-2*r)**2*np.heaviside(z-a-2*r, 0)
    return -theta

# Calculating the deflection as a function of z, to get actual deflection
# multiply this by stiffness E*I

def calculate_vEI(z):
    v = 0.5*MA*z**2 + 1/6*Ay*z**3 - 1/6*P1y*(z-a)**3*np.heaviside(z-a, 0) - 1/6*P2y*(z-a-2*r)**3*np.heaviside(z-a-2*r, 0)
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
