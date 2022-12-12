# Import modules
import math
import numpy as np

# Set variables
R = 0.442  # [m]
L = 2*R  # [m]
t_1 = 0.001  # [m]
t_2 = 0.001  # [m]
E = 70e9  # [Pa]
v = 0.342  # [-]
p = 17e5  # [Pa]
# F = 700  # [N], from WP4
rho = 4400  # [kg/m3]
mass_fuel = 880 * 0.362 + 1440 * 0.359

# Calculate the cross-sectional area moment of inertia
def moi(r, t_1):
    return math.pi/4*(r**4-(r-t_1)**4)

# Calculate the cross-sectional area
def area(r, t_1):
    return math.pi*(r**2-(r-t_1)**2)

# Calculate the Q factor needed for the shell buckling analysis
def Q(p, E, R, t_1):
    return p/E*(R/t_1)**2

# Calculate the mass of the fuel tank
def mass(R, L, t_1, t_2, rho):
    return rho * ((math.pi*(R**2-(R-t_1)**2)*(L-2*R)) + 4/3 * math.pi * (R**3 - (R-t_1)**3))

# Calculate the critical column buckling stress
def sigma_cr_column(E, r, t_1, L):
    return (math.pi**2*E*moi(r, t_1))/(area(r, t_1)*L**2)

# Calculate the critical shell buckling stress
def sigma_cr_shell(p, E, R, t_1, v):
    return (1+0.983*(1-math.e**(-23.14*Q(p, E, R, t_1)))) * (E*t_1)/(math.sqrt(3)*R*math.sqrt(1-v**2))

# Calculate the stresses on the fuel tank
def stress(F, r, t_1):
    return F / area(r, t_1)


# Add increments for the three dimensions for if the critical stresses are too low
t_incr = np.linspace(0, 100, 10000)
R_incr = np.linspace(0, 10000, 10000)
L_incr = np.linspace(0, 10000, 10000)

# List to save all possible values
changes = []

# Calculate the force on the fuel tank
F = (mass(R, L, t_1, t_2, rho) + mass_fuel) * 6 * 9.80665
print(F)

# Change the three parameters separately if sigma_cr's are too low to find a possible increase in 1 parameter
if stress(F, R, t_1) > sigma_cr_shell(p, E, R, t_1, v) or stress(F, R, t_1) > sigma_cr_column(E, R, t_1, L):
    for t_added in t_incr:
        t_1_b = t_1 + t_added
        if stress(F, R, t_1_b) < sigma_cr_shell(p, E, R, t_1_b, v) and stress(F, R, t_1_b) < sigma_cr_column(E, R, t_1_b, L):
            changes.append(['t_1', t_1_b, mass(R, L, t_1_b, t_2, rho)])
            print('found one')
            break
    for R_added in R_incr:
        R_b = R + R_added
        if stress(F, R_b, t_1) < sigma_cr_shell(p, E, R_b, t_1, v) and stress(F, R_b, t_1) < sigma_cr_column(E, R_b, t_1, L):
            changes.append(['R', R_b, mass(R_b, L, t_1, t_2, rho)])
            print('found one')
            break
    for L_added in L_incr:
        L_b = L + L_added
        if stress(F, R, t_1) < sigma_cr_shell(p, E, R, t_1, v) and stress(F, R, t_1) < sigma_cr_column(E, R, t_1, L_b):
            changes.append(['L', L_b, mass(R, L_b, t_1, t_2, rho)])
            print('found one')
            break

# Find the change with the lowest mass (increase) and change the dimension to its new value
if len(changes) > 0:
    new_config = changes[changes.index(min(i[2] for i in changes))]
    if new_config[0] == 't_1':
        t_1 += new_config[1]
    elif new_config[0] == 'R':
        R += new_config[1]
    elif new_config[0] == 'L':
        L += new_config[1]

    print(new_config)
# Calc stress again, rerun tests (while loop?)


# If sigma_cr's are too high (checked here with safety factor), possible decrease in dimensions
if stress(F, R, t_1) < sigma_cr_shell(p, E, R, t_1, v) and stress(F, R, t_1) < sigma_cr_column(E, R, t_1, L):
    SF = min(sigma_cr_shell(p, E, R, t_1, v)/stress(F, R, t_1), sigma_cr_column(E, R, t_1, L)/stress(F, R, t_1))
    print(f'The safety factor with the initial dimensions is {SF}')
    if SF > 1.5: # decrease parameters until viable SF is found (ESA ECSS standards?)
        print('more than 1.5')

print(stress(F,R,t_1)/1e6, 'MPa')
print(sigma_cr_shell(p, E, R, t_1, v)/1e6, 'MPa')
print(sigma_cr_column(E, R, t_1, L)/1e6,'MPa')
