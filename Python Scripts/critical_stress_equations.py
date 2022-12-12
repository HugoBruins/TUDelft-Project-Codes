# Import modules
import math
import numpy as np

# Set variables
R = 1 # m
L = 1 # m
t_1 = 1 # m
t_2 = 1 # m
E = 1e9 #Pa
v = 1 # -
p = 1 # Pa

# Calculate the cross-sectional area moment of inertia
def moi(r, t_1):
    return math.pi/4*((r+t_1)**4-r**4)

# Calculate the cross-sectional area
def area(r, t_1):
    return math.pi*((r+t_1)**2-r**2)

# Calculate the Q factor needed for the shell buckling analysis
def Q(p, E, R, t_1):
    return p/E*(R/t_1)**2

# Calculate the mass of the fuel tank
def mass():
    return

# Calculate the critical column buckling stress
def sigma_cr_column(E, r, t_1, L):
    return (math.pi**2*E*moi(r, t_1))/(area(r,t_1)*L**2)

# Calculate the critical shell buckling stress
def sigma_cr_shell(p, E, R, t_1, v):
    return (1+0.983*(1-math.e**(-23.14*Q(p, E, R, t_1))))* (E*t_1)/(math.sqrt(3)*R*math.sqrt(1-v**2))

# Calculate the stresses on the fuel tank
def stress():
    return

# Add increments for the three dimensions for if the critical stresses are too low
t_incr = np.linspace(0,0.01,10)
R_incr = np.linspace(0,1,10)
L_incr = np.linspace(0,1,10)

# List to save all possible values
changes = []

# Change the three parameters separately if sigma_cr's are too low to find a possible increase in 1 parameter
if stress() < sigma_cr_shell(p, E, R, t_1, v) or stress() < sigma_cr_column(E, R, t_1, L):
    for t in t_incr:
        t_1_b = t_1 + t
        if stress() > sigma_cr_shell(p, E, R, t_1_b, v) and stress() > sigma_cr_column(E, R, t_1_b, L):
            changes.append(['t_1', t_1_b, mass()])
    for R_added in R_incr:
        R_b = R + R_added
        if stress() > sigma_cr_shell(p, E, R_b, t_1, v) and stress() > sigma_cr_column(E, R_b, t_1, L):
            changes.append(['R', R_b, mass()])
    for L_added in L_incr:
        L_b = L + L_added
        if stress() > sigma_cr_shell(p, E, R, t_1, v) and stress() > sigma_cr_column(E, R, t_1, L_b):
            changes.append(['L', L_b, mass()])

# Find the change with lowest mass (increase) and change the dimension to its new value
new_config = changes[changes.index(min(i[2] for i in changes))]
if new_config[0] == 't_1':
    t_1 += new_config[1]
elif new_config[0] == 'R':
    R += new_config[1]
elif new_config[0] == 'L':
    L += new_config[1]

# Calc stress again, rerun tests (while loop?)


# If sigma_cr's are too high (checked here with safety factor), possible decrease in dimensions
if stress() > sigma_cr_shell(p, E, R, t_1, v) and sigma_cr_column(E, R, t_1, L):
    SF = min(stress()/sigma_cr_shell(), stress()/sigma_cr_column())
    print(f'The safety factor with the initial dimensions is {SF.2f}')
    if SF > 1.5:
        # decrease parameters until viable SF is found (ESA ECSS standards?)
