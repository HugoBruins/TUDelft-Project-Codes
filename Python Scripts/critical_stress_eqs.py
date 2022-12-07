import math
# import numpy as np

def moi(r, t_1):
    return math.pi/4*((r+t_1)**4-r**4)

def area(r, t_1):
    return math.pi*((r+t_1)**2-r**2)

def poisson(p, E, R, t_1):
    return(p/E)*(R/t_1)**2

def Q(p, E, R, t_1):
    return p/E*(R/t_1)**2

def sigma_cr_column(E, r, t_1, L):
    return (math.pi**2*E*moi(r, t_1))/(area(r,t_1)*L**2)

def sigma_cr_shell(p, E, R, t_1):
    return (1+0.983*(1-math.e**(-23.14*Q(p, E, R, t_1))))* (E*t_1)/(math.sqrt(3)*R*math.sqrt(1-poisson(p, E, R, t_1)**2))
