import math
import numpy as np

LF_Vol=0.362 #m^3
Ox_Vol=0.359 #m^3
SC_width=2 #m
SC_depth=2 #m
SC_height=3 #m
LF_mass=303.2 #kg
Ox_mass= 500.3 #kg
tank_pressure=17E5 #Pa

class Material():
    def __init__(self, E, sigma_y, density, v):
        # Save all the material properties in a container
        self.E = E
        self.sigma_y = sigma_y
        self.density = density
        self.v = v

class Fueltank():
    def __init__(self, p, volume, material):
        # Get the mass
        self.radius=(3*V/(4*math.pi))**(1/3)
        self.material=material
        self.thickness=(p*self.radius*SF)/(2*self.material.sigma_y)
        self.area=4*math.pi*self.radius**2
        self.vol=self.area*self.thickness
        self.mass=self.vol*self.material.density
        self.mass_tot=self.mass+Fuel_mass
        pass
    
  
# For now it is just a thin walled square for easy
class Beam():
    def __init__(self, b, t, L, material):
        # Calculate Ixx
        self.Ixx = 2/3 * b**3 * t
        self.A = 4*b*t
        
        # lateral spring constant
        self.klat = 3*material.E*self.Ixx/L**3
        
        # axial spring constant
        self.kaxial = self.A*material.E / L
        
class FuelTankAssembly():
    def __init__(self, material):
        pass
    
    def add_attachment(self):
        pass
    
    def add_fuel_tank(self):
        pass
    
    def get_natural_frequency(self):
        pass

    def get_total_mass(self):
        pass
            
hydrazine_fuel_tank = Fueltank()

hydrazine_fuel_tank.mass

AL6061T6 = Material(68e9, 276e6, 2700, 0.33)    # aluminium alloy
AL7075T6 = Material(71e9, 503e6, 2800, 0.33)    # aluminium alloy
TI6ALV4  = Material(110e9, 825e6, 4400, 0.33)   # commonly used Titanium alloy
MGA231B  = Material(45e9, 220e6, 1700, 0.33)    # magnesium alloy
BES65A   = Material(304e9, 207e6, 2000, 0.33)   # berrylium alloy (very non toxic!)
AM350    = Material(200e9, 1034e6, 7700, 0.33)  # ferrous alloy

LF_tank=Fueltank(tank_pressure,LF_Vol,AL6061T6,1.5,LF_mass)
Ox_tank=Fueltank(tank_pressure,Ox_Vol,AL6061T6,1.5,Ox_mass)

# a beam
square_beam = Beam(0.1, 1e-3, 1, AL6061T6)
print(square_beam.Ixx)
print(square_beam.klat)
