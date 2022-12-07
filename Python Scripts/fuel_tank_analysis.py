import math
import numpy as np

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
        pass
    

class BeamAttachment():
    def __init__(self, material):
        # Calculate Ixx
        
        # Calculate spring constant lateral
        
        # Calculate spring constant axial
        pass
        
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
        pass\
            
hydrazine_fuel_tank = Fueltank()

hydrazine_fuel_tank.mass

AL6061T6 = Material(68e9, 276e6, 2700, 0.33)    # aluminium alloy
AL7075T6 = Material(71e9, 503e6, 2800, 0.33)    # aluminium alloy
TI6ALV4  = Material(110e9, 825e6, 4400, 0.33)   # commonly used Titanium alloy
MGA231B  = Material(45e9, 220e6, 1700, 0.33)    # magnesium alloy
BES65A   = Material(304e9, 207e6, 2000, 0.33)   # berrylium alloy (very non toxic!)
AM350    = Material(200e9, 1034e6, 7700, 0.33)  # ferrous alloy

