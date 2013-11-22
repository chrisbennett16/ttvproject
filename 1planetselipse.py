# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 13:38:23 2013

@author: chrisbennett
"""
#period in this is not accurate as we have step errors, therefore will colculate it bases on the semi major axis
#problem with acceleration over or under powered by stepping
from pylab import*
import matplotlib.pylab as plt


from math import pi, cos, sin, acos, asin, atan
G=6.67384e-11

def effaccstar(angle1,ecc1,a0):
    distance1=a0/(1+ecc1*cos(angle1))
    x1=distance1*cos(angle1)
    y1=distance1*sin(angle1)
    if sin(angle1)==0 and x1>0:
        angle2=3*pi/2
    elif sin(angle1)==0 and x1<0:
        angle2=pi/2
    else:
        angle2=atan(-(cos(angle1)+ecc1)/sin(angle1))
    if angle1>pi:
        angle2=angle2-pi
    if angle2<0:
        angle2=angle2+2*pi
    angle3=angle2-angle1+pi
    #angle4 is the angle between velocity direction and perpendicular to the radius    
    angle4=pi/2-angle3
    
    effacc=cos(angle3)/(distance1**2)
    return effacc
    

# inouts, require masses periods radii and step number for 2nd planet orbit

star_mass_in_ms = raw_input('What is the stars mass in solar masses?')
star_mass=float(star_mass_in_ms)*1.98e30
star_radius=6.995e8*float(raw_input('star radius in sun radii'))

step_number=int(raw_input("How many steps?"))

#for planet1

smAxisERadii1=raw_input("What is the SemiMajorAxis in Earth Orbits?")
smAxis1=float(smAxisERadii1)*1.49e11
period1=(smAxis1**3*4*pi**2/(star_mass*G))**0.5
velocity1=2*pi*smAxis1/period1
planet_radius1=6.371e6*float(raw_input('planet radius 1 in earth radii?'))
planet_mass1=6e24*float(raw_input('planet mass 1 in earth masses?'))
ecc1=float(raw_input("What is the eccentricity of planet 1's orbit?"))
planet_x1=0.0

planet_y1=smAxis1

plot1=[]
    

for t in range(1,int(2*step_number)):
    if (planet_x1>0 and planet_y1>=0):    
        angle1=atan(planet_y1/planet_x1)
    elif (planet_x1<0):
        angle1=atan(planet_y1/planet_x1)+pi
    elif (planet_x1==0 and planet_y1>0):
        angle1=pi/2
    else:
        angle1=atan(planet_y1/planet_x1)+2*pi
    distance1=smAxis1/(1+ecc1*cos(angle1))
    planet_a1=-effaccstar(angle1,ecc1,smAxis1)*(G*star_mass)
    velocity1=velocity1+planet_a1*period1/step_number
    angle1=angle1+velocity1*period1/(step_number*distance1)
    distance1=smAxis1/(1+ecc1*cos(angle1))    
    planet_x1=distance1*cos(angle1)
    planet_y1=distance1*sin(angle1)
    
plot(plot1)
    
