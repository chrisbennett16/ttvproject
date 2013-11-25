# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 20:50:03 2013

@author: chrisbennett
"""
#need toadd planet interactions and transit detection
#planet masses dependent on already defined values so can retrieve from that however not representative
#travelling too fast due to step errors could correct by altering velocity or acceleration
#for a transit with 0.5 ecc the accel must be multiplied by 1.21459233575 to get the correct period, get this by minimisation
#will continue anyway to get interacting orbits in the same way assumin no path change
from pylab import*
import matplotlib.pylab as plt


from math import pi, cos, sin, acos, asin, atan
G=6.67384e-11

def effaccstar(angle,ecc,a0):
    distance=a0/(1+ecc*cos(angle))
    x1=distance*cos(angle1)
    y1=distance*sin(angle1)
    if sin(angle)==0 and x1>0:
        angle02=3*pi/2
    elif sin(angle)==0 and x1<0:
        angle02=pi/2
    else:
        angle02=atan(-(cos(angle)+ecc)/sin(angle))
    if angle>pi:
        angle02=angle02-pi
    if angle02<0:
        angle02=angle02+2*pi
    angle03=angle02-angle+pi
    #angle4 is the angle between velocity direction and perpendicular to the radius    
    angle04=pi/2-angle03
    
    effacc=cos(angle03)/(distance**2)
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

planet_radius1=6.371e6*float(raw_input('planet radius 1 in earth radii?'))
planet_mass1=6e24*float(raw_input('planet mass 1 in earth masses?'))
ecc1=float(raw_input("What is the eccentricity of planet 1's orbit?"))
velocity1=(G*star_mass*(2/(smAxis1*(1-ecc1**2))-1/smAxis1))**(0.5)


smAxisERadii2=raw_input("What is the SemiMajorAxis of planet 2 in Earth Orbits?")
smAxis2=float(smAxisERadii2)*1.49e11
period2=(smAxis2**3*4*pi**2/(star_mass*G))**0.5

planet_radius2=6.371e6*float(raw_input('planet radius 2 in earth radii?'))
planet_mass2=6e24*float(raw_input('planet mass 2 in earth masses?'))
ecc2=float(raw_input("What is the eccentricity of planet 2's orbit?"))
velocity2=(G*star_mass*(2/(smAxis2*(1-ecc2**2))-1/smAxis2))**(0.5)


planet_x2=0.0
latus2=smAxis2*(1-ecc2**2)
planet_y2=latus2
latus1=smAxis1*(1-ecc1**2)
planet_x1=0.0
planet_y1=latus1
plot2=[]
plot1=[]
print velocity1
    

for t in range(1,int(2*step_number)):
    if (planet_x1>0 and planet_y1>=0):    
        angle1=atan(planet_y1/planet_x1)
    elif (planet_x1<0):
        angle1=atan(planet_y1/planet_x1)+pi
    elif (planet_x1==0 and planet_y1>0):
        angle1=pi/2
    else:
        angle1=atan(planet_y1/planet_x1)+2*pi
    if (planet_x2>0 and planet_y2>=0):    
        angle2=atan(planet_y2/planet_x2)
    elif (planet_x2<0):
        angle2=atan(planet_y2/planet_x2)+pi
    elif (planet_x2==0 and planet_y2>0):
        angle2=pi/2
    else:
        angle2=atan(planet_y2/planet_x2)+2*pi
    distance1=latus1/(1+ecc1*cos(angle1))
    distance2=latus2/(1+ecc2*cos(angle2)) 
    planet_a1=-effaccstar(angle1,ecc1,latus1)*(G*star_mass)
    planet_a2=-effaccstar(angle2,ecc2,latus2)*(G*star_mass)
    velocity1=velocity1+planet_a1*period1/step_number
    velocity2=velocity2+planet_a2*period1/step_number
    angle1=angle1+velocity1*period1/(step_number*distance1)
    angle2=angle2+velocity2*period1/(step_number*distance2)
    distance1=latus1/(1+ecc1*cos(angle1))
    distance2=latus2/(1+ecc2*cos(angle2))    
    planet_x1=distance1*cos(angle1)
    planet_y1=distance1*sin(angle1)
    planet_x2=distance2*cos(angle2)
    planet_y2=distance2*sin(angle2)    
    velocity3=(G*star_mass*(2/distance1-1/smAxis1))**(0.5)
    if t==1:
        print planet_a1*period1/step_number
        print velocity1
        print velocity3
    if t==100000:
        print planet_x1
    plot2.append(velocity1)
    plot1.append(velocity3)
plot(plot1)
plot(plot2)
    