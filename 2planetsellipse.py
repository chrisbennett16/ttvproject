# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 20:50:03 2013

@author: chrisbennett
"""
#need toadd planet interactions and transit detection
#planet masses dependent on already defined values so can retrieve from that however not representative
#travelling too fast due to step errors could correct by altering velocity or acceleration
#for a transit with 0.5 ecc the accel must be multiplied by 1.21459233575 to get the correct period, get this by minimisation
#will continue anyway to get interacting orbits in the same way assuming no path change
#add in something to detect if paths cross as it will be unstable
#for noncircular orbits the velocity for each loop is slightly under
from pylab import*
import matplotlib.pylab as plt


from math import pi, cos, sin, acos, asin, atan
G=6.67384e-11

def effaccstar(angle,ecc,a0):
    distance=a0/(1+ecc*cos(angle))
    x1=distance*cos(angle)
    y1=distance*sin(angle)
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
    effacc=-cos(angle03)/(distance**2)
    return effacc

def effaccplanet(anglea, ecca, a0, x2, y2):
    distance=a0/(1+ecca*cos(anglea))
    x1=distance*cos(anglea)
    y1=distance*sin(anglea)
    if sin(anglea)==0 and x1>0:
        angle02=3*pi/2
    elif sin(anglea)==0 and x1<0:
        angle02=pi/2
    else:
        angle02=atan(-(cos(anglea)+ecca)/sin(anglea))
    if anglea>pi:
        angle02=angle02-pi
    if angle02<0:
        angle02=angle02+2*pi
    #now we have the correct gradient angle need to compare to position gradients
    if x1==x2 and y2>y1:
        angle03=pi/2
    elif x1==x2 and y1>y2:
        angle03=3*pi/2
    angle03=atan((y2-y1)/(x2-x1))
    if y2>y1 and angle03<0:
        angle03=angle03+pi
    if y2<y1 and angle03>0:
        angle03=angle03+pi
    elif y2<y1 and angle03<0:
        angle03=angle03+ 2*pi
    angle04=angle03-angle02    
    
    effaccp=-cos(angle04)/((x1-x2)**2+(y1-y2)**2)
    return effaccp
def angle5(x1,y1,ecc,angle):
    if (x1==0 and y1>0):
        angleo=pi-atan(-(cos(angle)+ecc)/sin(angle))
    elif (x1==0 and y1<0):
        angleo=-atan(-(cos(angle)+ecc)/sin(angle))
    elif (angle==0 or angle==pi):
        angleo=0
    else:
        angleo=atan(y1/x1)+pi/2-atan(-(cos(angle)+ecc)/sin(angle))
    return angleo
    

# inouts, require masses periods radii and step number for 2nd planet orbit

star_mass_in_ms = raw_input('What is the stars mass in solar masses?')
star_mass=float(star_mass_in_ms)*1.98e30
star_radius=6.995e8*float(raw_input('star radius in sun radii'))

step_number=int(raw_input("How many steps?"))

#for planet1

period1=365*24*3600*float(raw_input("What is the period of planet 1 in years?"))
smAxis1=(G*star_mass*period1**2/(4*pi**2))**(1.0/3)

planet_radius1=6.371e6*float(raw_input('planet radius 1 in earth radii?'))
planet_mass1=6e24*float(raw_input('planet mass 1 in earth masses?'))
ecc1=float(raw_input("What is the eccentricity of planet 1's orbit?"))
latus1=smAxis1*(1-ecc1**2)
velocity1=(G*star_mass*(2/latus1-1/smAxis1))**(0.5)


period2=365*24*3600*float(raw_input("What is the period of planet 2 in years?"))
smAxis2=(G*star_mass*period2**2/(4*pi**2))**(1.0/3)

planet_radius2=6.371e6*float(raw_input('planet radius 2 in earth radii?'))
planet_mass2=6e24*float(raw_input('planet mass 2 in earth masses?'))
ecc2=float(raw_input("What is the eccentricity of planet 2's orbit?"))
latus2=smAxis2*(1-ecc2**2)
velocity2=(G*star_mass*(2/latus2-1/smAxis2))**(0.5)


planet_x2=0.0
planet_y2=latus2
planet_x1=0.0
planet_y1=latus1
plot2=[]
plot1=[]
plot3=[]
transit_on=1
transit_on2=1
s=0
transit_counter1=0
transit_counter2=0
print period1, period2

for t in range(1,int(50*step_number)):
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
    planet_a1=effaccstar(angle1,ecc1,latus1)*(G*star_mass)+effaccplanet(angle1, ecc1, latus1, planet_x2, planet_y2)*G*planet_mass2
    planet_a2=effaccstar(angle2,ecc2,latus2)*(G*star_mass)+effaccplanet(angle2, ecc2, latus2, planet_x1, planet_y1)*G*planet_mass1
    velocity1=velocity1+planet_a1*period1/step_number
    velocity2=velocity2+planet_a2*period1/step_number
    angle1=angle1+velocity1*period1*abs(cos(angle5(planet_x1,planet_y1,ecc1,angle1)))/(step_number*distance1)
    angle2=angle2+velocity2*period1*abs(cos(angle5(planet_x2,planet_y2,ecc2,angle2)))/(step_number*distance2)
    distance1=latus1/(1+ecc1*cos(angle1))
    distance2=latus2/(1+ecc2*cos(angle2))    
    planet_x1=distance1*cos(angle1)
    planet_y1=distance1*sin(angle1)
    planet_x2=distance2*cos(angle2)
    planet_y2=distance2*sin(angle2)
    
    #transit detector required
    if planet_x1<star_radius+planet_radius1 and planet_x1>-star_radius-planet_radius1 and planet_y1>0 and transit_on==1:
        transit_on=0
        transit_counter1=transit_counter1+1
        print 'Transit for planet 1 '+ str(t)
        if transit_counter1!=1 and transit_counter1!=2:
            plot1.append(t)
        elif transit_counter1==2:
            first_transit=t
        s=t
    elif transit_on==0 and not (planet_x1<star_radius+planet_radius1 and planet_x1>-star_radius-planet_radius1 and planet_y1>0):
        transit_on=1
        print 'Transit finished ' + str(t)
 
    if planet_x2<star_radius+planet_radius2 and planet_x2>-star_radius-planet_radius2 and planet_y2>0 and transit_on2==1:
        print 'To2 '+ str(t)
        transit_on2=0
        transit_counter2=transit_counter2+1
    elif transit_on2==0 and not (planet_x2<star_radius+planet_radius2 and planet_x2>-star_radius-planet_radius2 and planet_y2>0):
        transit_on2=1
        print 'Tf' +str(t-1)

plot2=[0]*len(plot1)
plot3=[0]*(len(plot1)-1)
average_transit=(float(s)-float(first_transit))/(transit_counter1-2)
print average_transit
for a in range(0,len(plot1)):
    plot2[a]=plot1[a]-average_transit*(a+2)
    if a!=0:
        plot3[a-1]=plot1[a]-plot1[a-1]
    
plot(plot3)
print plot1
    