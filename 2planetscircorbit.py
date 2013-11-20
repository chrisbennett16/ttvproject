# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 22:12:30 2013

@author: chrisbennett
"""

#for two planets and a fixed sun in a circular orbit, transit times
#using a small shift, need to use the fact that forces affect accels therefore different velocities
#find position by previous position plus velocity times time
#might need to correct as it is small shifts
#have to extend to noncircular orbits


from math import pi, cos, sin, acos, asin
G=6.67384e-11

#this works out the acceleration in the direction of travel not including GM

def effacc(x1,y1,x2,y2,r1,r2):
    if y1>0:
        angle1=acos(x1/r1)
    else:
        angle1=2*pi-acos(x1/r1)
    if y2>0:
        angle2=acos(x2/r2)
    else:
        angle2=2*pi-acos(x2/r2)
    if angle2-angle1>0 and angle2-angle1<pi or angle2-angle1<-pi:
        distance=-r1*(1-((x1*x2+y1*y2)/((x1**2+y1**2)*(x2**2+y2**2)))**2)**(0.5)
    else:
        distance=r1*(1-((x1*x2+y1*y2)/((x1**2+y1**2)*(x2**2+y2**2)))**2)**(0.5)
    a=distance/((x1-x2)**2+(y1-y2)**2)**(1.5)
    return a


# inouts, require masses periods radii and step number for 2nd planet orbit

star_mass_in_ms = raw_input('What is the stars mass in solar masses?')
star_mass=float(star_mass_in_ms)*1.98e30
star_radius=6.995e8*float(raw_input('star radius in sun radii'))

step_number=int(raw_input("How many steps?"))
#for planet1
period_years1=raw_input("What is the period of planet 1 in years?")
period1=float(period_years1)*3600*24*365.25
orbital_radius1=(period1**2*star_mass*G/(4*pi**2))**(1.0/3)
velocity1=2*pi*orbital_radius1/period1
planet_radius1=6.371e6*float(raw_input('planet radius 1 in earth radii'))
planet_mass1=6e24*float(raw_input('planet mass 1 in earth masses'))
#for planet 2
period_years2=raw_input("What is the period of planet 2 in years?")
period2=float(period_years2)*3600*24*365.25
orbital_radius2=(period2**2*star_mass*G/(4*pi**2))**(1.0/3)
velocity2=2*pi*orbital_radius2/period2
planet_radius2=6.371e6*float(raw_input('planet radius 2 in earth radii'))
planet_mass2=6e24*float(raw_input('planet mass 2 in earth masses'))

#planets start aligned

planet_x1=orbital_radius1
planet_y1=0.0
planet_x2=orbital_radius2
planet_y2=0.0

#to detect a change for transit state

transit_on=1
transit_on2=1

#time goes for 10 periods of the second planet

for t in range(1,10*int(period_years2)*step_number):

#calculating new accelerations and velocities
    
    planet_a2=effacc(planet_x1,planet_y1,planet_x2,planet_y2,orbital_radius1,orbital_radius2)*(G*planet_mass1)
    planet_a1=effacc(planet_x2,planet_y2,planet_x1,planet_y1,orbital_radius2,orbital_radius1)*(G*planet_mass2)
    velocity1=velocity1+planet_a1*period2/step_number
    velocity2=velocity2+planet_a2*period2/step_number


#move position on circle around by the new velocity times time

    if planet_y1<0:
        planet_x1a=orbital_radius1*cos(velocity1*period2/(step_number*orbital_radius1)+2*pi-acos(planet_x1/orbital_radius1))
    else:
        planet_x1a=orbital_radius1*cos(velocity1*period2/(step_number*orbital_radius1)+acos(planet_x1/orbital_radius1))
    if planet_x1<0:
        planet_y1=orbital_radius1*sin(velocity1*period2/(step_number*orbital_radius1)+pi-asin(planet_y1/orbital_radius1))
    else:
        planet_y1=orbital_radius1*sin(velocity1*period2/(step_number*orbital_radius1)+asin(planet_y1/orbital_radius1))
    planet_x1=planet_x1a

    if planet_y2<0:
        planet_x2a=orbital_radius2*cos(velocity2*period2/(step_number*orbital_radius2)+2*pi-acos(planet_x2/orbital_radius2))
    else:
        planet_x2a=orbital_radius2*cos(velocity2*period2/(step_number*orbital_radius2)+acos(planet_x2/orbital_radius2))
    if planet_x2<0:
        planet_y2=orbital_radius2*sin(velocity2*period2/(step_number*orbital_radius2)+pi-asin(planet_y2/orbital_radius2))
    else:
        planet_y2=orbital_radius2*sin(velocity2*period2/(step_number*orbital_radius2)+asin(planet_y2/orbital_radius2))
    planet_x2=planet_x2a
 
#transit detecting    

    if planet_x1<star_radius+planet_radius1 and planet_x1>-star_radius-planet_radius1 and planet_y1>0 and transit_on==1:
        transit_on=0
        print 'Transit for planet 1 '+ str(t)
    elif transit_on==0 and not (planet_x1<star_radius+planet_radius1 and planet_x1>-star_radius-planet_radius1 and planet_y1>0):
        transit_on=1
        print 'Transit finished ' + str(t-1)
 
    if planet_x2<star_radius+planet_radius2 and planet_x2>-star_radius-planet_radius2 and planet_y2>0 and transit_on2==1:
        print 'To2 '+ str(t)
        transit_on2=0
    elif transit_on2==0 and not (planet_x2<star_radius+planet_radius2 and planet_x2>-star_radius-planet_radius2 and planet_y2>0):
        transit_on2=1
        print 'Tf' +str(t-1)