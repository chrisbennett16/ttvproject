# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 12:04:26 2013

@author: chrisbennett
"""
#making a program that can use n planets and trying to correct for the velocity problem
from pylab import*
import matplotlib.pylab as plt


from math import pi, cos, sin, acos, asin, atan
G=6.67384e-11
class Planet(object):
    def __init__(self,periodYears,radius,mass,ecc):
        self.periodYears=periodYears
        self.period=float(self.periodYears)*365.25*3600*24
        self.radius=radius
        self.mass=mass
        self.ecc=ecc
        self.smAxis=(G*starMass*self.period**2/(4*pi**2))**(1.0/3)
        self.latus=self.smAxis*(1-ecc**2)
        self.velocity=(G*starMass*(2/self.latus-1/self.smAxis))**(0.5)
        self.x=0.0
        self.y=self.latus
        self.angle=0.0
        self.distance=self.latus/(1+ecc*cos(self.angle))
        self.da=0.0
        self.dv=0.0
        self.transit=1
        self.lastTransit=0
        self.transitCounter=0
        self.vchange=0.0
        self.orbitCounter=0
        self.orbitTest=1
        
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
    

#star inputs
starMassInMs = raw_input('What is the stars mass in solar masses?')
starMass=float(starMassInMs)*1.98e30
starRadius=6.995e8*float(raw_input('star radius in sun radii'))

step_number=int(raw_input("How many steps?"))
planet=[0]
planet_number=int(raw_input('How many planets?'))

#adding planets

for a in range(0,planet_number):
    print 'Planet ' + str(a+1)
    planet.append(Planet(float(raw_input("What is the period in years?")),6.371e6*float(raw_input('radius in earth radii?')),6e24*float(raw_input('planet mass in earth masses?')),float(raw_input("What is the eccentricity?"))))
#will sort transit detectors later and plotting
plot1=[]
plot2=[]
plot3=[]
plot4=[]
plot5=[]
plot6=[]
plot7=[]
plot8=[]
plot9=[]
plot10=[]
plot11=[]
for a in range(1,planet_number+1):
    planet[a].vchange=0.0
#code to find average dv
testSteps=1
for a in range(1,planet_number+1):
    testSteps=testSteps*int(planet[a].periodYears)
precision=1000
print testSteps
for repeats in range(0,5):
    for a in range(1,planet_number+1):
        planet[a].x=0.0
        planet[a].y=planet[a].latus
        planet[a].dv=0.0
        planet[a].orbitCounter=0
        planet[a].orbitTest=1
        del plot3[:]
        del plot4[:]
    for t in range(1,(testSteps*10+1)*precision):
        for a in range(1,planet_number+1):
            planet[a].da=0.0
            if (planet[a].x>0 and planet[a].y>=0):    
                planet[a].angle=atan(planet[a].y/planet[a].x)
            elif (planet[a].x<0):
                planet[a].angle=atan(planet[a].y/planet[a].x)+pi
            elif (planet[a].x==0 and planet[a].y>0):
                planet[a].angle=pi/2
            elif (planet[a].x==0 and planet[a].y<0):
                planet[a].angle=3*pi/2
            else:
                planet[a].angle=atan(planet[a].y/planet[a].x)+2*pi     
            planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle))
            for b in range (1,planet_number+1):
                if b!=a:
                    # toavoid calculating the acc from itself
                    planet[a].da=planet[a].da+effaccplanet(planet[a].angle, planet[a].ecc, planet[a].latus, planet[b].x, planet[b].y)*G*planet[b].mass    
            planet[a].dv=planet[a].dv+planet[a].da*planet[1].period/precision
            planet[a].velocity=(G*starMass*(2/planet[a].distance-1/planet[a].smAxis))**(0.5)+planet[a].dv-planet[a].vchange    
            planet[a].angle=planet[a].angle+planet[a].velocity*planet[1].period*abs(cos(angle5(planet[a].x,planet[a].y,planet[a].ecc,planet[a].angle)))/(precision*planet[a].distance)
            planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle))
        for a in range(1,planet_number+1):
            planet[a].x=planet[a].distance*cos(planet[a].angle)
            planet[a].y=planet[a].distance*sin(planet[a].angle)
            if planet[a].x<0.0 and planet[a].orbitTest==0:
                planet[a].orbitCounter=planet[a].orbitCounter+1
                planet[a].orbitTest=1
            elif planet[a].orbitTest==1 and planet[a].x>0.0:
                planet[a].orbitTest=0
        if repeats==0:
            plot7.append(planet[1].x)
        if planet[1].orbitCounter==(10*testSteps):
            break
        plot3.append(planet[1].dv)
        plot4.append(planet[2].dv)
        
    planet[1].vchange=sum(plot3)/len(plot3)
    print planet[1].vchange
    planet[2].vchange=sum(plot4)/len(plot4)
    print planet[2].vchange




for a in range(1,planet_number+1):
    planet[a].x=0.0
    planet[a].y=planet[a].latus
    planet[a].dv=0.0


for t in range(1,300*step_number):
    for a in range(1,planet_number+1):
        planet[a].da=0.0
        if (planet[a].x>0 and planet[a].y>=0):    
            planet[a].angle=atan(planet[a].y/planet[a].x)
        elif (planet[a].x<0):
            planet[a].angle=atan(planet[a].y/planet[a].x)+pi
        elif (planet[a].x==0 and planet[a].y>0):
            planet[a].angle=pi/2
        elif (planet[a].x==0 and planet[a].y<0):
            planet[a].angle=3*pi/2
        else:
            planet[a].angle=atan(planet[a].y/planet[a].x)+2*pi     
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle))
        for b in range (1,planet_number+1):
            if b!=a:
                # toavoid calculating the acc from itself
                planet[a].da=planet[a].da+effaccplanet(planet[a].angle, planet[a].ecc, planet[a].latus, planet[b].x, planet[b].y)*G*planet[b].mass    
        planet[a].dv=planet[a].dv+planet[a].da*planet[1].period/step_number
        planet[a].velocity=(G*starMass*(2/planet[a].distance-1/planet[a].smAxis))**(0.5)+planet[a].dv-planet[a].vchange   
        planet[a].angle=planet[a].angle+planet[a].velocity*planet[1].period*abs(cos(angle5(planet[a].x,planet[a].y,planet[a].ecc,planet[a].angle)))/(step_number*planet[a].distance)
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle))
    for a in range(1,planet_number+1):
        planet[a].x=planet[a].distance*cos(planet[a].angle)
        planet[a].y=planet[a].distance*sin(planet[a].angle)
    #finding the transits of the planets
    for a in range(1,planet_number+1):
        if planet[a].x<starRadius+planet[a].radius and planet[a].x>-starRadius-planet[a].radius and planet[a].y>0 and planet[a].transit==0:
            planet[a].transitCounter=planet[a].transitCounter+1            
            planet[a].transit=1            
            print 'Transit start for planet ' + str(a) + ' at ' + str(t)
            if a==1:
                plot1.append(t)
            planet[a].lastTransit=t
        elif planet[a].transit==1 and not (planet[a].x<starRadius+planet[a].radius and planet[a].x>-starRadius-planet[a].radius and planet[a].y>0):
            planet[a].transit=0
            print 'Transit end for planet ' +str(a) + ' at ' + str(t)
    plot10.append(planet[1].velocity)
averagePeriod=float((planet[1].lastTransit-plot1[0]))/(planet[1].transitCounter-1)
print averagePeriod
for a in range (0, len(plot1)):
    plot2.append(plot1[a]-(step_number*(a+1)))
plot(plot10)
    
print plot1[0]