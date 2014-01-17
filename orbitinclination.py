# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 08:55:40 2014

@author: chrisbennett
"""

#making a program that can use n planets and trying to correct for the velocity problem
from pylab import*
import matplotlib.pylab as plt


from math import pi, cos, sin, acos, asin, atan
G=6.67384e-11
class Planet(object):
    def __init__(self,periodYears,radius,mass,ecc,inclination,azimuth):
        self.periodYears=periodYears
        self.period=float(self.periodYears)*365.25*3600*24
        self.radius=radius
        self.mass=mass
        self.ecc=ecc
        self.inclination=inclination
        self.azimuth=azimuth
        self.smAxis=(G*starMass*self.period**2/(4*pi**2))**(1.0/3)
        self.latus=self.smAxis*(1-ecc**2)
        self.velocity=(G*starMass*(2/self.latus-1/self.smAxis))**(0.5)
        self.x=0.0
        self.y=0.0
        self.angle=pi/2
        self.distance=self.latus/(1+ecc*cos(self.angle+self.azimuth))
        self.da=0.0
        self.dv=0.0
        self.transit=1
        self.lastTransit=0
        self.transitCounter=0
        self.vchange=0.0
        self.orbitCounter=0
        self.orbitTest=1
        self.accTotal=0.0
        self.accMean=0.0
        self.dt=0.0
        self.dvlist=[]
        self.transitTimes=[]
        self.averagePeriod=0.0
        self.ominusc=[]
def effaccplanet(anglea, ecca, a0,azimuth, x2, y2):
    distance=a0/(1+ecca*cos(anglea+azimuth))
    x1=distance*cos(anglea)
    y1=distance*sin(anglea)
    if ecca*(cos(anglea)*sin(anglea+azimuth)-sin(anglea)*cos(anglea+azimuth))-sin(anglea)==0 and x1>0:
        angle02=pi/2
        print 'hi'
    elif ecca*(cos(anglea)*sin(anglea+azimuth)-sin(anglea)*cos(anglea+azimuth))-sin(anglea)==0 and x1<0:
        angle02=3*pi/2
        print 'ho'
    else:
        angle02=atan((ecca*(sin(anglea)*sin(anglea+azimuth)+cos(anglea)*cos(anglea+azimuth))+cos(anglea))/(ecca*(cos(anglea)*sin(anglea+azimuth)-sin(anglea)*cos(anglea+azimuth))-sin(anglea)))
    if ecca*(cos(anglea)*sin(anglea+azimuth)-sin(anglea)*cos(anglea+azimuth))-sin(anglea)>0:
        angle02=angle02+pi

    if angle02<0:
        angle02=angle02+2*pi
    #now we have the correct gradient angle need to compare to position gradients
    if x1==x2 and y2>y1:
        angle03=3*pi/2
    elif x1==x2 and y1>y2:
        angle03=pi/2
    else:
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
def angle5(ecc,angle):
    x1=cos(angle)
    y1=sin(angle)
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
    planet.append(Planet(float(raw_input("What is the period in years?")),6.371e6*float(raw_input('radius in earth radii?')),6e24*float(raw_input('planet mass in earth masses?')),float(raw_input("What is the eccentricity?")),2*pi*float(raw_input('inlination in degrees?(between ±90)'))/360,2*pi*float(raw_input('azimuthal angle?(in degrees)'))/360))
#will sort transit detectors later and plotting
plot1=[]
plot2=[]
#code to find average dv
testSteps=1
for a in range(1,planet_number+1):
    testSteps=testSteps*int(planet[a].periodYears)
precision=1000
repeat=5
print testSteps
for repeats in range(0,repeat):
    for a in range(1,planet_number+1):
        planet[a].angle=pi/2
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azimuth))
        planet[a].x=planet[a].distance*cos(planet[a].angle)
        planet[a].y=planet[a].distance*sin(planet[a].angle)
        planet[a].dv=0.0
        planet[a].orbitCounter=0
        planet[a].orbitTest=1
        del planet[a].dvlist[:]
    for t in range(1,int(testSteps+1)*precision):
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
            planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azimuth))
            for b in range (1,planet_number+1):
                if b!=a:
                    # toavoid calculating the acc from itself
                    planet[a].da=planet[a].da+effaccplanet(planet[a].angle, planet[a].ecc, planet[a].latus,planet[a].azimuth, planet[b].x, planet[b].y)*G*planet[b].mass              
            if repeats==(repeat-1):
                planet[a].accTotal=planet[a].accTotal+abs(planet[a].da)
            planet[a].dv=planet[a].dv+planet[a].da*planet[1].period/precision
            planet[a].velocity=(G*starMass*(2/planet[a].distance-1/planet[a].smAxis))**(0.5)+planet[a].dv-planet[a].vchange       
            planet[a].angle=planet[a].angle+planet[a].velocity*planet[1].period*abs(cos(angle5(planet[a].ecc,planet[a].angle+planet[a].azimuth)))/(precision*planet[a].distance)
            planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azimuth))
        for a in range(1,planet_number+1):
            planet[a].x=planet[a].distance*cos(planet[a].angle)
            planet[a].y=planet[a].distance*sin(planet[a].angle)
            if planet[a].x<0.0 and planet[a].orbitTest==0:
                planet[a].orbitCounter=planet[a].orbitCounter+1
                planet[a].orbitTest=1
            elif planet[a].orbitTest==1 and planet[a].x>0.0:
                planet[a].orbitTest=0
        if planet[1].orbitCounter==int(testSteps):
            break
        for a in range(1,planet_number+1):
            planet[a].dvlist.append(planet[a].dv)
    for a in range(1,planet_number+1):
        planet[a].vchange=sum(planet[a].dvlist)/len(planet[a].dvlist)
        print a,planet[a].vchange


for a in range(1,planet_number+1):
    planet[a].angle=pi/2
    planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azimuth))
    planet[a].x=planet[a].distance*cos(planet[a].angle)
    planet[a].y=planet[a].distance*sin(planet[a].angle)
    planet[a].dv=0.0
    planet[a].accMean=planet[a].accTotal/t
dt=0.0
t=0.0
stepCounter=0
while t<(planet[1].period*2):
    stepCounter=stepCounter+1
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
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azimuth))
        for b in range (1,planet_number+1):
            if b!=a:
                # toavoid calculating the acc from itself
                planet[a].da=planet[a].da+effaccplanet(planet[a].angle, planet[a].ecc, planet[a].latus,planet[a].azimuth, planet[b].x, planet[b].y)*G*planet[b].mass    
        if planet[a].da==0:
            planet[a].dt=planet[1].period*10/(step_number)
        else:
            planet[a].dt=planet[a].accMean*planet[1].period/(step_number*abs(planet[a].da))
    dt=planet[1].dt
    for a in range(1,planet_number+1):
        if planet[a].dt<=dt:
            if planet[a].dt<planet[1].period*10/(step_number):
                dt=planet[a].dt
            else:
                dt=planet[1].period*10/(step_number)
        if planet[a].x<2*(starRadius+planet[a].radius) and planet[a].x>2*(-starRadius-planet[a].radius) and planet[a].y>0:
            #for 3 seconds precision stepnumber*1000 is 10000000            
            dt=planet[1].period/(10000000)
        elif planet[a].x<10*(starRadius+planet[a].radius) and planet[a].x>10*(-starRadius-planet[a].radius) and planet[a].y>0:
            #for 3 seconds precision stepnumber*1000 is 100000000            
            dt=planet[1].period/(1000000)
        
    t=t+dt
    
    for a in range(1,planet_number+1):
        planet[a].dv=planet[a].dv+planet[a].da*dt
        planet[a].velocity=(G*starMass*(2/planet[a].distance-1/planet[a].smAxis))**(0.5)+planet[a].dv-planet[a].vchange   
        planet[a].angle=planet[a].angle+planet[a].velocity*dt*abs(cos(angle5(planet[a].ecc,planet[a].angle+planet[a].azimuth)))/(planet[a].distance)
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azimuth))
    for a in range(1,planet_number+1):
        planet[a].x=planet[a].distance*cos(planet[a].angle)
        planet[a].y=planet[a].distance*sin(planet[a].angle)
    #finding the transits of the planets
    for a in range(1,planet_number+1):
        if planet[a].x<starRadius+planet[a].radius and planet[a].x>-starRadius-planet[a].radius and planet[a].y>0 and planet[a].transit==0:
            planet[a].transitCounter=planet[a].transitCounter+1            
            planet[a].transit=1            
            print 'Transit start for planet ' + str(a) + ' at ' + str(t) + ' with duration' + str(planet[a].velocity*starRadius+planet[a].radius*2)
            planet[a].transitTimes.append(t)            
            planet[a].lastTransit=t
        elif planet[a].transit==1 and not (planet[a].x<starRadius+planet[a].radius and planet[a].x>-starRadius-planet[a].radius and planet[a].y>0):
            planet[a].transit=0
    plot1.append(-planet[1].y)
    plot2.append(t)
for a in range(1,planet_number +1):
    if planet[a].transitCounter>1:
        planet[a].averagePeriod=float((planet[a].lastTransit-planet[a].transitTimes[0]))/(planet[a].transitCounter-1)
    else:
        print 'period is not calculated for planet ' + str(a)
        averagePeriod=planet[1].period
    for b in range (0, len(planet[a].transitTimes)):
        planet[a].ominusc.append(planet[a].transitTimes[b]-(planet[a].averagePeriod*(b+1)))
        
    print planet[a].averagePeriod-planet[a].period
plot(plot2,plot1)
    
print stepCounter