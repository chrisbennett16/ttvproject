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
        self.azproj=atan(tan(self.azimuth)*cos(self.inclination))
        if self.azimuth>pi/2 and self.azimuth<=3*pi/2:
            self.azproj=self.azproj+pi
        elif self.azimuth>3*pi/2:
            self.azproj=self.azproj+2*pi
        elif azimuth>=-3*pi/2 and azimuth<-pi/2:
            self.azproj=self.azproj-pi
        elif azimuth<-3*pi/2:
            self.azproj=self.azproj-2*pi
        self.smAxis=(G*starMass*self.period**2/(4*pi**2))**(1.0/3)
        self.latus=self.smAxis*(1-ecc**2)
        self.velocity=(G*starMass*(2/self.latus-1/self.smAxis))**(0.5)
        self.x=0.0
        self.y=0.0
        self.z=0.0
        self.angle=pi/2
        self.distance=self.latus/(1+ecc*cos(self.angle+self.azproj))
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
        
def effaccplanet(angle, ecc, distance,i,azproj,azimuth,latus,x1,y1,z1, x2, y2,z2):
    dx=distance*((cos(angle+azproj)*cos(i)*cos(azimuth)+sin(angle+azproj)*sin(azimuth))*distance*ecc*sin(angle+azproj)/latus-cos(i)*cos(azimuth)*sin(angle+azproj)+cos(angle+azproj)*sin(azimuth))
    dy=distance*((-cos(angle+azproj)*cos(i)*sin(azimuth)+sin(angle+azproj)*cos(azimuth))*distance*ecc*sin(angle+azproj)/latus+cos(i)*sin(azimuth)*sin(angle+azproj)+cos(angle+azproj)*cos(azimuth))
    dz=-distance*sin(i)*(distance*ecc*cos(angle+azproj)*sin(angle+azproj)/latus-sin(angle+azproj))
    r=((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**(0.5)    
    effaccplanet=((x2-x1)*dx+(y2-y1)*dy+(z2-z1)*dz)/((dx**2+dy**2+dz**2)**(0.5)*r**3)
    return effaccplanet
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
def transitStart(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,precision,angle1):
    solution=0.0
    distance=((latus/(1+ecc*cos(angle1)))*(cos(angle1)*cos(inc)*cos(azimuth)+sin(angle1)*sin(azimuth))-xs)**2+((latus/(1+ecc*cos(angle1)))*sin(inc)*cos(angle1)-zs)**2-(rp+rs)**2
    if distance==abs(distance) or (latus/(1+ecc*cos(angle1-azproj)))*(-cos(angle1-azproj)*cos(inc)*sin(azimuth)+sin(angle1-azproj)*sin(azimuth))<0:
        start=1
    else:
        start=0
    for a in range(0,360):
        angle = a*2*pi/(precision*360)+angle1    
        distance=((latus/(1+ecc*cos(angle)))*(cos(angle)*cos(inc)*cos(azimuth)+sin(angle)*sin(azimuth))-xs)**2+((latus/(1+ecc*cos(angle)))*sin(inc)*cos(angle)-zs)**2-(rp+rs)**2
        if distance==abs(distance) or (latus/(1+ecc*cos(angle1-azproj)))*(-cos(angle1-azproj)*cos(inc)*sin(azimuth)+sin(angle1-azproj)*sin(azimuth))<0:
            start=1
        if distance!=abs(distance) and start==1:
            if precision<10000000000:
                solution=transitStart(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,precision*10,angle-2*pi/(precision*360))
            else:
                return angle1
            break
    if solution==0.0:
        return 0.0
    else:
        return solution
def transitEnd(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,precision,angle1):
    solution=0.0    
    distance=((latus/(1+ecc*cos(angle1)))*(cos(angle1)*cos(inc)*cos(azimuth)+sin(angle1)*sin(azimuth))-xs)**2+((latus/(1+ecc*cos(angle1)))*sin(inc)*cos(angle1)-zs)**2-(rp+rs)**2
    if distance==abs(distance) or (latus/(1+ecc*cos(angle1-azproj)))*(-cos(angle1-azproj)*cos(inc)*sin(azimuth)+sin(angle1-azproj)*sin(azimuth))<0:
        start=1
    else:
        start=0
    for a in range(0,360):
        angle = a*2*pi/(precision*360)+angle1    
        distance=((latus/(1+ecc*cos(angle)))*(cos(angle)*cos(inc)*cos(azimuth)+sin(angle)*sin(azimuth))-xs)**2+((latus/(1+ecc*cos(angle)))*sin(inc)*cos(angle)-zs)**2-(rp+rs)**2
        if distance!=abs(distance):
            start=0
        if (distance==abs(distance) or (latus/(1+ecc*cos(angle1-azproj)))*(-cos(angle1-azproj)*cos(inc)*sin(azimuth)+sin(angle1-azproj)*sin(azimuth))<0) and start==0:
            if precision<10000000000:
                solution=transitEnd(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,precision*10,angle-2*pi/(precision*360))
            else:
                return angle1
            break
    if solution==0.0:
        return 0.0
    else:
        return solution
def transitLength(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs):
    angleEnd=transitEnd(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,1,0)
    angleStart=transitStart(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,1,0)
    length=(angleEnd-angleStart)*(latus/(1+ecc*cos(angleEnd-azproj)))/abs(cos(angle5(ecc,angleEnd-azproj)))
    print angleEnd, angleStart    
    return length

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
    planet.append(Planet(float(raw_input("What is the period in years?")),6.371e6*float(raw_input('radius in earth radii?')),6e24*float(raw_input('planet mass in earth masses?')),float(raw_input("What is the eccentricity?")),2*pi*float(raw_input('inclination in degrees?(between ±90)'))/360,2*pi*float(raw_input('azimuthal angle?(in degrees)'))/360))
#will sort transit detectors later and plotting
plot1=[]
plot2=[]
#code to find average dv
testSteps=1
for a in range(1,planet_number+1):
    testSteps=testSteps*int(planet[a].periodYears)
precision=1000
repeat=2
print testSteps
for repeats in range(0,repeat):
    for a in range(1,planet_number+1):
        planet[a].angle=pi/2
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
        planet[a].x=planet[a].distance*(cos(planet[a].angle+planet[a].azproj)*cos(planet[a].inclination)*cos(planet[a].azimuth)+sin(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth))
        planet[a].y=planet[a].distance*(sin(planet[a].angle+planet[a].azproj)*cos(planet[a].azimuth)-cos(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth)*cos(planet[a].inclination))
        planet[a].z=planet[a].distance*sin(planet[a].inclination)*cos(planet[a].angle+planet[a].azproj)
        if ((planet[a].angle+planet[a].azproj)%(2*pi))>pi/2 or ((planet[a].angle+planet[a].azproj)%(2*pi))<3*pi/2:
            planet[a].z=-planet[a].z
        planet[a].dv=0.0
        planet[a].orbitCounter=0
        planet[a].orbitTest=1
        del planet[a].dvlist[:]
    for t in range(1,int(testSteps+1)*precision):
        for a in range(1,planet_number+1):
            planet[a].da=0.0
            planet[a].angle=acos((planet[a].x*cos(planet[a].azimuth)-planet[a].y*sin(planet[a].azimuth))/(planet[a].distance*cos(planet[a].inclination)))
            if planet[a].x*sin(planet[a].azimuth)+planet[a].y*cos(planet[a].azimuth)<0:
                planet[a].angle=2*pi-planet[a].angle
            planet[a].angle=planet[a].angle-planet[a].azproj
            if planet[a].angle<0.0:
                planet[a].angle=planet[a].angle+2*pi
            planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
            for b in range (1,planet_number+1):
                if b!=a:
                    # toavoid calculating the acc from itself
                    planet[a].da=planet[a].da+effaccplanet(planet[a].angle, planet[a].ecc,planet[a].distance,planet[a].inclination, planet[a].azproj,planet[a].azimuth,planet[a].latus, planet[a].x,planet[a].y,planet[a].z,planet[b].x, planet[b].y,planet[b].z)*G*planet[b].mass              
            if repeats==(repeat-1):
                planet[a].accTotal=planet[a].accTotal+abs(planet[a].da)
            planet[a].dv=planet[a].dv+planet[a].da*planet[1].period/precision
            planet[a].velocity=(G*starMass*(2/planet[a].distance-1/planet[a].smAxis))**(0.5)+planet[a].dv-planet[a].vchange       
            planet[a].angle=planet[a].angle+planet[a].velocity*planet[1].period*abs(cos(angle5(planet[a].ecc,planet[a].angle+planet[a].azproj)))/(precision*planet[a].distance)
            planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
        for a in range(1,planet_number+1):
            planet[a].x=planet[a].distance*(cos(planet[a].angle+planet[a].azproj)*cos(planet[a].inclination)*cos(planet[a].azimuth)+sin(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth))
            planet[a].y=planet[a].distance*(sin(planet[a].angle+planet[a].azproj)*cos(planet[a].azimuth)-cos(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth)*cos(planet[a].inclination))
            planet[a].z=planet[a].distance*sin(planet[a].inclination)*cos(planet[a].angle+planet[a].azproj)
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
    planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
    planet[a].x=planet[a].distance*(cos(planet[a].angle+planet[a].azproj)*cos(planet[a].inclination)*cos(planet[a].azimuth)+sin(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth))
    planet[a].y=planet[a].distance*(sin(planet[a].angle+planet[a].azproj)*cos(planet[a].azimuth)-cos(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth)*cos(planet[a].inclination))
    planet[a].z=planet[a].distance*sin(planet[a].inclination)*cos(planet[a].angle+planet[a].azproj)
    planet[a].dv=0.0
    planet[a].accMean=planet[a].accTotal/t
dt=0.0
t=0.0
stepCounter=0
while t<planet[1].period*10:
    stepCounter=stepCounter+1
    for a in range(1,planet_number+1):
        planet[a].da=0.0
        planet[a].angle=acos((planet[a].x*cos(planet[a].azimuth)-planet[a].y*sin(planet[a].azimuth))/(planet[a].distance*cos(planet[a].inclination)))
        if planet[a].x*sin(planet[a].azimuth)+planet[a].y*cos(planet[a].azimuth)<0:
            planet[a].angle=2*pi-planet[a].angle
        planet[a].angle=planet[a].angle-planet[a].azproj
        if planet[a].angle<0.0:
            planet[a].angle=planet[a].angle+2*pi
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
        for b in range (1,planet_number+1):
            if b!=a:
                # toavoid calculating the acc from itself
                planet[a].da=planet[a].da+effaccplanet(planet[a].angle, planet[a].ecc,planet[a].distance,planet[a].inclination, planet[a].azproj,planet[a].azimuth,planet[a].latus, planet[a].x,planet[a].y,planet[a].z,planet[b].x, planet[b].y,planet[b].z)*G*planet[b].mass              
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
                
#may have a problem for dt if one transiting and other not quite
        if planet[a].x**2+planet[a].z**2<(2*(starRadius+planet[a].radius))**2 and planet[a].y>0:            
            #for 3 seconds precision stepnumber*1000 is 10000000            
            dt=planet[1].period/(10000000)
        elif planet[a].x**2+planet[a].z**2<(10*(starRadius+planet[a].radius))**2 and planet[a].y>0:            
            #for 3 seconds precision stepnumber*1000 is 100000000            
            dt=planet[1].period/(1000000)
    t=t+dt
    
    for a in range(1,planet_number+1):
        planet[a].dv=planet[a].dv+planet[a].da*dt
        planet[a].velocity=(G*starMass*(2/planet[a].distance-1/planet[a].smAxis))**(0.5)+planet[a].dv-planet[a].vchange   
        planet[a].angle=planet[a].angle+planet[a].velocity*dt*abs(cos(angle5(planet[a].ecc,planet[a].angle+planet[a].azproj)))/(planet[a].distance)
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
    sunX=0.0
    sunZ=0.0
    for a in range(1,planet_number+1):
        planet[a].x=planet[a].distance*(cos(planet[a].angle+planet[a].azproj)*cos(planet[a].inclination)*cos(planet[a].azimuth)+sin(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth))
        planet[a].y=planet[a].distance*(sin(planet[a].angle+planet[a].azproj)*cos(planet[a].azimuth)-cos(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth)*cos(planet[a].inclination))
        planet[a].z=planet[a].distance*sin(planet[a].inclination)*cos(planet[a].angle+planet[a].azproj)
        sunX=sunX-planet[a].x*planet[a].mass/starMass
        sunZ=sunZ-planet[a].z*planet[a].mass/starMass

    #finding the transits of the planets
    for a in range(1,planet_number+1):
        if ((planet[a].x-sunX)**2+(planet[a].z-sunZ)**2)<(starRadius+planet[a].radius)**2 and planet[a].y>0 and planet[a].transit==0:
            planet[a].transitCounter=planet[a].transitCounter+1            
            planet[a].transit=1    
#the duration will be different as it doesn't necessarily transit thewhoel diameter
            duration =transitLength(planet[a].latus,planet[a].ecc,planet[a].inclination,planet[a].azimuth,planet[a].azproj,planet[a].radius,sunX,sunZ,starRadius)/planet[a].velocity
            print 'Transit start for planet ' + str(a) + ' at ' + str(t) + ' with duration ' +str(duration)
            planet[a].transitTimes.append(t)            
            planet[a].lastTransit=t
        elif planet[a].transit==1 and not (((planet[a].x-sunX)**2+(planet[a].z-sunZ)**2)<(starRadius+planet[a].radius)**2 and planet[a].y>0):
            planet[a].transit=0
    plot2.append(t)
for a in range(1,planet_number +1):
    if planet[a].transitCounter>1:
        planet[a].averagePeriod=float((planet[a].lastTransit-planet[a].transitTimes[0]))/(planet[a].transitCounter-1)
        print 'change in period ' + str(planet[a].averagePeriod-planet[a].period)
    else:
        print 'period is not calculated for planet ' + str(a)
        averagePeriod=planet[1].period
    for b in range (0, len(planet[a].transitTimes)):
        planet[a].ominusc.append(planet[a].transitTimes[b]-(planet[a].averagePeriod*(b+1)))
        

plot(planet[a].ominusc)
    
print stepCounter