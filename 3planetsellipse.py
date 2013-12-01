# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 20:14:10 2013

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

def Effaccstar(orbitAngle,ecc,latus):
    orbitRadius=latus/(1+ecc*cos(orbitAngle))
    orbitX=orbitRadius*cos(orbitAngle)
    orbitY=orbitRadius*sin(orbitAngle)
    if sin(orbitAngle)==0 and orbitX>0:
        gradientAngle=3*pi/2
    elif sin(orbitAngle)==0 and orbitX<0:
        gradientAngle=pi/2
    else:
        gradientAngle=atan(-(cos(orbitAngle)+ecc)/sin(orbitAngle))
    if orbitAngle>pi:
        gradientAngle=gradientAngle-pi
    if gradientAngle<0:
        gradientAngle=gradientAngle+2*pi
    angleDifference=gradientAngle-orbitAngle+pi
    effectiveAcc=-cos(angleDifference)/(orbitRadius**2)
    return effectiveAcc

def Effaccplanet(orbitAngle1, ecc, latus, orbitX2, orbitY2):
    orbitRadius=latus/(1+ecc*cos(orbitAngle1))
    orbitX1=orbitRadius*cos(orbitAngle1)
    orbitY1=orbitRadius*sin(orbitAngle1)
    if sin(orbitAngle1)==0 and orbitX1>0:
        gradientAngle=3*pi/2
    elif sin(orbitAngle1)==0 and orbitX1<0:
        gradientAngle=pi/2
    else:
        gradientAngle=atan(-(cos(orbitAngle1)+ecc)/sin(orbitAngle1))
    if orbitAngle1>pi:
        gradientAngle=gradientAngle-pi
    if gradientAngle<0:
        gradientAngle=gradientAngle+2*pi
    #now we have the correct gradient angle need to compare to position gradients
    if orbitX1==orbitX2 and orbitY2>orbitY1:
        planetAngle=pi/2
    elif orbitX1==orbitX2 and y1>orbitY2:
        planetAngle=3*pi/2
    planetAngle=atan((orbitY2-orbitY1)/(orbitX2-orbitX1))
    if orbitY2>orbitY1 and planetAngle<0:
        planetAngle=planetAngle+pi
    if orbitY2<orbitY1 and planetAngle>0:
        planetAngle=planetAngle+pi
    elif orbitY2<orbitY1 and planetAngle<0:
        planetAngle=planetAngle+ 2*pi
    relativeAngle=planetAngle-gradientAngle    
    
    effAcc=-cos(relativeAngle)/((orbitX1-orbitX2)**2+(orbitY1-orbitY2)**2)
    return effAcc
def AnglePlanetAway(orbitX,orbitY,ecc,orbitAngle):
    if (orbitX==0 and orbitY>0):
        angleAway=pi-atan(-(cos(orbitAngle)+ecc)/sin(orbitAngle))
    elif (orbitX==0 and orbitY<0):
        angleAway=-atan(-(cos(orbitAngle)+ecc)/sin(orbitAngle))
    elif (orbitAngle==0 or orbitAngle==pi):
        angleAway=0
    else:
        angleAway=atan(orbitY/orbitX)+pi/2-atan(-(cos(orbitAngle)+ecc)/sin(orbitAngle))
    return angleAway
    

# inouts, require masses periods radii and step number for 2nd planet orbit

starMass=float(raw_input('What is the stars mass in solar masses?'))*1.98e30
starRadius=6.995e8*float(raw_input('star radius in sun radii'))
stepNumber=int(raw_input("How many steps?"))

#for planet1

period1=365*24*3600*float(raw_input("What is the period of planet 1 in years?"))
smAxis1=(G*starMass*period1**2/(4*pi**2))**(1.0/3)

planetRadius1=6.371e6*float(raw_input('planet radius 1 in earth radii?'))
planetMass1=6e24*float(raw_input('planet mass 1 in earth masses?'))
ecc1=float(raw_input("What is the eccentricity of planet 1's orbit?"))
latus1=smAxis1*(1-ecc1**2)
velocity1=(G*starMass*(2/latus1-1/smAxis1))**(0.5)


period2=365*24*3600*float(raw_input("What is the period of planet 2 in years?"))
smAxis2=(G*starMass*period2**2/(4*pi**2))**(1.0/3)

planetRadius2=6.371e6*float(raw_input('planet radius 2 in earth radii?'))
planetMass2=6e24*float(raw_input('planet mass 2 in earth masses?'))
ecc2=float(raw_input("What is the eccentricity of planet 2's orbit?"))
latus2=smAxis2*(1-ecc2**2)
velocity2=(G*starMass*(2/latus2-1/smAxis2))**(0.5)

period3=365*24*3600*float(raw_input("What is the period of planet 3 in years?"))
smAxis3=(G*starMass*period3**2/(4*pi**2))**(1.0/3)

planetRadius3=6.371e6*float(raw_input('planet radius 3 in earth radii?'))
planetMass3=6e24*float(raw_input('planet mass 3 in earth masses?'))
ecc3=float(raw_input("What is the eccentricity of planet 3's orbit?"))
latus3=smAxis3*(1-ecc3**2)
velocity3=(G*starMass*(2/latus3-1/smAxis3))**(0.5)



planetX3=0.0
planetY3=latus3
planetX2=0.0
planetY2=latus2
planetX1=0.0
planetY1=latus1
plot2=[]
plot1=[]
plot3=[]
transitOn=1
transitOn2=1
transitOn3=1
lastTransitStart=0
transitCounter1=0
transitCounter2=0
transitCounter3=0

for t in range(1,int(4*stepNumber)):
    if (planetX1>0 and planetY1>=0):    
        orbitAngle1=atan(planetY1/planetX1)
    elif (planetX1<0):
        orbitAngle1=atan(planetY1/planetX1)+pi
    elif (planetX1==0 and planetY1>0):
        orbitAngle1=pi/2
    else:
        orbitAngle1=atan(planetY1/planetX1)+2*pi
    if (planetX2>0 and planetY2>=0):    
        orbitAngle2=atan(planetY2/planetX2)
    elif (planetX2<0):
        orbitAngle2=atan(planetY2/planetX2)+pi
    elif (planetX2==0 and planetY2>0):
        orbitAngle2=pi/2
    else:
        orbitAngle2=atan(planetY2/planetX2)+2*pi
    if (planetX3>0 and planetY3>=0):    
        orbitAngle3=atan(planetY3/planetX3)
    elif (planetX3<0):
        orbitAngle3=atan(planetY3/planetX3)+pi
    elif (planetX3==0 and planetY3>0):
        orbitAngle3=pi/2
    else:
        orbitAngle3=atan(planetY3/planetX3)+2*pi
    
    distance1=latus1/(1+ecc1*cos(orbitAngle1))
    distance2=latus2/(1+ecc2*cos(orbitAngle2))
    distance3=latus3/(1+ecc3*cos(orbitAngle3))

    planet_a1=Effaccstar(orbitAngle1,ecc1,latus1)*(G*starMass)+Effaccplanet(orbitAngle1, ecc1, latus1, planetX2, planetY2)*G*planetMass2+Effaccplanet(orbitAngle1, ecc1, latus1, planetX3, planetY3)*G*planetMass3
    planet_a2=Effaccstar(orbitAngle2,ecc2,latus2)*(G*starMass)+Effaccplanet(orbitAngle2, ecc2, latus2, planetX1, planetY1)*G*planetMass1+Effaccplanet(orbitAngle2, ecc2, latus2, planetX3, planetY3)*G*planetMass3
    planet_a3=Effaccstar(orbitAngle3,ecc3,latus3)*(G*starMass)+Effaccplanet(orbitAngle3, ecc3, latus3, planetX1, planetY1)*G*planetMass1+Effaccplanet(orbitAngle3, ecc3, latus3, planetX2, planetY2)*G*planetMass2
    
    velocity1=velocity1+planet_a1*period1/stepNumber
    velocity2=velocity2+planet_a2*period1/stepNumber
    velocity3=velocity3+planet_a3*period1/stepNumber

    orbitAngle1=orbitAngle1+velocity1*period1*abs(cos(AnglePlanetAway(planetX1,planetY1,ecc1,orbitAngle1)))/(stepNumber*distance1)
    orbitAngle2=orbitAngle2+velocity2*period1*abs(cos(AnglePlanetAway(planetX2,planetY2,ecc2,orbitAngle2)))/(stepNumber*distance2)
    orbitAngle3=orbitAngle3+velocity3*period1*abs(cos(AnglePlanetAway(planetX3,planetY3,ecc3,orbitAngle3)))/(stepNumber*distance3)
    
    distance1=latus1/(1+ecc1*cos(orbitAngle1))
    distance2=latus2/(1+ecc2*cos(orbitAngle2))    
    distance3=latus3/(1+ecc3*cos(orbitAngle3))    

    planetX1=distance1*cos(orbitAngle1)
    planetY1=distance1*sin(orbitAngle1)
    planetX2=distance2*cos(orbitAngle2)
    planetY2=distance2*sin(orbitAngle2)
    planetX3=distance3*cos(orbitAngle3)
    planetY3=distance3*sin(orbitAngle3)    
    plot3.append(velocity1)
    #transit detector required
    if planetX1<starRadius+planetRadius1 and planetX1>-starRadius-planetRadius1 and planetY1>0 and transitOn==1:
        transitOn=0
        transitCounter1=transitCounter1+1
        print 'Transit for planet 1 '+ str(t)
        if transitCounter1!=1 and transitCounter1!=2:
            plot1.append(t)
        elif transitCounter1==2:
            firstTransit=t
        lastTransitStart=t
    elif transitOn==0 and not (planetX1<starRadius+planetRadius1 and planetX1>-starRadius-planetRadius1 and planetY1>0):
        transitOn=1
        print 'Transit finished ' + str(t)
 
    if planetY2<starRadius+planetRadius2 and planetY2>-starRadius-planetRadius2 and planetY2>0 and transitOn2==1:
        print 'To2 '+ str(t)
        transitOn2=0
        transitCounter2=transitCounter2+1
    elif transitOn2==0 and not (planetX2<starRadius+planetRadius2 and planetX2>-starRadius-planetRadius2 and planetY2>0):
        transitOn2=1
        print 'Tf' +str(t-1)
    if planetY3<starRadius+planetRadius3 and planetY3>-starRadius-planetRadius3 and planetY3>0 and transitOn3==1:
        print 'To3 '+ str(t)
        transitOn3=0
        transitCounter3=transitCounter3+1
    elif transitOn3==0 and not (planetX3<starRadius+planetRadius3 and planetX3>-starRadius-planetRadius3 and planetY3>0):
        transitOn3=1
        print 'Tf' +str(t-1)


    
plot(plot3)

    