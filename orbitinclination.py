# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 08:55:40 2014

@author: chrisbennett
"""

#pylab allows us to use plotting tools
from pylab import*
import matplotlib.pylab as plt
from math import pi, cos, sin, acos, asin, atan
G=6.673848e-11
#class containing all of thevariables associated with each individual planet
class Planet(object):
    def __init__(self,periodYears,radius,mass,ecc,inclination,azimuth):
        self.periodYears=periodYears
        self.period=float(self.periodYears)*365.25*3600*24
        self.radius=radius
        self.mass=mass
        self.ecc=ecc
        #inclination is the angle between the xy plane and theorbital plane
        self.inclination=inclination
        #azimuth is the angle of rotation of the system wrt the x axis, 0 means nearest point of orbit is on the x axis  
        self.azimuth=azimuth
        #the azimuthal angle projected on the orbit is the azproj
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
        #angle is in the orbital plane wrt the x axis
        self.angle=pi/2
        self.distance=self.latus/(1+ecc*cos(self.angle+self.azproj))
        #acceleration due to other planets is da and it affects dv, the change in velocity from single planet orbit
        self.da=0.0
        self.dv=0.0
        #transit is the detector for if the planet was transiting on its last step
        self.transit=1
        #gives the final transit time for use in average time calculation
        self.lastTransit=0
        #counts transits
        self.transitCounter=0
        #for use calculating the velocity shift correction, orbit counter makes sure it stops calculating on an integer number of orbits, orbit test makes it add one for each orbit
        self.vchange=0.0
        self.orbitCounter=0
        self.orbitTest=1
        #accTotal and accMean is used to calculate how big steps we should use
        self.accTotal=0.0
        self.accMean=0.0
        #dt is thestep size that should be taken according to that planet
        self.dt=0.0
        #for calculating vchange
        self.dvlist=[]
        #a list of when the transits occur
        self.transitTimes=[]
        self.averagePeriod=0.0
        #observed minus the calculated times for a linear emidermis, this is what is plotted eventually
        self.ominusc=[]
#effaccplanet calculates the acceleration* due to each planet in the direction of travel. *The acceleration but without G or mass included
def effaccplanet(angle, ecc, distance,i,azproj,azimuth,latus,x1,y1,z1, x2, y2,z2):
    #calculated by differentiation
    dx=distance*((cos(angle+azproj)*cos(i)*cos(azimuth)+sin(angle+azproj)*sin(azimuth))*distance*ecc*sin(angle+azproj)/latus-cos(i)*cos(azimuth)*sin(angle+azproj)+cos(angle+azproj)*sin(azimuth))
    dy=distance*((-cos(angle+azproj)*cos(i)*sin(azimuth)+sin(angle+azproj)*cos(azimuth))*distance*ecc*sin(angle+azproj)/latus+cos(i)*sin(azimuth)*sin(angle+azproj)+cos(angle+azproj)*cos(azimuth))
    dz=-distance*sin(i)*(distance*ecc*cos(angle+azproj)*sin(angle+azproj)/latus-sin(angle+azproj))
    r=((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**(0.5)    
    effaccplanet=((x2-x1)*dx+(y2-y1)*dy+(z2-z1)*dz)/((dx**2+dy**2+dz**2)**(0.5)*r**3)
    return effaccplanet
#this gives the angle between the travel direction and a tangent if it were a circle, used for calculating the change in angle and the transit length, best used with abs(cos(angle5)) as it is only tested for that so far and sign errors may occur
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
#this returns the angle that the transit starts( this angle is in the orbit plane wrt the nearest point) This may miss a transit if the radius to orbit ratio is too small, will write a dtectetor for this sometime
def transitStart(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,precision,angle1):
    #solution remains 0.0 if no transit occurs
    solution=0.0
    #distance is actually distance from planet to sun squared minus the the radii sum squared so has units m^2
    distance=((latus/(1+ecc*cos(angle1)))*(cos(angle1)*cos(inc)*cos(azimuth)+sin(angle1)*sin(azimuth))-xs)**2+((latus/(1+ecc*cos(angle1)))*sin(inc)*cos(angle1)-zs)**2-(rp+rs)**2
    #sets start to be correct initially, start=1 if transit is off    
    if distance==abs(distance) or (latus/(1+ecc*cos(angle1-azproj)))*(-cos(angle1-azproj)*cos(inc)*sin(azimuth)+sin(angle1-azproj)*sin(azimuth))<0:
        start=1
    else:
        start=0
    #calcualted distance for every angle, used 100000 to initially calculate for the whole circle, 10000 is the semsitivity it can detect to, if the ratio of radius to orbit is smaller than this the transit may be missed. This will be noticable as the transit will say it has o length but is still detected 
    for a in range(0,10000):
        angle = a*2*pi/(precision*10000)+angle1   
        #if it is in the orbit take a step back and do it again 10 times as precise until the precision is reached, then it returns the angle throuh the loops
        distance=((latus/(1+ecc*cos(angle)))*(cos(angle)*cos(inc)*cos(azimuth)+sin(angle)*sin(azimuth))-xs)**2+((latus/(1+ecc*cos(angle)))*sin(inc)*cos(angle)-zs)**2-(rp+rs)**2
        if distance==abs(distance) or (latus/(1+ecc*cos(angle1-azproj)))*(-cos(angle1-azproj)*cos(inc)*sin(azimuth)+sin(angle1-azproj)*sin(azimuth))<0:
            start=1
        if distance!=abs(distance) and start==1:
            if precision<10000000000:
                solution=transitStart(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,precision*10,angle-2*pi/(precision*10000))
            else:
                return angle1
            break
    #0 is returned if the transit is not detected
    if solution==0.0:
        return 0.0
    else:
        return solution
#almost thesame code but it detects the end of the tranist angle
def transitEnd(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,precision,angle1):
    solution=0.0    
    distance=((latus/(1+ecc*cos(angle1)))*(cos(angle1)*cos(inc)*cos(azimuth)+sin(angle1)*sin(azimuth))-xs)**2+((latus/(1+ecc*cos(angle1)))*sin(inc)*cos(angle1)-zs)**2-(rp+rs)**2
    if distance==abs(distance) or (latus/(1+ecc*cos(angle1-azproj)))*(-cos(angle1-azproj)*cos(inc)*sin(azimuth)+sin(angle1-azproj)*sin(azimuth))<0:
        start=1
    else:
        start=0
    for a in range(0,10000):
        angle = a*2*pi/(precision*10000)+angle1    
        distance=((latus/(1+ecc*cos(angle)))*(cos(angle)*cos(inc)*cos(azimuth)+sin(angle)*sin(azimuth))-xs)**2+((latus/(1+ecc*cos(angle)))*sin(inc)*cos(angle)-zs)**2-(rp+rs)**2
        if distance!=abs(distance):
            start=0
        if (distance==abs(distance) or (latus/(1+ecc*cos(angle1-azproj)))*(-cos(angle1-azproj)*cos(inc)*sin(azimuth)+sin(angle1-azproj)*sin(azimuth))<0) and start==0:
            if precision<10000000000:
                solution=transitEnd(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,precision*10,angle-2*pi/(precision*10000))
            else:
                return angle1
            break
    if solution==0.0:
        return 0.0
    else:
        return solution
#initiates the code for finding the start and end of transit angles and uses this to give the length of the transit in meters
def transitLength(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs):
    angleEnd=transitEnd(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,1,0)
    angleStart=transitStart(latus,ecc,inc,azimuth,azproj,rp,xs,zs,rs,1,0)
    length=(angleEnd-angleStart)*(latus/(1+ecc*cos(angleEnd-azproj)))/abs(cos(angle5(ecc,angleEnd-azproj)))
    print angleEnd, angleStart    
    return length

#star inputs, could make this into a class but not really necessary
starMassInMs = raw_input('What is the stars mass in solar masses?')
starMass=float(starMassInMs)*1.98e30
starRadius=6.995e8*float(raw_input('star radius in sun radii'))

#system inputs, step number is now an approximation as alot of the time it is predetermined
step_number=int(raw_input("How many steps?"))
#makes the array tht stores each planet, i added a 0 so the planets could be referenced starting at 1
planet=[0]
planet_number=int(raw_input('How many planets?'))

#adding planets
for a in range(0,planet_number):
    print 'Planet ' + str(a+1)
    planet.append(Planet(float(raw_input("What is the period in years?")),6.371e6*float(raw_input('radius in earth radii?')),6e24*float(raw_input('planet mass in earth masses?')),float(raw_input("What is the eccentricity?")),2*pi*float(raw_input('inclination in degrees?(between ±90)'))/360,2*pi*float(raw_input('azimuthal angle?(in degrees)'))/360))
#plot1 and 2 are used for plotting graphs x af y axis
plot1=[]
plot2=[]

#code to find average dv

#teststeps gives the number of orbits the 1st planet that are required to give a complete cycle of the ttv as to avoid the velocity shift, it is an integer though so an error coud occur when using fractional orbits. possibly need to find lowest integer multiple
testSteps=1
for a in range(1,planet_number+1):
    testSteps=testSteps*int(planet[a].periodYears)
print testSteps
#precision gives the number of steps used per planet 1 orbit
precision=1000
#repeat is thenumber of iterations used to get the velocity shift, sometimes 2 is enouh, sometimes we need many more
repeat=2
for repeats in range(0,repeat):
    #reset all the values to original ones, with the plnets all at pi/2
    for a in range(1,planet_number+1):
        planet[a].angle=pi/2
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
        planet[a].x=planet[a].distance*(cos(planet[a].angle+planet[a].azproj)*cos(planet[a].inclination)*cos(planet[a].azimuth)+sin(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth))
        planet[a].y=planet[a].distance*(sin(planet[a].angle+planet[a].azproj)*cos(planet[a].azimuth)-cos(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth)*cos(planet[a].inclination))
        planet[a].z=planet[a].distance*sin(planet[a].inclination)*cos(planet[a].angle+planet[a].azproj)
        #deleted some code here that reversed the z coordinate, i think it was just out of date and is now inline with the othercode        
        planet[a].dv=0.0
        planet[a].orbitCounter=0
        planet[a].orbitTest=1
        del planet[a].dvlist[:]
    for t in range(1,int(testSteps+1)*precision):
        for a in range(1,planet_number+1):
            planet[a].da=0.0
            #calculates the angle from the positions            
            planet[a].angle=acos((planet[a].x*cos(planet[a].azimuth)-planet[a].y*sin(planet[a].azimuth))/(planet[a].distance*cos(planet[a].inclination)))
            if planet[a].x*sin(planet[a].azimuth)+planet[a].y*cos(planet[a].azimuth)<0:
                planet[a].angle=2*pi-planet[a].angle
            planet[a].angle=planet[a].angle-planet[a].azproj
            if planet[a].angle<0.0:
                planet[a].angle=planet[a].angle+2*pi
            #calculates the distance from the angle
            planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
            #add the acceleration from all other planets to give thetotal acceleration on the planet in the direction of its orbit     
            for b in range (1,planet_number+1):
                if b!=a:
                    # toavoid calculating the acc from itself
                    planet[a].da=planet[a].da+effaccplanet(planet[a].angle, planet[a].ecc,planet[a].distance,planet[a].inclination, planet[a].azproj,planet[a].azimuth,planet[a].latus, planet[a].x,planet[a].y,planet[a].z,planet[b].x, planet[b].y,planet[b].z)*G*planet[b].mass              
            #for the last repeat of the process calculate the total accelreation, used to calculate step sizes
            if repeats==(repeat-1):
                planet[a].accTotal=planet[a].accTotal+abs(planet[a].da)
            #calculate the new velocity change (from the original orbit) for a tim step planetperiod1/precision
            planet[a].dv=planet[a].dv+planet[a].da*planet[1].period/precision
            #vchange is taken from this, it is the velocity shift calcualted on theprevious repeat
            planet[a].velocity=(G*starMass*(2/planet[a].distance-1/planet[a].smAxis))**(0.5)+planet[a].dv-planet[a].vchange       
            #recalculateangleand distance for each planet            
            planet[a].angle=planet[a].angle+planet[a].velocity*planet[1].period*abs(cos(angle5(planet[a].ecc,planet[a].angle+planet[a].azproj)))/(precision*planet[a].distance)
            planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
        for a in range(1,planet_number+1):
            #calculate new x,y,z for each planet. Needs to be done in a seperate loop to before as these affect the accelerations
            planet[a].x=planet[a].distance*(cos(planet[a].angle+planet[a].azproj)*cos(planet[a].inclination)*cos(planet[a].azimuth)+sin(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth))
            planet[a].y=planet[a].distance*(sin(planet[a].angle+planet[a].azproj)*cos(planet[a].azimuth)-cos(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth)*cos(planet[a].inclination))
            planet[a].z=planet[a].distance*sin(planet[a].inclination)*cos(planet[a].angle+planet[a].azproj)
            #counts the number of orbits for each planet, orbit test works out which half of thetransit it is in so adds one to the counter when it changes            
            if planet[a].x<0.0 and planet[a].orbitTest==0:
                planet[a].orbitCounter=planet[a].orbitCounter+1
                planet[a].orbitTest=1
            elif planet[a].orbitTest==1 and planet[a].x>0.0:
                planet[a].orbitTest=0
        #when the planets have done the full ttvperiod the system stops an run again
        if planet[1].orbitCounter==int(testSteps):
            break
        #add the dv value to thelist for use in calculating average
        for a in range(1,planet_number+1):
            planet[a].dvlist.append(planet[a].dv)
    for a in range(1,planet_number+1):
        #vchange is the average velocity change, if itbecoms constant it is correct
        planet[a].vchange=sum(planet[a].dvlist)/len(planet[a].dvlist)
        print a,planet[a].vchange

#resetting all of the terms
for a in range(1,planet_number+1):
    planet[a].angle=pi/2
    planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
    planet[a].x=planet[a].distance*(cos(planet[a].angle+planet[a].azproj)*cos(planet[a].inclination)*cos(planet[a].azimuth)+sin(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth))
    planet[a].y=planet[a].distance*(sin(planet[a].angle+planet[a].azproj)*cos(planet[a].azimuth)-cos(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth)*cos(planet[a].inclination))
    planet[a].z=planet[a].distance*sin(planet[a].inclination)*cos(planet[a].angle+planet[a].azproj)
    planet[a].dv=0.0
    # the t involved here is thetotla number of steps used in the section calculating vchange
    planet[a].accMean=planet[a].accTotal/t
#setting dt and t to 0 as they give the time and timestep in the next section
dt=0.0
t=0.0
#stepcounter gives thetotal number of steps used, it is just for interest and can be used to give a time estimate for code completion
stepCounter=0
#the main code simulating the transit, runs until a certain number of planet1 periods are complete. 
while t<planet[1].period*10:
    stepCounter=stepCounter+1
    for a in range(1,planet_number+1):
        planet[a].da=0.0
        #calculate the angle dependent on the positions (i think it give between 0 and 360 in degrees)
        planet[a].angle=acos((planet[a].x*cos(planet[a].azimuth)-planet[a].y*sin(planet[a].azimuth))/(planet[a].distance*cos(planet[a].inclination)))
        if planet[a].x*sin(planet[a].azimuth)+planet[a].y*cos(planet[a].azimuth)<0:
            planet[a].angle=2*pi-planet[a].angle
        planet[a].angle=planet[a].angle-planet[a].azproj
        if planet[a].angle<0.0:
            planet[a].angle=planet[a].angle+2*pi
        #calcualting the distance
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
        #calculates the acceleration due to each planet in the orbit direction        
        for b in range (1,planet_number+1):
            if b!=a:
                # toavoid calculating the acc from itself
                planet[a].da=planet[a].da+effaccplanet(planet[a].angle, planet[a].ecc,planet[a].distance,planet[a].inclination, planet[a].azproj,planet[a].azimuth,planet[a].latus, planet[a].x,planet[a].y,planet[a].z,planet[b].x, planet[b].y,planet[b].z)*G*planet[b].mass              
        #calculating the step size, for no acceleration use the given precision/10
        if planet[a].da==0:
            planet[a].dt=planet[1].period*10/(step_number)
        else:
            #else dt is the normal step number time the acceleration/accmean
            planet[a].dt=planet[a].accMean*planet[1].period/(step_number*abs(planet[a].da))
    #dt set to the first planets one for starters, will change
    dt=planet[1].dt
    for a in range(1,planet_number+1):
        #finds the smallest dt, allows no bigger steps stepnumber/10
        if planet[a].dt<=dt:
            if planet[a].dt<planet[1].period*10/(step_number):
                dt=planet[a].dt
            else:
                dt=planet[1].period*10/(step_number)
                
#may have a problem for dt if one transiting and other not quite
        #if the planet is near a transit dt should be smaller for all, could get rid of this once thetransit has started
        if planet[a].x**2+planet[a].z**2<(2*(starRadius+planet[a].radius))**2 and planet[a].y>0:            
            #for 3 seconds precision stepnumber*1000 is 10000000            
            dt=planet[1].period/(10000000)
        elif planet[a].x**2+planet[a].z**2<(10*(starRadius+planet[a].radius))**2 and planet[a].y>0:            
            #for 3 seconds precision stepnumber*1000 is 100000000            
            dt=planet[1].period/(1000000)
    #step on in time by the calculated dt
    t=t+dt
    for a in range(1,planet_number+1):
        #step on all variables by the dt and recalculate angle and distance
        planet[a].dv=planet[a].dv+planet[a].da*dt
        planet[a].velocity=(G*starMass*(2/planet[a].distance-1/planet[a].smAxis))**(0.5)+planet[a].dv-planet[a].vchange   
        planet[a].angle=planet[a].angle+planet[a].velocity*dt*abs(cos(angle5(planet[a].ecc,planet[a].angle+planet[a].azproj)))/(planet[a].distance)
        planet[a].distance=planet[a].latus/(1+planet[a].ecc*cos(planet[a].angle+planet[a].azproj))
    #for calculating the sun position for transit detection, didn't do y as its not relevant but could easily
    sunX=0.0
    sunZ=0.0
    for a in range(1,planet_number+1):
        #recalculate the positions usin the new angles and use these to find thenew sun position
        planet[a].x=planet[a].distance*(cos(planet[a].angle+planet[a].azproj)*cos(planet[a].inclination)*cos(planet[a].azimuth)+sin(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth))
        planet[a].y=planet[a].distance*(sin(planet[a].angle+planet[a].azproj)*cos(planet[a].azimuth)-cos(planet[a].angle+planet[a].azproj)*sin(planet[a].azimuth)*cos(planet[a].inclination))
        planet[a].z=planet[a].distance*sin(planet[a].inclination)*cos(planet[a].angle+planet[a].azproj)
        sunX=sunX-planet[a].x*planet[a].mass/starMass
        sunZ=sunZ-planet[a].z*planet[a].mass/starMass

    #finding the transits of the planets
    for a in range(1,planet_number+1):
        #detects if it is transiting and is inpositive y and transit detects if we have already noticed it
        if ((planet[a].x-sunX)**2+(planet[a].z-sunZ)**2)<(starRadius+planet[a].radius)**2 and planet[a].y>0 and planet[a].transit==0:
            #counts thenumber of transits, not including th first one as transit is set to 1            
            planet[a].transitCounter=planet[a].transitCounter+1            
            planet[a].transit=1   
            #uses the code i wrote to give the length of transit, using only the initial sun position and velocity
            duration =transitLength(planet[a].latus,planet[a].ecc,planet[a].inclination,planet[a].azimuth,planet[a].azproj,planet[a].radius,sunX,sunZ,starRadius)/planet[a].velocity
            print 'Transit start for planet ' + str(a) + ' at ' + str(t) + ' with duration ' +str(duration)
            #crates a list of the beginning of the transit times            
            planet[a].transitTimes.append(t)
            # this variable gives the time of start of the last transit detected, used togive average periods            
            planet[a].lastTransit=t
        #detects it if nottransit is happening
        elif planet[a].transit==1 and not (((planet[a].x-sunX)**2+(planet[a].z-sunZ)**2)<(starRadius+planet[a].radius)**2 and planet[a].y>0):
            #variable to say no transit occured on the previous step            
            planet[a].transit=0
    #pltting the x axis of time
    plot2.append(t)
    
#data handlers
for a in range(1,planet_number +1):
    if planet[a].transitCounter>1:
        #calculates theaverage period using the transit start times
        planet[a].averagePeriod=float((planet[a].lastTransit-planet[a].transitTimes[0]))/(planet[a].transitCounter-1)
        #this output should give 1 if vchange has worked properly, if not the period is inorrect
        print 'period fraction ' + str(planet[a].averagePeriod/planet[a].period)
    else:
        #lets us know if not enough transits were detected to give a result and takes the kepler period as true
        print 'period is not calculated for planet ' + str(a)
        averagePeriod=planet[1].period
    for b in range (0, len(planet[a].transitTimes)):
        #calculates the difference in the observed and calculated transit times using the average period (should be able to use the kepler period as they should be the same)
        planet[a].ominusc.append(planet[a].transitTimes[b]-(planet[a].averagePeriod*(b+1)))
        
#plots a graph
plot(planet[a].ominusc)
#prints thenumber of steps taken (approx stepnumber + transit steps)
print stepCounter