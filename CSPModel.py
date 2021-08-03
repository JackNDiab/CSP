# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 08:37:44 2021

@author: jdiab
"""

#Stability of Colloids based on charge-based properties (gravity independant)
import numpy as np
import matplotlib.pyplot as plt

global phpzc_avg, ci, soln_ph, x1,x2, p_h2o, planet, t

p_h20=1.0273
kb=1.3806 * 10**-23
SE=[]

t=float(input("Temperature (Values outside the range of 0 to 95 degrees C are extrapolated): "))
c=float(input("Effective Salt Concentration (counterion) in mg ions/g of Colloid *Does not represent total salt concentration* *Acuurate upper limit ~ 5mg/g*: "))
pc=int(input("Composition (0=silicate based, 1=magnetite based, 2=Hematite based, 3= Rutile based, 4=carbonate based): "))
ci=int(input("Confidence interval of collod stability (75, 80, 90, 95) *80 reccomended*: "))
soln_ph=float(input("pH of solution: "))
planet= int(input("Planetary body; Ceres (0), Enceladus(1), Earth(2): "))

#comp returns a integer that accounts for differences in the zero point charge pH (phzpc) between particles of differing composition
#comp correction is based on literature values of 'Schoonen 1994' as well as other emperical phzpc data.
def comp(p):
    if pc==0:
        return -3
    elif pc==1:
        return 2
    elif pc==2:
        return 4
    elif pc==3:
        return 0
    elif pc==4:
        return 4.2

#Magnetite pzc w/ respect to temp based on 'Schoonen 1994'.
#Experimental values do not exceed 100 degrees C, up to 350 degrees C has been extrapolated
phpzc_magT=7.4318 - 0.02299*t + 7.32e-5 * (t**2) - 1.24e-7 * (t**3) + (comp(pc))
phpzc_magT2=7.2806 - 0.03248*t + 0.000152*(t**2) - 2.72e-7 * (t**3) + (comp(pc))

#Silica pzc w/ respect to "salt concentration" (concentration of counter ions i.e. cations - anions or vice versa) based on 'Milonjic 2007'.
#these ions destabilize the colloid pushing the phzpc up. This equation is based on very limited data and is not apllicable at concs. higher than around 5%
phpzc_siC=1.854 * c + 3.2

#Magnetite formula shape applies generaly to many colloidal particles and therfore is used as the basis model.
#Emperically-based scaling of phzpc with respect to ion concentration effects is applied using an average.
#Its a rudimentary implementation that increases the accuracy of the phzpc calculations within reasonable temperature range (-20 to 100).
phpzc_avg=(phpzc_magT+phpzc_magT2+phpzc_siC)/3
print(f"\n\033[4m*** Output ***\033[0m \n\npH of Point of Zero Charge is {phpzc_avg}")


class Stability: 
    #I created the class so that I can call and print individual funcs if needed to check output.
    #It can also provide more functionality if program is expanded
    #
    #---Output---
    #Program will generate potential size ranges in nm based on Temp, pH, Colloid Composition, and Salt Concentration
    
    def stat():
    #stat returns x1 and x2 values based on Confidence Interval (C.I.)
    #x1 and x2 represent the pH values bounding the "region of instability", higher C.I. = increased probability that the estimated size range output will bestable
        X=(2.01, 2.24, 2.88, 3.43)
        
        def substat(y):
            #sub func for making simple raw x1 and x2 calculations 
            x1= phpzc_avg-y
            x2= phpzc_avg+y
            
            return x1,x2
        
        if ci == 75:
            return substat(X[0]) 
        
        elif ci == 80:
            return substat(X[1])
        
        elif ci == 90:
            return substat(X[2])
        
        elif ci == 95:
            return substat(X[3])
            
            
    def is_stable():
        #Based on pH bounds x1 and x2, is_stable calculates the "raw" size of the colloid
        #xe represents the distance from the "instability boundary" which is used to extrapolate zeta potential based on a emperically based curve (Ostolska and Wiśniewska 2014; Berg et al. 2009)
        x1,x2=Stability.stat()
        
        if 0<=soln_ph<=x1:
            xe=x1-soln_ph
            ZP_magnitude=-40*(pow(0.8, xe))+70
            size_mag_const=-220*np.arctan((ZP_magnitude/10)-4) + 310
            return size_mag_const
        
        elif x2<=soln_ph<=14:

            xe=soln_ph-x2
            ZP_magnitude=-40*(pow(0.8, xe))+70
            size_mag_const=-220*np.arctan((ZP_magnitude/10)-4) + 310
            return size_mag_const
    
        else:
            i=input(print("Colloid unstable. Start instability analysis? (y/n): "))
            return i
        
        
    def estimate_size():
        #Uses the size magnitue constant from is_stable to calculate corrections based on user defined composition
        #Range is represented as an arbitrary uncertainty based on characteristics of each composition
        if Stability.is_stable() == float():
            smc=Stability.is_stable()
        else:
            return
        
        if pc==0:
            sr0=smc/1.3
            usr0=sr0/5
            srl0=[sr0,usr0]
            print(f"\nSilica size range estimate is {sr0} +/- {usr0}")
            return SE.append(srl0)
            
        elif pc==1:
            sr1=smc/4.5
            #Based on Thapa et al 2004
            usr1=sr1/5
            srl1=[sr1, usr1]
            print(f"\nMagnetite size range estimate is {sr1} +/- {usr1}")
            return SE.append(srl1)
            
        elif pc==2:
            sr2=smc/4
            #Based on Supattarasakda et al 2013
            usr2=sr2/5
            srl2=[sr2,usr2]
            print(f"\nHematite particle range estimate is {sr2} +/- {usr2}")
            if (smc/4)<30:
                print("\nWhile parameters are stable, particle size is too small to represent a stable colloid for this composition (<30 nm)")
            return SE.append(srl2)

        elif pc==3:
            sr3=smc*1.67
            usr3=sr3/4
            srl3=[sr3, usr3]
            print(f"\nRulite size range estimate is {sr3} +/- {usr3} ")
            return SE.append(srl3)
        
        elif pc==4:
            sr4=smc*1.67
            usr4=sr4/3
            srl4=[sr4, usr4]
            print(f"\nCarbonate size range estimate is {sr4} +/- {usr4} ***Note*** Carbonate crystals form very rapidly and therefore ")
            return SE.append(srl4)
        
    def flocculation():
        #Calculates flocculation parameters and returns statements regarding the nature of the flocculation behavior
        def p_particle():
            #Returns density based on user composition input
            if pc==0:
                return 2.6
            elif pc==1:
                return 5.3
            elif pc==2:
                return 5.2
            elif pc==3:
                return 4.23
            elif pc==4:
                return 2.8
            
        def grav():
            #returns gravitational constants in terms of user inputed planetary body
            if planet == 0:
                return 0.36
            if planet == 1:
                return 0.113
            if planet == 2:
                return 9.8
            
        dp=p_particle() - p_h20
        g= grav()
        radii = [a for a in range(int((float(SE[0][0]) - float(SE[0][1]))/2), int((float(SE[0][0]) + float(SE[0][1]))/2), 1)]
        #radii is a list of radii determined based on sizes from estimate_size
        
        def peclet(r):
            #Calculates peclet number
            return (2 * np.pi * g * dp * pow((r/(1e-9)), 4))/(3 * kb * (t + 0.000001))
        
        list1=list(map(peclet, radii))
        print(list1)
        bools = []
        #Determines if peclet numbers calculated based on a range of radii are >> or <<1 (indicating ortho or perikinetic)
        #creates a 'bool-like' list of 0,1,or 2 based on peclet number characteristics and then calculates percentage of each characteristic
        for j in list1:
            if j >1000:
                bools.append(1)
            elif j < 0.001:
                bools.append(2)
            else:
                bools.append(0)
                
        periK = bools.count(1)
        orthoK = bools.count(2)
        diff = bools.count(0)
        percent_pk = periK/(len(bools))
        percent_diff = diff/(len(bools))
        percent_ok = orthoK/(len(bools))
        if percent_pk > 0.5:
            print(f"Flocculation is approximately {percent_pk * 100}% perikenetic")
            
        elif percent_ok > 0.5:
            print(f"Flocculation is approximately {percent_ok * 100}% orthokenetic")
        else:
            print(f"Flocculation is approximately {percent_diff * 100}% diffusive. Niether sedimentation or mass transport domiates ({percent_pk * 100}% PK, {percent_ok * 100}% OK)")
            
#size simply controls the path of functions executed based on necessity and user choice
    def size():
        Stability.is_stable()
        
        #This just allows the program to skip estimate_size calculations if there is no need
        #this will be more helpful as more taxing calculations are added and functionality for characterizing flocculation is added
        if Stability.is_stable() == 'y':
            Stability.estimate_size()
            Stability.flocculation()
        #last 2 statements do the same thing. In the future this might be differentiated.
        elif Stability.is_stable() == 'n':
            Stability.estimate_size()
        else:
            Stability.estimate_size()
            print("\nIF destabilization under these conditions/parameters occurs:")
            Stability.flocculation()

Stability.size()
            