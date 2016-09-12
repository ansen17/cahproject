#!/usr/bin/env python

import subprocess
import multiprocessing
import select
import fcntl, os
import re

import glob
import sys

import argparse
import math
from mike_math import step_function, sign
from dprint import dprint

import numpy as np

from matplotlib import cm
from matplotlib import colorbar
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import time

import evolver
import logging


evolver.logger.setLevel(logging.DEBUG)
evolver.logger.addHandler(logging.StreamHandler())


def main():
    """The main routine."""

    with evolver.Evolver() as E:

        #Where user can specify physical parameters of surface through command line with default values right below

        parser = argparse.ArgumentParser(description='Create a Surface Evolver script.')
        parser.add_argument("-A", action="store", dest="amplitude", type=float,
                    metavar="<amp>", help="Amplitude")
        parser.add_argument("-l", action="store", dest="lambda1", type=float,
                    metavar="<lambda>", help="lambda1")
        parser.add_argument("-n", action="store", dest="grooves", type=int,
                    metavar="<grooves>", help="number of grooves")
        parser.add_argument("-M", action="store", dest="mode",
                    choices=["C","W"],
                    metavar="<wetting mode>", help="wetting mode")
        parser.add_argument("-b", action="store", dest="Bo", type=float,
                    metavar="<bond number>", help="bond number")
        parser.add_argument("-e", action="store", dest="thetaE", type=float,
                    metavar="<thetaE>", help="thetaE in degrees")
        args = parser.parse_args()

        if args.thetaE == None: thetaE = 115.0
        else: thetaE = args.thetaE

        if args.amplitude == None: amplitude = 0.5
        else: amplitude = args.amplitude

        if args.lambda1 == None: lambda1 = 1.0
        else: lambda1 = args.lambda1

        if args.grooves == None: grooves = 3
        else: grooves = args.grooves

        if args.mode == None: mode = "Wenzel"
        else:
            if args.mode == "C":
                mode = "Cassie"
            else:
                mode = "Wenzel"

        if args.Bo == None: Bo = 0.0
        else: Bo = args.Bo

        #-----------------------------Variable Instantiation ------------------------------------------

        flagger = 0                     # Flag variable to indicate if droplet interface crosses sinusoidal surface
        breaker = False

        tangle = 0.00                   # Sinusoidal Surface's Tilt Angle
        thetaE = 115.00                 # Internal contact angle of droplet with the surface
        

        increment = 1.00                # combination of n and increment to decide how slowly to tilt angle of the sinusoidal surface         
        n=0
        
    
        initialleftx = -200             # Right & Left vertex coordinates and apparent contact angles
        initialrightx= -200
        carapp = -300
        calapp = -400
        rightsliding = False            
        leftsliding = False


        ampperlam = amplitude/lambda1   # Amplitude/Wavelength values for test cases
        
        params = np.zeros([7, int(100.0/increment)])
        grooves= int(2.1/lambda1)

        CAH = -1.00
        _nrg_condition_for_E = 7.00
        _nrg_condition_for_e = 8.00
        
        delta_e = 1500                  #Energy constraint values
        delta_E = 1500                  
        E_0 = 1500                      
        e_0 = 100
        e_1 = 200
        E_1 = 900

        minimum_e_for_tilt_angle = -7000    #Placeholder for initial energy value

        '''Main evolution algorithm'''

        while (tangle < 90.00 and breaker==False):

            '''Incrementing Tilt Angle '''
            tangle = n*increment

            
            '''PRERUN 1: Opening multiple dropSinusoidalTest.fe files through .sinusoidal.py file'''

            if (tangle == 0.00):
                os.system("python /storage/home/azs5619/work/scripts/evolverset/sinusoidal.py -A {d} -l {e} -M W -b {bond} -n {f}".format(d=amplitude, e= lambda1, f= grooves, bond = Bo)) 
                while True:
                    try: 
                        E.open_file('dropSinusoidalTest.fe')
                        break
                    except:
                        print 'trying to open the .fe file'

                ''' PRERUN 1: Setting Tilt Angle & Obtaining Vertex Coordinates of Initial Drop'''

                x = E.run_command("tilt_angle:={a}".format(a=tangle))
                E.run_command('print "starting point of fun"')
                t= 0 
                while (t <= 5):
                    E.run_command('{g 10; V 100;} 15;')
                    E.run_command(' ')
                    E.run_command('rz')
                    t = t+1

                E.run_command('')
                E.run_command('print "stalling"')
                while True:
                    try:
                        initialleftx = float(E.run_command('left_vx'))
                        E.run_command('')
                        initialrightx = float(E.run_command('right_vx'))
                        break
                    except:
                        print 'Trying to get Inital Left Vx and Initial Right Vx'

                E.run_command('')
                E.close_file()

            '''Real Run: Same as above with slight modifications and obtaining the initial energy value readings'''

            os.system("python /storage/home/azs5619/work/scripts/evolverset/sinusoidal.py -A {d} -l {e} -M W -b {bond} -n {f}".format(d=amplitude, e= lambda1, f= grooves, bond = Bo)) 
            while True:
                try: 
                    E.open_file('dropSinusoidalTest.fe')
                    break
                except:
                    print 'trying to open the .fe file'

            x = E.run_command("tilt_angle:={a}".format(a=tangle))                                                   
            E.run_command('print "Pretreatment of droplet evolution"')
            E.run_command('{g 1500; V 100;}')
            E.run_command('')
            while True:
                try:
                    averageLength = float(E.run_command('averageLength'))
                    E.run_command('')
                    break
                except:
                    print 'Trying to get Average Length'

            E.run_command('')

            while (averageLength > 0.10 * lambda1):
                E.run_command('{g 10; V 100;}')
                E.run_command('rz')
                E.run_command('')
                while True:
                    try:
                        averageLength = float(E.run_command('averageLength'))
                        E.run_command('')
                        break
                    except:
                        print 'Trying to get Average Length in the while (rz) check loop'
            E.run_command('')
            while True:
                try:
                    thetaright = float(E.run_command('car'))
                    E.run_command('')
                    thetaleft = float(E.run_command('cal'))
                    E.run_command('')
                    e_0 = float(E.evolve())
                    E.run_command('')
                    break
                except:
                    print 'Trying to get e_0, thetaright and thetaleft before entering nested while loops'

            E.run_command('')

            delta_E =100

            ''' Main core of algorithm: check to meet energy and sliding constraints '''
            while ((abs(delta_E) > 10**(-_nrg_condition_for_E) or abs(thetaright-thetaE)>1.00 or abs(thetaleft-thetaE)>1.00) and breaker == False):
                while True:
                    try:
                        e_0 = float(E.evolve())
                        E.run_command('')
                        break
                    except:
                        print 'Trying to get e_0 in the outer loop'
                no_of_g = 0
                delta_e = 100
                while(abs(delta_e) > 10**(-_nrg_condition_for_e)):
                    no_of_g = no_of_g + 1

                    while True:
                        try:
                            e_1 = float(E.evolve())
                            E.run_command('')
                            break
                        except:
                            print 'Trying to get e_0 in the inner loop'
                    E.run_command('')
                        
                    E.run_command('even_out')

                    delta_e = float((e_1 - e_0)/(e_0))
                    E.run_command('')
                    while True:
                        try:
                            leftx = float(E.run_command('left_vx'))
                            print 'Got current Left vx'
                            rightx = float(E.run_command('right_vx'))
                            E.run_command('')
                            break
                        except:
                            print 'Trying to get the current left and right vertex values'
                    E.run_command('')
                    rightdifference = rightx - initialrightx
                    leftdifference = leftx - initialleftx
                    if  rightdifference > lambda1 and leftdifference > lambda1:
                        rightsliding = True
                        leftsliding = True
                        breaker = True
                        break
                    if (no_of_g % 10):
                        check = E.run_command('am_I_close')
                        if (check == "True"):
                            flagger = -1
                            breaker = True
                            break   
                    delta_e = float((e_1 - e_0)/(e_0))           
                    if (delta_e > 0):
                        flagger=1
                        break

                    e_0 = e_1
                E_1 = e_1
                while True:
                    try:
                        thetaright = float(E.run_command('car'))
                        E.run_command('')
                        thetaleft = float(E.run_command('cal'))
                        E.run_command('')
                        break
                    except:
                        print 'Trying to get thetaright and thetaleft before leaving outer nested while loop'

                E.run_command('')
                delta_E = (E_1 - E_0)/(E_0)
                if (delta_E > 0):
                    flagger=1
                E_0 = E_1
                E.run_command('rz')


            ''' Recording all the relevant measurements including the CAH '''

            while True:
                try:
                    E.run_command('')
                    calapp = float(E.run_command('cal_app'))
                    print 'Got calapp'
                    E.run_command('')
                    carapp = float(E.run_command('car_app'))
                    print 'Got carapp'
                    E.run_command('')
                    break
                except:
                    print 'Trying to get APPARENT contact angles'
            E.run_command('')

            CAH = abs(carapp - calapp)
            
            while True:
                try:
                    E.run_command('')
                    params[0, n] = tangle
                    E.run_command('')
                    params[1, n] = carapp
                    E.run_command('')
                    params[2, n] = thetaright
                    E.run_command('')
                    params[3, n] = calapp
                    E.run_command('')
                    params[4, n] = thetaleft
                    E.run_command('')
                    params[5, n] = minimum_e_for_tilt_angle
                    E.run_command('')
                    params[6, n] = CAH
                    print 'Got the param values in '
                    break
                except:
                    print 'Trying so so hard to get these last param values in but they are being difficult'

           
            E.close_file()

            if (rightsliding == False or leftsliding == False):
                minimum_e_for_tilt_angle = e_1

            n = n + 1

        carappholder = params[1,:n]
        calappholder = params[3,:n]

        index_max_carapp = np.argmax(carappholder)
        index_min_calapp = np.argmin(calappholder)

        carappvalue = params[1, index_max_carapp]
        calappvalue = params[3, index_min_calapp]

        CAH = carappvalue - calappvalue
    

        if (tangle >= 90):
            tangle = 90 

        ''' Below, I log the data and graph the data '''

        '''Below, there are various scenarios that have to be independently run through the LionX Clusters. The commented blocks
           show all the possible scenarios that were run'''          

        ''' Checking the contact angles with the thetaE values within 0.5 degrees; if not BREAK from entire program

        if (abs (float(E.run_command('car')) - thetaE) > 7.50):
            print 'car has broken thetaE'
            print 'CAR is ' + str(E.run_command('car'))
            print 'thetaE is ' + str(E.run_command('print thetaE'))
            break
        if(abs (float(E.run_command('cal')) - thetaE) > 7.50):
            print 'cal has broken thetaE'
            print 'CAL is ' + str(E.run_command('cal'))
            print 'thetaE is ' + str(E.run_command('print thetaE'))
            break '''

                
        

        '''Setting amp/lam as 0.1 '''
        
        '''
        
        txtFile = open('LogFileGraphBondNumber_Wavelength_Tangle.txt', 'a') 
        txtFile.write(str(Bo) + '             ')
        txtFile.write(str(lambda1) + '              ')
        txtFile.write(str(tangle) + '            ')
        txtFile.write(str(flagger) + '     \n ')

        txtFile = open('LogFileGraphBondNumber_Amplitude_Tangle.txt', 'a') 
        txtFile.write(str(Bo) + '               ')
        txtFile.write(str(amplitude) + '           ')
        txtFile.write(str(tangle) + '                ')
        txtFile.write(str(flagger) + '       \n ') '''

        

        '''Phase Diagram between bond no. 0 to 1.00 and ampperlam from 0 to 0.3; change amplitude from 0.01 to 0.30 with lambda being fixed at 1.0 '''
        '''
        txtFile = open('Lastbitofampperlamandtanlge.txt', 'a') 
        txtFile.write(str(Bo) + '         ')
        txtFile.write(str(ampperlam) + '              ')
        txtFile.write(str(tangle) + '            ')
        txtFile.write(str(flagger) + '\n ')         '''



        '''Phase Diagram between wavlength from 0 to 1.00 and ampperlam from 0.0 to 0.40; use ampperlam to find corresponding amplitude for each point '''

        txtFile = open('LogFileGraphWavelength_AmpLam_CAH.txt', 'a') 
        txtFile.write(str(lambda1) + '                         ')
        txtFile.write(str(ampperlam) + '                         ')
        txtFile.write(str(CAH) + '                                  ')
        txtFile.write(str(flagger) + '                      \n ')

        txtFile = open('LogFileGraphWavelength_AmpLam_Energy.txt', 'a') 
        txtFile.write(str(lambda1) + '              ')
        txtFile.write(str(ampperlam) + '                ')
        txtFile.write(str(minimum_e_for_tilt_angle) + '                      ')
        txtFile.write(str(flagger) + '      \n ')
        

        print('Bond Number ' + str(Bo) + ' \n' )
        print('Wavelength ' + str(lambda1) + ' \n ')
        print ('Amplitude ' + str(amplitude) + ' \n ')
        print ('Tilt Angle ' + str(tangle))
        #print ('Order of Magnitude that the threshold is on' + str(_nrg_) + ' \n ')
        print ('Contact Angle Hysteresis ' + str(CAH) + '\n')
        print ('Flagger ' + str(flagger) + ' \n ')
        #txtFile.write('Amplitude '+ str(amplitude)+ '\n ')
        txtFile.close()

        E.run_command('')
        #E.run_command('dump "amp{a}bo{b}"'.format(a=amplitude, b=Bo))

        E.run_command('')
        E.close_file()

        print 'To Plot'
        #print 'amp/lambda array', AMPperLAMparams
        #print 'lambda array', LAMparams
        #print 'cah array', CAHparams
        plt.plot(params[0, :], params[5, :])
        plt.ylabel('Energy')
        plt.xlabel('Tilt Angle')
        plt.title('Energy vs. Tilt Angle')
        #print 'tilt angle', params[0,:]


if __name__ == "__main__":
    main()
