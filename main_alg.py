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


        flagger = 0  # If droplet has crossed surface (use the number 1)//se 0 if the droplet has not crossed surface

        tangle = 0.00
        thetaE = 115.0
        increment = 1.00

        breaker = False
    
        initialleftx = -200
        initialrightx= -200
        ampperlam = amplitude/lambda1
        n=0
        rightsliding = False
        leftsliding = False
        params = np.zeros([7, int(100.0/increment)])
        grooves= int(2.1/lambda1)
        CAH = -1.00
        _nrg_condition_for_E = 7.00
        _nrg_condition_for_e = 8.00
        completebreaker = False
        delta_e = 100
        delta_E = 100
        E_0 = 100
        e_0 = 100
        e_1 = 200
        E_1 = 900
        carapp = -300
        calapp = -400
        minimum_e_for_tilt_angle = -7000

        while (tangle < 90.00 and breaker==False):

            tangle = n*increment

            if (tangle == 0.00):
                os.system("python /storage/home/azs5619/work/scripts/evolverset/sinusoidal.py -A {d} -l {e} -M W -b {bond} -n {f}".format(d=amplitude, e= lambda1, f= grooves, bond = Bo)) 
                while True:
                    try: 
                        E.open_file('dropSinusoidalTest.fe')
                        print 'opened the .fe file'
                        break
                    except:
                        print 'trying to open the .fe file'
                print 'Opening .fe file in Surface Evolver'

                x = E.run_command("tilt_angle:={a}".format(a=tangle))
                E.run_command('print "starting point of fun"')
                t= 0 
                while (t <= 5):
                    E.run_command('{g 10; V 100;} 15;')
                    E.run_command(' ')
                    E.run_command('rz')
                    t = t+1
                #E.run_command('dump "firstdump{b}"'.format(a=ampperlam, b = Bo))
                print 'stop'
                E.run_command('')
                E.run_command('print "stalling"')
                while True:
                    try:
                        initialleftx = float(E.run_command('left_vx'))
                        print 'Got Initialleftx'
                        E.run_command('')
                        initialrightx = float(E.run_command('right_vx'))
                        print 'Got InitialRightx'
                        break
                    except:
                        print 'Trying to get Inital Left Vx and Initial Right Vx'

                E.run_command('')
                print 'Initial Left vx ', initialleftx
                print 'Initial right vx ', initialrightx
                print 'Leaving...Got the initial values'
                E.run_command('')
                E.close_file()

            os.system("python /storage/home/azs5619/work/scripts/evolverset/sinusoidal.py -A {d} -l {e} -M W -b {bond} -n {f}".format(d=amplitude, e= lambda1, f= grooves, bond = Bo)) 
            while True:
                try: 
                    E.open_file('dropSinusoidalTest.fe')
                    print 'opened the .fe file'
                    break
                except:
                    print 'trying to open the .fe file'

            print 'The number of time of opening .fe file in Surface Evolver is', n

            x = E.run_command("tilt_angle:={a}".format(a=tangle))                                                   
            print 'start'
            E.run_command('print "Pretreatment of droplet evolution"')
            E.run_command('{g 1500; V 100;}')
            E.run_command('')
            while True:
                try:
                    averageLength = float(E.run_command('averageLength'))
                    print 'Got averageLength'
                    E.run_command('')
                    break
                except:
                    print 'Trying to get Average Length'

            print 'First Average Length is ', averageLength
            E.run_command('')

            while (averageLength > 0.10 * lambda1):
                E.run_command('{g 10; V 100;}')
                E.run_command('rz')
                E.run_command('')
                while True:
                    try:
                        averageLength = float(E.run_command('averageLength'))
                        print 'Got averageLength in the while (rz) check loop'
                        E.run_command('')
                        break
                    except:
                        print 'Trying to get Average Length in the while (rz) check loop'
                print 'averageLength is ', averageLength
            E.run_command('')
            while True:
                try:
                    thetaright = float(E.run_command('car'))
                    print 'Got thetaright'
                    E.run_command('')
                    thetaleft = float(E.run_command('cal'))
                    print 'Got thetaleft'
                    E.run_command('')
                    e_0 = float(E.evolve())
                    print 'Got e_0 (initial energy value)'
                    E.run_command('')
                    break
                except:
                    print 'Trying to get e_0, thetaright and thetaleft before entering nested while loops'

            E.run_command('')

            print 'Theta right and Theta Left values are respectively ', thetaright, ' and ', thetaleft
            print 'Initial e_0 (initial energy for iteration energy) is ', e_0

            E.run_command('')

            delta_E =100
            while ((abs(delta_E) > 10**(-_nrg_condition_for_E) or abs(thetaright-thetaE)>1.00 or abs(thetaleft-thetaE)>1.00) and breaker == False):
                while True:
                    try:
                        e_0 = float(E.evolve())
                        print 'Got e_0 (initial energy value)'
                        E.run_command('')
                        break
                    except:
                        print 'Trying to get e_0 in the outer loop'
                no_of_g = 0
                print "inside outer loop"
                delta_e = 100
                while(abs(delta_e) > 10**(-_nrg_condition_for_e)):
                    no_of_g = no_of_g + 1

                    while True:
                        try:
                            e_1 = float(E.evolve())
                            print 'Got e_1 (current energy value)'
                            E.run_command('')
                            break
                        except:
                            print 'Trying to get e_0 in the inner loop'
                    E.run_command('')
                    print 'even out'
                        
                    E.run_command('even_out')

                    delta_e = float((e_1 - e_0)/(e_0))
                    E.run_command('')
                    while True:
                        try:
                            leftx = float(E.run_command('left_vx'))
                            print 'Got current Left vx'
                            E.run_command('')
                            rightx = float(E.run_command('right_vx'))
                            print 'Got current right vx'
                            E.run_command('')
                            break
                        except:
                            print 'Trying to get the current left and right vertex values'
                    E.run_command('')
                    rightdifference = rightx - initialrightx
                    leftdifference = leftx - initialleftx
                    print 'right difference (has to be greater than', lambda1, ')', rightdifference
                    print 'left difference (has to be greater than', lambda1, ')', leftdifference
                    print 'delta_e ', delta_e
                    print 'delta_E ', delta_E
                    print 'e_1 ', e_1
                    print 'e_0 ', e_0
                    print 'Tilt Angle ', tangle
                    print 'This is an left cal contact angle measurement  ', thetaleft
                    print 'This is an right car contact angle measurement  ', thetaright
                    if  rightdifference > lambda1 and leftdifference > lambda1:
                        print 'SLIDING CONDITION IS BROKEN'
                        print 'rightdifference', rightdifference
                        print 'leftdifference', leftdifference
                        print 'This is FINAL left contact angle measurement  ', thetaleft
                        print 'This is FINAL right contact angle measurement  ', thetaright
                        rightsliding = True
                        leftsliding = True
                        breaker = True
                        break
                    if (no_of_g % 10):
                        print "Entering the am_i_Close block"
                        check = E.run_command('am_I_close')
                        print 'Check is ', check
                        if (check == "True"):
                            print 'flag changed'
                            flagger = -1
                            breaker = True
                            break   
                    delta_e = float((e_1 - e_0)/(e_0))             #it is a negative value
                    if (delta_e > 0):
                        flagger=1
                        print "e_1", e_1
                        print "e_0", e_0
                        break

                    e_0 = e_1
                E_1 = e_1
                while True:
                    try:
                        thetaright = float(E.run_command('car'))
                        print 'Got thetaright under E_1 recording'
                        E.run_command('')
                        thetaleft = float(E.run_command('cal'))
                        print 'Got thetaleft under E_1 recording'
                        E.run_command('')
                        break
                    except:
                        print 'Trying to get thetaright and thetaleft before leaving outer nested while loop'

                E.run_command('')

                print 'This is an updated left cal contact angle measurement  ', thetaleft
                print 'This is an updated right car contact angle measurement  ', thetaright
                delta_E = (E_1 - E_0)/(E_0)
                if (delta_E > 0):
                    flagger=1
                    #breaker = True
                    #break
                E_0 = E_1
                print 'Current Energy Value ', E_1
                print 'e_0 ', e_0
                print 'e_1 ', e_1
                print 'E_0' , E_0
                print 'E_1' , E_1
                E.run_command('rz')

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
            '''
            while True:
                try:
                    E.run_command('dump "bo{a}ampperlam{b}"'.format(a = Bo, b = amplitude))
                    print 'GOT IT! the dump a file for troubleshooting'
                    break
                except:
                    print 'trying to get the dump file for troubleshooting'

            '''
            CAH = abs(carapp - calapp)
            
            print ' '
            print ' '
            print ' '
            print '     Tilt Angle is ', tangle
            print 'Right Contact Angle Difference is ', abs(thetaright-thetaE)
            print 'Left Contact Angle Difference is  ', abs(thetaleft-thetaE)
            print "     CARAPP IS ", carapp
            print "     CALAPP IS ", calapp
            print "     CAH IS ", CAH
            print 'breaker is ', breaker
            print 'condition 1 is ', abs(delta_E) > 10**(-_nrg_condition_for_E)
            print 'condition 2 is ', abs(thetaright-thetaE)>1.00
            print 'condition 3 is ', abs(thetaleft-thetaE)>1.00
            print 'Outer While Loop Condition is ', ((abs(delta_E) > 10**(-_nrg_condition_for_E) or abs(thetaright-thetaE)>1.00 or abs(thetaleft-thetaE)>1.00) and breaker == False)
            print 'delta_e is ', delta_e
            print 'delta_E is ', delta_E
            print 'e_1 is ', e_1
            print 'e_0 is ', e_0

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

            
            print "Preparing to close"
           


            #E.run_command('dump "firstdump{b}"'.format(a=ampperlam, b = Bo))
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
        

        print 'Wavelength ', lambda1
        print 'Amplitude ', amplitude
        print 'Tilt Angle of the Plane:  ', tangle
        print 'delta_E ', delta_E
        print 'delta_e ', delta_e
        print 'breaker', breaker
        print 'Number of times through the inner loop', no_of_g
        print 'CAH ', CAH  
        print 'right apparent contact angle ', carapp
        print 'left apparent contact angle ', calapp
        print 'right contact angle  ', thetaleft
        print 'left contact angle ', thetaright
        print 'Current Energy Value ', minimum_e_for_tilt_angle
        #E.run_command('dump "dump{a}"'.format(a=tangle)) 

        if (tangle >= 90):
            tangle = 90 

                
                      

        ''' Checking the contact angles with the thetaE values within 0.5 degrees; if not BREAK from entire program

        if (abs (float(E.run_command('car')) - thetaE) > 7.50):
            print 'car has broken thetaE'
            print 'CAR is ' + str(E.run_command('car'))
            print 'thetaE is ' + str(E.run_command('print thetaE'))
            completebreaker = True
            break
        if(abs (float(E.run_command('cal')) - thetaE) > 7.50):
            print 'cal has broken thetaE'
            print 'CAL is ' + str(E.run_command('cal'))
            print 'thetaE is ' + str(E.run_command('print thetaE'))
            completebreaker = True
            break '''

                
        

        '''Setting amp/lam as 0.1 '''
        
        '''
        
        txtFile = open('LogFileGraphBondNumber_Wavelength_Tangle.txt', 'a') #;remember to delete file everytime or adjust program, because it will keep on appending values to the textfile
        txtFile.write(str(Bo) + '             ')
        txtFile.write(str(lambda1) + '              ')
        txtFile.write(str(tangle) + '            ')
        txtFile.write(str(flagger) + '     \n ')

        txtFile = open('LogFileGraphBondNumber_Amplitude_Tangle.txt', 'a') #;remember to delete file everytime or adjust program, because it will keep on appending values to the textfile
        txtFile.write(str(Bo) + '               ')
        txtFile.write(str(amplitude) + '           ')
        txtFile.write(str(tangle) + '                ')
        txtFile.write(str(flagger) + '       \n ') '''

        

        '''Phase Diagram between bond no. 0 to 1.00 and ampperlam from 0 to 0.3; change amplitude from 0.01 to 0.30 with lambda being fixed at 1.0 '''
        '''
        txtFile = open('Lastbitofampperlamandtanlge.txt', 'a') #;remember to delete file everytime or adjust program, because it will keep on appending values to the textfile
        txtFile.write(str(Bo) + '         ')
        txtFile.write(str(ampperlam) + '              ')
        txtFile.write(str(tangle) + '            ')
        txtFile.write(str(flagger) + '\n ')         '''



        '''Phase Diagram between wavlength from 0 to 1.00 and ampperlam from 0.0 to 0.40; use ampperlam to find corresponding amplitude for each point '''

        txtFile = open('LogFileGraphWavelength_AmpLam_CAH.txt', 'a') #;remember to delete file everytime or adjust program, because it will keep on appending values to the textfile
        txtFile.write(str(lambda1) + '                         ')
        txtFile.write(str(ampperlam) + '                         ')
        txtFile.write(str(CAH) + '                                  ')
        txtFile.write(str(flagger) + '                      \n ')

        txtFile = open('LogFileGraphWavelength_AmpLam_Energy.txt', 'a') #;remember to delete file everytime or adjust program, because it will keep on appending values to the textfile
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


        #counter = counter + 1

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
