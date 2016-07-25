# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 10:25:03 2016

@author: sam
"""

###imports
from __future__ import division #enables default float division
import numpy as np #used for signal processing
import scipy.signal as signal #used for signal processing
import matplotlib.pyplot as plt #used for plotting, mainly for debugging
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd #has better CSV capabilities than numpy
import serial #used to communicate with the sensors
import time #used for timing sections of code, and will be used to meter sensor communications
from otherFunctions import * #populate name space with otherFunctions code
from navFunctions import * #populate name space with navFunctions code

##debug function, reacts to nan values returned from the Kalman filter
def nanDebug(writeValues):
    np.savetxt('valueDump.csv', writeValues, delimiter=',')
    raise ValueError('NaN values returned from the filter, outputs written to valueDump.csv')

#################################################################################################################
##################################### Entering test #############################################################
#################################################################################################################

print 'Running test...'

navTime = np.arange(3600*100. - 3000.)
lat     = np.zeros(len(navTime)+1)
lon     = np.zeros(len(navTime)+1)
alti    = np.zeros(len(navTime)+1)
vel     = np.zeros((len(navTime)+1,3))

roll    = np.zeros(len(navTime)+1)
pitch   = np.zeros(len(navTime)+1)
yaw     = np.zeros(len(navTime)+1)
velNErr = np.zeros((len(navTime)+1))
velEErr = np.zeros(len(navTime)+1)
velDErr = np.zeros(len(navTime)+1)
latErr  = np.zeros(len(navTime)+1)
lonErr  = np.zeros((len(navTime)+1))
altiErr = np.zeros((len(navTime)+1))
biasGx  = np.zeros((len(navTime)+1))
biasGy  = np.zeros((len(navTime)+1))
biasGz  = np.zeros((len(navTime)+1))
biasAx  = np.zeros((len(navTime)+1))
biasAy  = np.zeros((len(navTime)+1))
biasAz  = np.zeros((len(navTime)+1))

lat[0] = convert2decimalDeg(3.852355371e+03)*(np.pi/180)
lon[0] = convert2decimalDeg(1.044188916e+04)*(np.pi/180)
alti[0] = 2004
control = kalmanControl(sim='stat')
control.startup()
for i in navTime.astype(int):
    #check for samples
    lat[i+1],lon[i+1],alti[i+1],vel[i+1,0],vel[i+1,1],vel[i+1,2], \
    roll[i+1],pitch[i+1],yaw[i+1],velNErr[i+1],velEErr[i+1],velDErr[i+1], \
    latErr[i+1],lonErr[i+1],altiErr[i+1],biasAx[i+1],biasAy[i+1],biasAz[i+1], \
    biasGx[i+1],biasGy[i+1],biasGz[i+1] = control.navigate()

    if np.isnan(vel).any():
        writeValues = np.column_stack((lat,lon,alti,vel,roll,pitch,yaw,velNErr,velEErr,velDErr,latErr,lonErr, \
        altiErr,biasAx,biasAy,biasAz,biasGx,biasGy,biasGz))
        nanDebug(writeValues)
    if i % 10000 == 0:
        print 'GPS position in deg:\n', control.gpsPos*180./np.pi
        print
        print 'GPS velocity in NED:\n',control.gpsVelNed
        print
        print 'NED velocity:\n' ,control.velNed
        print
        print 'Measured accelartion:\n', control.accelMeas
        print
        print 'Accel bias: ',control.prevState[12],' ',control.prevState[13],' ',control.prevState[14]
        print
        print 'Gyro bias: ',control.prevState[9],' ',control.prevState[10],' ',control.prevState[11]
        print

lat *= 180/np.pi
lon *= 180/np.pi

# lat = convert2DegMM(lat)
# lon = convert2DegMM(lon)
print np.min(vel[:,1])
print np.max(vel[:,1])
print lat[0]
print lon[0]
print lat[-1]
print lon[-1]
print alti
print 'Error between simulated end latitude and actual latitude is: ', lat[-1] - convert2decimalDeg(3.952355371e+03)
print 'Error between simulated end longitude and actual longitude is: ', lon[-1] - convert2decimalDeg(1.044188916e+04)
latTrack = np.linspace(convert2decimalDeg(3.852355371e+03),convert2decimalDeg(3.952355371e+03),len(lat))
print 'Mean sqruared error for latitude is: ', np.mean((latTrack - lat)**2)/float(len(lat))
# fid = open('simCircNoiseBias.csv','a')
# data = np.column_stack((lat.T.flatten(),lon.T.flatten(),alti.T.flatten()))
# np.savetxt(fid,data,delimiter = ',')
# fid.close()
# plt.figure()
# plt.plot(np.arange(len(lat)),lat)
# plt.figure()
# plt.plot(np.arange(len(lon)),lon)
# plt.figure()
# plt.plot(np.arange(len(alti)),alti)
plt.figure()
# plt.plot(lon,lat)
plt.plot(lon[0],lat[0],'rs')
#plot the 'true' track here
#    plt.plot(lon[0]*np.ones(len(lat)),lat,'-r^')
plt.plot(lon[-1],lat[-1],'r*')
plt.plot(lon,lat)
plt.title('Plot of Latitude and Longitude: No measurement error')
plt.xlabel('Longitude (deg)')
plt.ylabel('Latitude (deg)')
plt.legend(('start','end','track'))
plt.xlim([np.min(lon) - 100e-3,np.max(lon) + 100e-3])
plt.ylim([np.min(lat) - 1000e-3,np.max(lat) + 1000e-3])
plt.grid()
plt.figure()
plt.plot(np.arange(len(lat)),lat)
plt.title('latitude')
plt.ylabel('latitude')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(lon)),lon)
plt.title('longitude')
plt.ylabel('longitude')
plt.xlabel('n')
plt.figure()
plt.plot(vel[:,1],vel[:,0])
plt.title('Plots of velocity in the north and east direction')
plt.xlabel('east')
plt.ylabel('north')
plt.figure()
plt.plot(np.arange(len(vel[:,0])),vel[:,0])
plt.title('Velocity in the north direction')
plt.ylabel('north')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(vel[:,1])),vel[:,1])
plt.title('Velocity in the east direction')
plt.ylabel('east')
plt.xlabel('n')
##################plot errors#############################
plt.figure()
plt.plot(np.arange(len(roll)),roll)
plt.title('Roll')
plt.ylabel('rad')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),pitch)
plt.title('Pitch')
plt.ylabel('rad')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),yaw)
plt.title('Yaw')
plt.ylabel('rad')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),velNErr)
plt.title('North Velocity Error')
plt.ylabel('m/s')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),velEErr)
plt.title('East Velocity Error')
plt.ylabel('m/s')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),velDErr)
plt.title('Down Velocity Error')
plt.ylabel('m/s')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),180*latErr/np.pi)
plt.title('Latitude Error')
plt.ylabel('deg')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),180*lonErr/np.pi)
plt.title('Longitude Error')
plt.ylabel('deg')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),altiErr)
plt.title('Altitude Error')
plt.ylabel('m')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),biasAx)
plt.title('Accel x bias')
plt.ylabel('m/s^2')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),biasAy)
plt.title('Accel y bias')
plt.ylabel('m/s^2')
plt.xlabel('n')
plt.figure()
plt.plot(np.arange(len(roll)),biasAz)
plt.title('Accel z bias')
plt.ylabel('m/s^2')
plt.xlabel('n')
##show plots
plt.show()