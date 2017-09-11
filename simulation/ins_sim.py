# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 10:22:31 2016

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

#class to simulate IMU readings
class imu_sim(object):
    def __init__(self,latStart,lonStart,latEnd,lonEnd,alti = 2004,sim = 'lin',simTime = 3600,updateRate = 100,noise = False,bias = False):
        #possible values for sim are 'lin', 'circ' or 'stat'
        if(sim != 'lin') and (sim != 'circ') and (sim != 'stat'):
            raise ValueError('sim can only be lin, circ, or stat')
        self.sim          = sim
        self.noiseFlag    = noise
        self.sampleTime   = 0
        self.ts           = 1/float(updateRate)
        self.Tlap         = simTime              # [s] lap time
        self.sLat         = np.sin(latStart)
        self.cLat         = np.cos(latStart)
        self.omegaEarth   = 0.7292115e-4*np.array([self.cLat,0,-self.sLat]) # Earthrate in NED coordinate
        self.lat          = latStart
        self.alti         = alti
        #Earth's radius
        self.R0 = 6378137 #metres
        #eccentricity of the WGS84 ellipsoid
        self.eccenSquare = 6.69437999014e-3
        #flattening of the wgs84 ellipsoid
        self.flat = 1/298.257223563
        self.omegaE   = 0.7292115e-4
        #Earth's gravitational constant
        self.mu = 3.986004418e14

        #handle the simulation type
        if sim == 'circ':
            #-Assume latEnd and lonEnd are points on the circle parallel to the starting latitude and
            #use Haversine's formula to find the radius of the circle
            self.circRadius = 2*self.R0*np.arcsin(np.sqrt(np.sin((latEnd-latStart)/2.)**2 +
                                                          np.cos(latEnd)*np.cos(latStart)*np.sin((lonEnd-lonStart)/2.)**2))
            self.curviLinR  = (latEnd-latStart)/2.
        else:
            self.gpsScale = ((latEnd - latStart)/(float(self.Tlap)))
        #get IMU bias values
        if bias:
            self.biasA   = np.sqrt(1e-2)*np.array([np.random.randn(1)[0],np.random.randn(1)[0],np.random.randn(1)[0]])
            self.biasG   = np.sqrt(1e-1)*np.array([np.random.randn(1)[0],np.random.randn(1)[0],np.random.randn(1)[0]])
            self.biasMag = np.sqrt(1e-1)*np.array([np.random.randn(1)[0],np.random.randn(1)[0],np.random.randn(1)[0]])
        else:
            self.biasA   = 0
            self.biasG   = 0
            self.biasMag = 0

        #init GPS reading
        self.gpsData = np.hstack((float(latStart),float(lonStart),float(alti),0.0))
        pass

    def accel(self):
        g0 = 9.7803253359*((1 + 0.001931853*(np.sin(self.lat))**2)/np.sqrt(1 - self.eccenSquare*(np.sin(self.lat))**2))
        gN = (-8.08e-9)*self.alti*np.sin(2*self.lat)
        gD = g0*(1 - (1 + self.flat*(1-2*(np.sin(self.lat)**2)) + ((self.omegaE*self.R0)**2)*(6356752.3142)/self.mu)*
                 (2*self.alti/self.R0) + (3/(self.R0**2))*(self.alti**2))
        if self.sim == 'lin':
            if self.sampleTime < 100:
                #accelarate north at a rate of .1 m/s^2
                self.accelMeas = np.array([10/100.-gN,0,-gD]) + self.noiseA + self.biasA
                #accelarate east at a rate of .03 m/s^2
                # self.accelMeas = np.array([-gN,3/100.,-gD])
                #accelarate north and east at a rate of .03 m/s^2
                # self.accelMeas = np.array([3/100.-gN,3/100.,-gD])
            else:
                self.accelMeas = np.array([-gN,0,-gD]) + self.noiseA + self.biasA
        elif self.sim == 'circ':
            self.accelMeas = np.array([-(self.circRadius*((2*np.pi/self.Tlap)**2))*np.cos((2*np.pi/self.Tlap)*self.sampleTime) - gN,
                                       -(self.circRadius*((2*np.pi/self.Tlap)**2))*np.sin((2*np.pi/self.Tlap)*self.sampleTime),-gD])
            self.accelMeas += self.noiseA + self.biasA
        else:
            #stationary
            self.accelMeas = np.array([-gN,0,-gD]) + self.noiseA + self.biasA
        pass

    def gyro(self):
        if self.sim == 'lin' or self.sim == 'stat':
            # self.RE = self.R0/(np.sqrt(1-self.eccenSquare*(np.sin(self.lat))**2)) #transverse radius
            # self.RN = (self.R0*(1-self.eccenSquare))/((1 - self.eccenSquare*(np.sin(self.lat))**2)**(1.5))
            # rotEarthandNed = np.array([[0/(self.RE + self.alti)],[-3/(self.RN + self.alti)],[-np.tan(self.lat)*self.velNed[1]/(self.RE + self.alti)]])
            # self.rotEarthandNedSkew = np.array([[0,-rotEarthandNed[2],rotEarthandNed[1]],[rotEarthandNed[2],0,-rotEarthandNed[0]],[-rotEarthandNed[1],rotEarthandNed[0],0]])
            self.gyroMeas = .00000000000001*self.omegaEarth #+ self.biasG#+ self.noiseG + self.biasG
        else:
            #The gyro measurments are going to be alot more complicated for the circular track
            self.gyroMeas = .00000000000001*self.omegaEarth #+ self.noiseG + self.biasG
        pass

    def mag(self):
        #probably won't touch the magnetometer, but a simplifying assumption is to only point north...
        self.magMeas = np.zeros(3) + self.noiseMag + self.biasMag
        pass

    def gps(self):
        if self.sim == 'lin':
            #use a linear interpolation
            if self.sampleTime < 100:
                self.gpsData[0] += self.gpsScale + self.noiseGPS[0]
                self.gpsData[1] = convert2decimalDeg(1.044188916e+04)*np.pi/180. +  self.noiseGPS[1]
                self.gpsData[3] += 10/100. + np.sqrt(7.3755e-03)*np.random.randn(1)[0] #increasing velocity
                self.gpsData[2] = 2004
            else:
                self.gpsData[0] += self.gpsScale + self.noiseGPS[0]
                self.gpsData[1] = convert2decimalDeg(1.044188916e+04)*np.pi/180. + self.noiseGPS[1]
                self.gpsData[3] = 10 + np.sqrt(7.3755e-03)*np.random.randn(1)[0]
                self.gpsData[2] = 2004
        elif self.sim == 'circ':
            self.gpsData[0] = self.curviLinR*np.sin(2*np.pi*self.sampleTime/self.Tlap) + self.noiseGPS[0]
            self.gpsData[1] = self.curviLinR*np.cos(2*np.pi*self.sampleTime/self.Tlap) + self.noiseGPS[0]
            self.gpsData[2] = self.alti
            self.gpsData[3] = self.curviLinR*2*np.pi/self.Tlap
        else:
            self.gpsData[0] = convert2decimalDeg(3.852355371e+03)*np.pi/180. + self.noiseGPS[0]
            self.gpsData[1] = convert2decimalDeg(1.044188916e+04)*np.pi/180. + self.noiseGPS[1]
            self.gpsData[2] = 2004 + np.sqrt(3.0849e+01)*np.random.randn(1)[0]
            self.gpsData[3] = 0 + np.sqrt(7.3755e-03)*np.random.randn(1)[0]
        pass

    def sample(self,lat,alti):
        if self.noiseFlag:
            self.noiseA   = np.array([np.sqrt(3.22026513267e-05)*np.random.randn(1)[0]
                                     ,np.sqrt(3.23082291483e-05)*np.random.randn(1)[0]
                                     ,np.sqrt(4.26628116972e-05)*np.random.randn(1)[0]])
            self.noiseG   = np.array([np.sqrt(0.080784168197)*np.random.randn(1)[0]
                                     ,np.sqrt(0.16687196328)*np.random.randn(1)[0]
                                     ,np.sqrt(0.088688910856)*np.random.randn(1)[0]])
            self.noiseGPS = np.array([np.sqrt(8.6486101901e-15)*np.random.randn(1)[0]
                                     ,np.sqrt(1.89549072165e-15)*np.random.randn(1)[0]])
            self.noiseMag = np.sqrt(.000005)*np.array([np.random.randn(1)[0],np.random.randn(1)[0],np.random.randn(1)[0]])
        else:
            self.noiseA    = np.zeros(3)
            self.noiseG    = np.zeros(3)
            self.noiseMag  = np.zeros(3)
            self.noiseGPS  = np.zeros(4)

        self.lat = lat
        self.alti = alti
        self.accel()
        self.gyro()
        self.mag()
        self.accelMeas /= 9.8
        self.gyroMeas = (180./np.pi)*self.gyroMeas
        if (self.sampleTime % 1 < 1e-5) or ((np.abs((self.sampleTime % 1) - 1)) < 1e-5):
            self.gps()
        self.sampleTime += self.ts
        return np.hstack((self.accelMeas,self.gyroMeas,self.magMeas,self.gpsData)).flatten()