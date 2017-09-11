'''
    Author: Sam Schmidt
    Purpose: navigation functions
    Date: 02/06/15 - 03/01/15
    Dependencies: numpy, scipy, matplotlib, pandas, serial, time, os
'''
###imports
from __future__ import division #enables default float division
import numpy as np #used for signal processing
import scipy.signal as signal #used for signal processing
import matplotlib.pyplot as plt #used for plotting, mainly for debugging
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd #has better CSV capabilities than numpy
import serial #used to communicate with the sensors
import time #used for timing sections of code, and will be used to meter sensor communications
import Queue #used to communicate between threads
import imu_comm_new as imu #used to communicate with the IMU
import threading #used to make threads
from otherFunctions import * #populate name space with otherFunctions code
from ins_sim import * #populate name space with ins_sim code
# import os #could be used for parallel threads


#################################################################################################################
######################## Entering the INS functions #############################################################
#################################################################################################################

##IIR filter class
#Filter used to remove high frequency (anything above 2 Hz) noise from the accelerometer
#currently configured for a Chebychev Type 2 IIR lowpass filter
class digitalFilter(object):
    def __init__(self,fc = [2,3],fs = 100,passGain = 1, stopGain = 70, analog = False, ftype = 'cheby2'):
        #filter coefficients
        self.b,self.a = signal.iirdesign(2*fc[0]/fs,2*fc[1]/fs,passGain,stopGain,analog,ftype)
        #filter memory
        self.zi = np.zeros(max([len(self.b),len(self.a)]))
        pass

    def iirFilt(self,x,reset = False):
        #reset is here incase somthing bad happens elsewhere
        if reset:
            self.zi = np.zeros(max([len(self.b),len(self.a)]))
        #implement direct form transpose II
        self.zi[0] = x - (self.a[1:]*self.zi[1:]).sum()
        y = (self.b*self.zi).sum()
        self.zi = np.roll(self.zi,1)
        self.zi[0] = x
        return y

#class used to navigate
#handles talking to the sensors and implements all navigation equations, it also controls the Kalman filter
class kalmanControl(object):
    def __init__(self,sim='lin',leverArm = np.zeros((3,1)),measNoise = np.zeros((15,15))):
        #set up the accel filters
#        self.axFilter = digitalFilter()
#        self.ayFilter = digitalFilter()
        # self.azFilter = digitalFilter()
#        ###############################
#        #test class for the imu
#        self.sim = sim
#        self.imu_comm = imu_sim(convert2decimalDeg(3.852355371e+03)*np.pi/180.,
#        convert2decimalDeg(1.044188916e+04)*np.pi/180.,convert2decimalDeg(3.952355371e+03)*np.pi/180.,\
#        convert2decimalDeg(1.044188916e+04)*np.pi/180.,alti = 2004,sim = self.sim,simTime = 3600, \
#        updateRate = 100,noise = True,bias = True)
#        print self.imu_comm.biasA
#        print self.imu_comm.biasG
#        ###############################
        ################# Constants #########################
        #sample rate
        self.fs = 100.00
        self.ts = 1/self.fs
        #degrees to radian conversion
        self.d2r = np.pi/180.
        #Earth's radius
        self.R0 = 6378137 #metres
        #eccentricity of the WGS84 ellipsoid
        self.eccenSquare = 6.69437999014e-3
        #Earth's rotational velocity
        self.omegaE = 7.2921159e-5
        #flattening of the wgs84 ellipsoid
        self.flat = 1/298.257223563
        #Earth's gravitational constant
        self.mu = 3.986004418e14
        ################# Kalman filter matrices ###############################
        #helper matrices for the measurement noise matrix
        A = np.vstack((np.hstack((.0000008*np.eye(3),np.zeros((3,3)))),np.hstack((np.zeros((3,3)),1e-3*np.eye(3))),np.zeros((3,6))))
        D = np.vstack((np.hstack((np.zeros((3,3)),1e-14*np.eye(3),np.zeros((3,3)))),np.hstack((np.zeros((3,6)),5e-8*np.eye(3)))))
        #measurement noise matrix: commonly called Q
        self.measNoise = np.vstack((np.hstack((A,np.zeros((9,9)))),np.hstack((np.zeros((6,6)),D))))*self.ts
        #noise covarience of the aiding system: commonly called R
        self.noiseCovar = np.array([[2.5066e-10,0,0,0,0,0],
                                    [0,2.9692e-10,0,0,0,0],
                                    [0,0,.30849e-5,0,0,0],
                                    [0,0,0,1e-3,0,0],
                                    [0,0,0,0,1e-3,0],
                                    [0,0,0,0,0,1e-3]])
        #curvilinear position scaling matrix
        self.cLinMat = np.array([[10e3,0,0],
                                 [0,10e3,0],
                                 [0,0,1]]) #100 on element 3,3 worked well
        #lever arm from the IMU to the aiding system, default is a column of zeros to simplify stuff
        self.leverArm = leverArm
        #helper matrices for the measurement matrix
        A = np.vstack((np.zeros((3,6)),np.hstack((np.zeros((3,3)),-1*np.eye(3)))))
        B = np.vstack((-self.cLinMat,np.zeros((3,3))))
        #measurement matrix commonly called: H
        self.measMatrix = np.hstack((A,B,np.zeros((6,6))))
        self.biasA = np.zeros((3,1))
        self.biasG = np.zeros(3) #might need to add in earths rotation
        self.biasM = np.zeros(3)
        #--------------------------Setup sensor communications----------------#
        '''
        Sensors comm is handled by through python threads, see information on
        python's global interpreter lock (GIL) if the kalman filter starts
        bottle necking things
        '''        
                
        #init imu
        check = imu.init_imu()
        assert(check == 0)
        
        #configure accel for default op
        fullScaleA = "2G"
        rate = "100HZ"
        check = imu.configAccel(rate, fullScaleA)
        assert(check == 0)
        
        #configure Mag for default op
        fullScaleM = "4GAUSS"
        rate = "100HZ"
        check = imu.configMag(rate, fullScaleM)
        assert(check == 0)
        
        #configure gyro for default op
        fullScaleG = "245DPS"
        check = imu.configGyro(fullScaleG)
        assert(check == 0)
        
        self.q = Queue.Queue()
        
        t1 = threading.Thread(target=imu.returnIMU, args=(self.q,"2G","245DPS","4GAUSS"))
        t2 = threading.Thread(target=imu.recordSerial, args=(self.q,'COM3'))
        t2.start()
        t1.start()        
        
#        for item in source():
#            q.put(item)
        
        pass

    def getIMUGPS(self,gravConst = 0):
            #this function will pull a value from the IMU, filter it and return stuff
            #inputs:
            #-gravConst is the current estimate of the positions gravity
            #grab the data, assume imu data is structured as (ax,ay,az,gx,gy,gz,mx,my,mz)
            #big impressive routine
            if not self.q.empty():
                data = self.q.get()
                #check GPS
                if '$' in data:
                    try:
                        #GPS
                        data = data.split(',')
                        self.gpsPos = np.array([[float(data[11])],[float(data[13])],[float(data[18])]])
                        self.gpsVelNed = float(data[7])*np.array([[1],[1],[1]])*(1000/3600.)
                    except:
                        print 'Bad Data! noob'
                    #Then look for new IMU
                data = self.q.get()
                while '$' in data:
                    data = self.q.get()
    
            #unpack data
            self.accelMeas = np.asarray([[data[0]],[data[1]],[data[2]]])
            self.rotMeas = np.asarray([data[3],data[4],data[5]])
            self.magMeas = np.asarray([data[6],data[7],data[8]])
            #scale accel by g
            self.accelMeas *= 9.8
            #remove accel Bias
            self.accelMeas -= self.biasA
            #pass accel through the filter
            # self.accelMeas[0] = self.axFilter.iirFilt(self.accelMeas[0])
            # self.accelMeas[1] = self.ayFilter.iirFilt(self.accelMeas[1])
            # self.accelMeas[2] = self.azFilter.iirFilt(self.accelMeas[2])
            #convert rotaions from dps to rad/s
            self.rotMeas *= self.d2r
            #remove gyro bias
            self.rotMeas -= self.biasG
            #store the rotaions as a skew sym matrix
            self.rotMeasSkew = np.array([[0,-self.rotMeas[2],self.rotMeas[1]],[self.rotMeas[2],0,self.rotMeas[0]],[-self.rotMeas[1],self.rotMeas[0],0]])
            #get magnitude of gyro
            self.magRotMeas = np.sqrt(np.dot(self.rotMeas,self.rotMeas))*self.ts
            #-get GPS measurements
    #        self.gpsPos = np.array([[imu[9]],[imu[10]],[imu[11]]])
    #        self.gpsVelNed = imu[12]*np.array([[1],[1],[1]])#*np.array([[np.sin(self.state[2])*np.cos(self.state[1])],[np.cos(self.state[2])*np.cos(self.state[1])],[np.sin(self.state[1])]])


    def startup(self):
        #this function will determine the starting orientation, heading, and position
        #it will also instansiate the kalman filter class
        #init latitude, longitude, and altitude
        self.lat = convert2decimalDeg(3.852355371e+03)*self.d2r
        self.lon = convert2decimalDeg(1.044188916e+04)*self.d2r
        self.alti = 2004

        #get grave constant
        self.g0 = 9.7803253359*((1 + 0.001931853*(np.sin(self.lat))**2)/np.sqrt(1 - self.eccenSquare*(np.sin(self.lat))**2))

        #init velocity
        if self.sim == 'lin' or self.sim == 'stat':
            self.velNed = np.array([[0],[0],[0]])
        else:
            self.velNed = np.array([[0],[500*(2*np.pi/3600.)],[0]])

        #init state vector
        self.state = 10e-18*np.ones(15)
        #init state transition matrix
        self.PHI = 10e-24*np.eye(15)
        #init covarience matrix
        self.covar = 10e-24*np.eye(15)
        #get the earths rotation in NED coordinates
        self.earthRotINNed = (self.omegaE)*np.array([[0,np.sin(self.lat),0],
                                                     [-np.sin(self.lat),0,-np.cos(self.lat)],
                                                     [0,np.cos(self.lat),0]])
        #estimate the bias values, use 30 seconds of data to do the estimation
        biasEstimatorA = np.zeros((3,1))
        biasEstimatorG = np.zeros(3)
        biasEstimatorM = np.zeros(3)
#        print biasEstimator + self.imu_comm.sample(self.lat,self.alti)
        sampleCount = 0
        while(not sampleCount == 3000):
            self.getIMUGPS()
            biasEstimatorA += self.accelMeas
            biasEstimatorG += self.rotMeas
            biasEstimatorM += self.magMeas
            sampleCount += 1
        biasEstimator /= 3000

        #init bias values
        self.biasG = np.asarray([biasEstimator[3],biasEstimator[4],biasEstimator[5]])
        self.biasA = np.asarray([[biasEstimator[0]],[biasEstimator[1]],[biasEstimator[2]]])
        print 'estimated gyro bias: ',self.biasG
        print 'estimated Accel bias: ',self.biasA
        #init body to nav frame rotation matrix
        self.rotB2N = self.earthRotINNed
        #calculate the earths radius
        self.RE = self.R0/(np.sqrt(1-self.eccenSquare*(np.sin(self.lat))**2)) #transverse radius
        self.RN = (self.R0*(1-self.eccenSquare))/((1 - self.eccenSquare*(np.sin(self.lat))**2)**(1.5))
        #get the rotation matrix to resolve to the difference between a earth centered frame and the local NED frame
        rotEarthandNed = np.array([[self.velNed[1]/(self.RE + self.alti)],[-self.velNed[0]/(self.RN + self.alti)],[-np.tan(self.lat)*self.velNed[1]/(self.RE + self.alti)]])
        self.rotEarthandNedSkew = np.array([[0,-rotEarthandNed[2],rotEarthandNed[1]],[rotEarthandNed[2],0,-rotEarthandNed[0]],[-rotEarthandNed[1],rotEarthandNed[0],0]])
        pass


    def getRot(self):
        #method to update the rotation matrix from the gyro readings
        #this is the rotation from the body frame to the navigation frame
        self.prevRotB2N = self.rotB2N
        # print self.rotB2N
        #updated gyro matrix
        rotBP2BM = np.eye(3) + (np.sin(self.magRotMeas)/self.magRotMeas)*self.rotMeasSkew + ((1-np.cos(self.magRotMeas))/self.magRotMeas)*(np.dot(self.rotMeasSkew,self.rotMeasSkew))
        # print self.magRotMeas
        #matrix to account for the Earth's rotation in a local reference frame
        self.earthRotINNed = (self.omegaE)*np.array([[0,np.sin(self.lat),0],[-np.sin(self.lat),0,-np.cos(self.lat)],[0,np.cos(self.lat),0]])
        #calculate teh earfs radius
        self.RE = self.R0/(np.sqrt(1-self.eccenSquare*(np.sin(self.lat))**2)) #transverse radius
        self.RN = (self.R0*(1-self.eccenSquare))/((1 - self.eccenSquare*(np.sin(self.lat))**2)**(1.5))
        rotEarthandNed = np.array([[self.velNed[1]/(self.RE + self.alti)],[-self.velNed[0]/(self.RN + self.alti)],[-np.tan(self.lat)*self.velNed[1]/(self.RE + self.alti)]])
        self.rotEarthandNedSkew = np.array([[0,-rotEarthandNed[2],rotEarthandNed[1]],[rotEarthandNed[2],0,-rotEarthandNed[0]],[-rotEarthandNed[1],rotEarthandNed[0],0]])
        #update rotation matrix
        self.rotB2N = np.dot(self.prevRotB2N,rotBP2BM) - np.dot((self.earthRotINNed + self.rotEarthandNedSkew),self.prevRotB2N*self.ts)


    def getVel(self):
        #method to get the velocity from the accels
        self.prevVelNed = self.velNed
        #update rotation matrix
        rotBB2BM = np.eye(3) + ((1-np.cos(self.magRotMeas))/(self.magRotMeas)**2)*self.rotMeasSkew + (1/(self.magRotMeas)**2)*(1 - (np.sin(self.magRotMeas)/self.magRotMeas))*(np.dot(self.rotMeasSkew,self.rotMeasSkew))
        rotB2NBar = np.dot(self.prevRotB2N,rotBB2BM) - 0.5*np.dot((self.earthRotINNed + self.rotEarthandNedSkew),self.prevRotB2N*self.ts)
        rotAccel = np.dot(rotB2NBar,self.accelMeas) #NOTE this is commented out, since I don't simulate gyro stuff....
        #rotAccel = self.accelMeas
        #get teh gravsssz
        self.g0 = 9.7803253359*((1 + 0.001931853*(np.sin(self.lat))**2)/np.sqrt(1 - self.eccenSquare*(np.sin(self.lat))**2))
        gN = (-8.08e-9)*self.alti*np.sin(2*self.lat)
        gD = self.g0*(1 - (1 + self.flat*(1-2*(np.sin(self.lat)**2)) + ((self.omegaE*self.R0)**2)*(6356752.3142)/self.mu)*(2*self.alti/self.R0) + (3/(self.R0**2))*(self.alti**2))
        self.gVec = np.array([[gN],[0],[gD]])
        #update velocity
        self.velNed = self.prevVelNed + (rotAccel + self.gVec - np.dot((self.earthRotINNed + 2*self.rotEarthandNedSkew),self.prevVelNed))*self.ts


    def getPos(self):
        #this method will get the current position of the system resolved in geodetic coordinates
        #get altitude
        newAlti = self.alti - (self.ts/2.)*(self.prevVelNed[2] + self.velNed[2])
        #get latitude
        #-first calculate the meridian radius of curvature
        RNPrevLat = (self.R0*(1-self.eccenSquare))/((1 - self.eccenSquare*(np.sin(self.lat))**2)**(1.5))
        newLat = self.lat + (self.ts/2.)*((self.prevVelNed[0]/(RNPrevLat + self.alti)) + (self.velNed[0]/(RNPrevLat + newAlti)))
        #get new longitude
        #-first calculate the new and old transverse radius
        REPrevLat = self.R0/(np.sqrt(1-self.eccenSquare*(np.sin(self.lat))**2))
        RENewLat = self.R0/(np.sqrt(1-self.eccenSquare*(np.sin(newLat[0]))**2))
        newLon = self.lon + (self.ts/2.)*((self.prevVelNed[1]/((REPrevLat + self.alti)*np.cos(self.lat))) + (self.velNed[1]/((RENewLat + newAlti[0])*np.cos(newLat[0]))))
        #update position
        self.alti = newAlti[0]
        self.lat = newLat[0]
        self.lon = newLon[0]

    '''
    ------------------------Begin Kalman Filter---------------------------------
    '''
    ################# calculate the transition matrix ##########################
    def buildStateTranMatrix(self):
        #need to recalulate the RE and RN with the newest corrected latitude, longitude, and altitude
        self.prevPHI = self.PHI
        #get submatrices
        self.f11 = -(self.earthRotINNed + self.rotEarthandNedSkew)

        self.f12 = np.array([[0,-1/(self.RE + self.alti),0],
                             [1/(self.RE + self.alti),0,0],
                             [0,np.tan(self.lat)/(self.RE + self.alti),0]])

        self.f13 = np.array([[self.omegaE*np.sin(self.lat),0,self.velNed[1]/(self.RE + self.alti)**2],
                             [0,0,-self.velNed[0][0]/(self.RN + self.alti)**2],
                             [self.omegaE*np.cos(self.lat) + self.velNed[1][0]/((self.RE + self.alti)*(np.cos(self.lat)**2)),0,-np.tan(self.lat)*self.velNed[1][0]/(self.RE + self.alti)**2]])

        self.f21 = -1*np.array([[0,-np.dot(self.rotB2N[2,:],self.accelMeas),np.dot(self.rotB2N[1,:],self.accelMeas)],
                              [np.dot(self.rotB2N[2,:],self.accelMeas),0,-np.dot(self.rotB2N[0,:],self.accelMeas)],
                              [-np.dot(self.rotB2N[1,:],self.accelMeas),np.dot(self.rotB2N[0,:],self.accelMeas),0]])

        self.f22 = np.array([[self.velNed[2][0]/(self.RN + self.alti),-2*np.tan(self.lat)*self.velNed[1][0]/(self.RE + self.alti) - 2*self.omegaE*np.sin(self.lat),self.velNed[0][0]/(self.RN + self.alti)],
                             [np.tan(self.lat)*self.velNed[1][0]/(self.RE + self.alti) + 2*self.omegaE*np.sin(self.lat),(self.velNed[0][0]*np.tan(self.lat) + self.velNed[2][0])/(self.RE + self.alti),self.velNed[1][0]/(self.RE + self.alti) - 2*self.omegaE*np.sin(self.lat)],
                             [-2*self.velNed[0][0]*(self.RN + self.alti),-2*self.velNed[1][0]*(self.RE + self.alti) - 2*self.omegaE*np.cos(self.lat),0]])

        self.f23 = np.array([[-(((self.velNed[1][0]**2)*((1/np.cos(self.lat))**2))/(self.RE + self.alti)) - 2*self.velNed[1][0]*self.omegaE*np.cos(self.lat),0,((((self.velNed[1][0]**2)*((np.tan(self.lat))))/((self.RE + self.alti)**2))) - ((self.velNed[0][0]*self.velNed[2][0])/(((self.RN + self.alti)**2)))],
                             [(((self.velNed[1][0]*self.velNed[0][0])*((1/np.cos(self.lat))**2))/(self.RE + self.alti)) + 2*self.velNed[0][0]*self.omegaE*np.cos(self.lat) - 2*self.velNed[2][0]*self.omegaE*np.sin(self.lat),0,-(((self.velNed[0][0]*self.velNed[1][0])*((np.tan(self.lat))) + (self.velNed[1][0]*self.velNed[2][0])))/((self.RE + self.alti)**2)],
                             [2*self.velNed[1][0]*self.omegaE*np.sin(self.lat),0,(self.velNed[1][0]**2)/((self.RE + self.alti)**2) + (self.velNed[0][0]**2)/((self.RN + self.alti)**2) - (2*self.g0)/(self.RE*np.sqrt(1 + (-2*np.sqrt(self.eccenSquare) + self.eccenSquare)*(np.sin(self.lat))**2))]])

        self.f32 = np.array([[1/(self.RN + self.alti),0,0],
                             [0,1/((self.RE + self.alti)*np.cos(self.lat)),0],
                             [0,0,-1]])

        self.f33 = np.array([[0,0,self.velNed[0][0]/(self.RN + self.alti)**2],
                             [self.velNed[1][0]*np.sin(self.lat)/((self.RE + self.alti)*(np.cos(self.lat)**2)),0,-self.velNed[1][0]*np.sin(self.lat)/(((self.RE + self.alti)**2)*(np.cos(self.lat)))],
                             [0,0,0]])

        #print self.f22
        #I feel like there will be a lot of dimension problems here....
        self.PHI = np.eye(15) + (np.vstack((np.hstack((self.f11,self.f12,self.f13,np.zeros((3,3)),self.rotB2N)),
                                          np.hstack((self.f21,self.f22,self.f23,self.rotB2N,np.zeros((3,3)))),
                                          np.hstack((np.zeros((3,3)),self.f32,self.f33,np.zeros((3,3)),np.zeros((3,3)))),
                                          np.zeros((6,15)))))*self.ts

    ################# calculate covarience #####################################
    def getCovar(self):
        #self.prevCovar = self.covar
        self.covar = np.dot(self.prevPHI,np.dot(self.covar,self.prevPHI.T)) + self.measNoise

    ################# Update states ############################################
    def update(self):
        #map curve pertubations from Cartesian to curvilinear
        T = np.array([[1/(self.RN + self.alti),0,0],[0,1/((self.RE + self.alti)*np.cos(self.lat)),0],[0,0,-1]])
        #update states
        newState = np.dot(self.PHI,self.state)
        #need to think about the implementation of this, might need a least squares solution
        self.kalmanGain = np.dot(self.covar,np.dot(self.measMatrix.T,np.linalg.inv(np.dot(self.measMatrix,np.dot(self.covar,self.measMatrix.T)) + self.noiseCovar)))
        #form the innovation
        #print self.gpsPos
        self.measInnov = np.vstack((np.dot(self.cLinMat,(self.gpsPos - np.array([[self.lat],[self.lon],[self.alti]]) - np.dot(T,np.dot(self.rotB2N,self.leverArm))))
                                  ,self.gpsVelNed-self.velNed - np.dot(self.rotB2N,np.dot(self.rotMeasSkew,self.leverArm)) + np.dot(self.earthRotINNed,np.dot(self.rotB2N,self.leverArm)))).flatten()
        #update the states
        #self.state = self.state + np.dot(self.kalmanGain,self.measInnov)
        self.state = newState + np.dot(self.kalmanGain,self.measInnov)
        #update the covarience
        self.covar = np.dot((np.eye(15) - np.dot(self.kalmanGain,self.measMatrix)),self.covar)

    '''
    -------------------------End Kalman Filter---------------------------------
    '''

    #this method is called by navigate to estimate the errors
    def estimateError(self):
        self.buildStateTranMatrix()
        self.getCovar()
        self.update()


    """
    ''Public'' methods, these are called to generate the navigation solution
    """

    #this method brings everything together and navigates
    def navigate(self):
        #get imu and gps data
        self.getIMUGPS()
        #advance the kalman filter
        self.estimateError()
        #print self.state
        #get new rotation
        self.getRot()
        #get new velocity
        self.getVel()
        #get new position
        self.getPos()
        #correct the measurements
        self.lat -= self.state[6]
        self.lon -= self.state[7]
        self.alti -= self.state[8]
        self.velNed[0] -= self.state[3]
        self.velNed[1] -= self.state[4]
        self.velNed[2] -= self.state[5]
        self.rotB2N = np.dot((np.eye(3) - np.array([[0,-self.state[2],self.state[1]],
                                                    [self.state[2],0,-self.state[0]],
                                                    [-self.state[1],self.state[0],0]])),self.rotB2N)
        self.biasA -= np.asarray([[self.state[9]],[self.state[10]],[self.state[11]]])
        self.biasG -= np.asarray([self.state[12],self.state[13],self.state[14]])
        #zero the states
        self.prevState = self.state
        self.state = np.zeros_like(self.state)
        #recalulate the Earths radius
        self.RE = self.R0/(np.sqrt(1-self.eccenSquare*(np.sin(self.lat))**2)) #transverse radius
        self.RN = (self.R0*(1-self.eccenSquare))/((1 - self.eccenSquare*(np.sin(self.lat))**2)**(1.5))
        #return
        return np.hstack((self.lat,self.lon,self.alti,self.velNed.T.flatten(),self.prevState.flatten())).flatten()
##end kalmonControl class