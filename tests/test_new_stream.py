# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 11:59:31 2016

@author: sam
"""

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

q = Queue.Queue()

t1 = threading.Thread(target=imu.returnIMU, args=(q,"2G","245DPS","4GAUSS"))
t2 = threading.Thread(target=imu.recordSerial, args=(q,'COM3'))
t2.start()
t1.start()        

f = open('newData.txt','w')

stopTime = 100
startTime = time.time()

while time.time() - startTime < stopTime:
    #useless loop
    #print q.qsize()
    if not q.empty():
        data = q.get()
        if '$' in data:
            data = data.split(',')
            print data[7],' ',data[11], ' ', data[13], ' ', data[18]
        f.write(str(data))
#t2.join()
#t1.join()        
f.close()