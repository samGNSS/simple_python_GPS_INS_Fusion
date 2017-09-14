try:
    import Adafruit_GPIO.FT232H as FT232H
    FT232H.use_FT232H()
    ft232h = FT232H.FT232H()
except ImportError:
    print("Missing Adafruit library!")
import numpy as np
import serial
import time
from registers import XM,GYRO

def init_imu():
    #check if the Accel/Mag is there
    #-set i2c to point at accelerometer & Mag
    i2c = FT232H.I2CDevice(ft232h,XM.ADDRESS)
    check = i2c.readU8(XM.WHO_AM_I)
    if check != XM.WHO_AM_I_OK:
        print "ACCEL INIT FAIL..."
        return 1
    #check if the GYRO is there
    #-set i2c to point at XM
    i2c = FT232H.I2CDevice(ft232h,GYRO.ADDRESS)
    check = i2c.readU8(GYRO.WHO_AM_I)
    if check != GYRO.WHO_AM_I_OK:
        print "GYRO INIT FAIL..."
        return 1
    #success!
    return 0

def configAccel(rate="100HZ",fullScale="2G"):
    #-set i2c to point at accelerometer & Mag
    i2c = FT232H.I2CDevice(ft232h,XM.ADDRESS)
    #set accel data to 100Hz and enable all axes
    i2c.write8(XM.CTRL_REG1_XM, 0b01101111) #0x67) #0b01010111)  #0b01100111)
    #set the accel scale
    scaleVal = XM.RANGE_A[fullScale]
    regVal = 0x00 | scaleVal
    i2c.write8(XM.CTRL_REG2_XM,regVal)
    return 0

def configMag(rate="100HZ",fullScale="4GAUSS"):
    #-set i2c to point at accelerometer & Mag
    i2c = FT232H.I2CDevice(ft232h,XM.ADDRESS)
    #set magnetometer to 100Hz
    i2c.write8(XM.CTRL_REG5_XM,0xF4)
    #Range 4GAUSS
    i2c.write8(XM.CTRL_REG6_XM,0x20)
    i2c.write8(XM.CTRL_REG7_XM,0x00)
    #set the accel scale
    scaleVal = XM.RANGE_M[fullScale]
    regVal = 0x00 | scaleVal
    i2c.write8(XM.CTRL_REG2_XM,0x00)
    return 0


def configGyro(fullScale = "245DPS"):
    i2c = FT232H.I2CDevice(ft232h,GYRO.ADDRESS)
    i2c.write8(GYRO.CTRL_REG1_G, 0x0F) #normal mode, all axes
    #set gyro scale
    scaleVal = GYRO.RANGE_G[fullScale]
    regVal = 0x00 | scaleVal
    i2c.write8(GYRO.CTRL_REG4_G,regVal)
    print regVal
    return 0

def getAccel(fullScale):
    i2c = FT232H.I2CDevice(ft232h,XM.ADDRESS)
    cal = XM.CAL_A[fullScale]
    ax_L = i2c.readU8(0x28)
    ax_H = i2c.readU8(0x29)
    ay_L = i2c.readU8(0x2A)
    ay_H = i2c.readU8(0x2B)
    az_L = i2c.readU8(0x2C)
    az_H = i2c.readU8(0x2D)

    ax = np.int16(ax_L | (ax_H <<8))*cal
    ay = np.int16(ay_L | (ay_H <<8))*cal
    az = np.int16(az_L | (az_H <<8))*cal
    return (ax,ay,az)

def getMag(fullScale):
    i2c = FT232H.I2CDevice(ft232h,XM.ADDRESS)
    cal = XM.CAL_M[fullScale]
    mx_L = i2c.readU8(0x08)
    mx_H = i2c.readU8(0x09)
    my_L = i2c.readU8(0x0A)
    my_H = i2c.readU8(0x0B)
    mz_L = i2c.readU8(0x0C)
    mz_H = i2c.readU8(0x0D)

    mx = np.int16(mx_L | (mx_H <<8))*cal
    my = np.int16(my_L | (my_H <<8))*cal
    mz = np.int16(mz_L | (mz_H <<8))*cal
    return (mx,my,mz)

def getRotation(fullScale):
    i2c = FT232H.I2CDevice(ft232h,GYRO.ADDRESS)
    cal = GYRO.CAL_G[fullScale]
    gx_L = i2c.readU8(0x28)
    gx_H = i2c.readU8(0x29)
    gy_L = i2c.readU8(0x2A)
    gy_H = i2c.readU8(0x2B)
    gz_L = i2c.readU8(0x2C)
    gz_H = i2c.readU8(0x2D)

    gx = np.int16(gx_L | (gx_H <<8))*cal
    gy = np.int16(gy_L | (gy_H <<8))*cal
    gz = np.int16(gz_L | (gz_H <<8))*cal
    return (gx,gy,gz)

def returnIMU(queue,fullScaleA = "2G",fullScaleG = "245DPS",fullScaleM = "4GAUSS"):
    '''
    Thread dont stop, use a flag to raise exception to kill thread. Python cant stop wont stop
    '''
    while True:
        #-read accel
        ax,ay,az = getAccel(fullScaleA)
        #-read gyro
        gx,gy,gz = getRotation(fullScaleG)
        #-read mag
        mx,my,mz = getMag(fullScaleM)
        #-save values
        queue.put(np.asarray([ax,ay,az,gx,gy,gz,mx,my,mz]))


def recordSerial(queue,portGPS = 'COM3'):
    '''
    too lagit 2 quit, use a flag to raise exception to kill thread. Python cant stop wont stop
    '''
    ser = serial.Serial(port='COM3',baudrate=9600)
    while True:
        #-read messages from the serial terminal
        m = ''
        ser.flushInput()
        ser.flushOutput()
        time.sleep(1)
        while ser.inWaiting() > 0:
            #-print messages line by line
            ck = ser.readline()
            if 'GPGGA' in ck:
                m += ck
            elif 'GPVTG' in ck:
                m += ck
            #-wait for new buffer
        #flush input and output buffers
        queue.put(m)
        ser.flushInput()
        ser.flushOutput()
