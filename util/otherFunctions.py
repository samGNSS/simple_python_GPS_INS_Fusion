from __future__ import division #enables default float division
import numpy as np #used for signal processing
import serial #used to communicate with the sensors
import time #used for timing sections of code, and will be used to meter sensor communications

#Function record_GPS:
# -records NMEA over a serial connection for a specified time
# -inputs:
# --recording time, in seconds
# --baud rate, in symbols per second
# --update rate, in Hz
# --serial port, example: /dev/ttyUSB0
# --output file name
# --fileIO, used to overwrite or append to a file
# -outputs:
# --1 for success
# --(-1) for serial port error
# --(-2) for file i/o error
def record_GPS(recordTime = 1, baudRate = 9600, updateRate = 1, serialPort = '/dev/ttyUSB0',\
               fileName = 'nmeaRecord1sec.txt',fileIO = 'w'):
    #set up serial port
    try:
        ser = serial.Serial(
        port= serialPort,
        baudrate= baudRate,
    )
    except:
        print "Error setting up serial port"
        return -1

    #-open file
    try:
        fid = open(fileName,'w')
    except:
        print "Error opening file"
        return -2

    #start recording
    #-define empty string to hold incoming data
    serialReceived = ''
    #-flush serial buffers, and wait 1/updateRate seconds before continuing
    #ser.flushInput()
    #ser.flushOutput()
    #time.sleep(1/updateRate)
    #-get starting time
    startTime = time.time()
    #-loop for the required time
    while time.time() - startTime < recordTime:
        #-while the serial buffer is transmitting: get data
#        timeTestS = time.time()
        while ser.inWaiting() > 0:
            #-print messages line by line
            serialReceived += ser.readline()
#            print serialReceived
#            print '\n'
#        print time.time() - timeTestS
        #write data to the file
        fid.write(serialReceived)
        #-clear serialReceived
        serialReceived = ''
        #-flush serial buffers, and wait 1/updateRate seconds before continuing
        ser.flushInput()
        ser.flushOutput()
        #time.sleep(1/updateRate)
    print "record over"
    fid.close()
    return 1
#end function: recordGPS

#Function parseNMEA:
#function parses a file containing NMEA data. Currently only parses: GGA,VTG, and GSV
#inputs:
#-fileName, file to parse
#-numpy flag, if one returns numpy arrays
#outputs:
#-stuff...
def parseNMEA(fileName = 'nmeaRecord4.txt',numpyFlag = 1):
    #read entire file, and store the contents in an list
    with open(fileName,'r') as fid:
        fileContents = fid.readlines()

    #define returns as empty lists
    #-latitude
    lat = []

    #-longitude
    lon = []

    #-altitude
    alti = []

    #-velocity
    vel = []

    #-satellite PRN
    prn = []

    #-satellite elevation angle
    ele = []

    #-satellite azimuth angle
    azi = []

    #-satellite signal to noise ratio
    snr = []

    #loop over fileContentents and fill the lists
    for line in xrange(len(fileContents)):
        #-Look for GGA message, if found extract latitude, longitude, and altitude
        if 'GPGGA' in fileContents[line]:
            parse = fileContents[line].strip().split(',')
            try:
                lat.append(float(parse[2]))
                lon.append(float(parse[4]))
                alti.append(float(parse[9]))
            except:
                print "GGA err"
        #-Look for VTG message, if found extract velocity
        elif 'GPVTG' in fileContents[line]:
            parse = fileContents[line].strip().split(',')
            try:
                vel.append(float(parse[7]))
            except:
                print "VTG err"
        #-Look for GSV, if found extract satellite PRN, elevation angle, and azimuth angle
        elif 'GPGSV' in fileContents[line]:
            parse = fileContents[line].strip().split(',')
            indx = [i for i, j in enumerate(parse) if j == '']
            for i,j in enumerate(indx):
                parse[j] = '0'
            numSatsInMes = (len(parse) - 4)//4
            try:
                for i in xrange(numSatsInMes):
                    prn.append(int(parse[4 + 4*i]))
                    ele.append(int(parse[5 + 4*i]))
                    azi.append(int(parse[6 + 4*i]))
                    if i != numSatsInMes-1:
                        snr.append(int(parse[7 + 4*i]))
                    else:
                        snr.append(int(parse[7 + 4*i][0:2]))
            except:
                if i == numSatsInMes-1:
                    snr.append(0)
                else:
                    print "GSV err"

    if numpyFlag == 1:
        return ((np.asarray(lat), np.asarray(lon), np.asarray(alti)), np.asarray(vel) , \
                 (np.asarray(prn), np.asarray(ele), np.asarray(azi), np.asarray(snr)))
    else:
        return ((lat, lon, alti), vel,(prn, ele, azi, snr))
#end function: parseNMEA

#Function: convert2decimalDeg
#function to convert from degrees minutes and decimal minutes, to decimal degrees
#-inputs
#-dd mm.mmmm
#-outputs
#-dd.dddddd
def convert2decimalDeg(degMin):
    mins = degMin % 100
    deg = np.floor(degMin/100)
    return deg + mins/60
#end function: convert2decimalDeg

#Function: convert2DegMM
#function to convert from decimal degrees to degrees and decimal minutes
#-inputs
#-dd.dddddd
#-outputs
#-dd mm.mmmm
def convert2DegMM(decDeg):
    deg = np.floor(decDeg)
    mins = (decDeg - deg)*60
    return deg*100 + mins
#end function: convert2decimalDeg

#Function: geodetic2ecef
#function converts data from geodetic to earth centered earth fixed coordinates
#inputs:
#-latitude
#-longitude
#-altitude
#-spheroid model
#-convert flag, convert lat and lon to decimal degrees
#outputs
#-X
#-Y
#-Z
def geodetic2ecef(lat,lon,alti,model = 'wgs84',convert = 1):
    #-used wgs84 as the ellipsoid model
    #-set squared eccentricity
    eccenSquare = 6.69437999014e-3
    #-set semi-major axis
    semiMajorAxis = 6378137 #meters

    #-check convert flag
    if convert == 1:
        lat = np.pi*convert2decimalDeg(lat)/180
        lon = np.pi*convert2decimalDeg(lon)/180

    #compute distance to the Z-axis
    latitudeNormal = semiMajorAxis/np.sqrt(1 - eccenSquare*(np.sin(lat)**2))

    #-compute x coordinate
    xCoor = (latitudeNormal + alti)*np.cos(lat)*np.cos(lon)
    yCoor = (latitudeNormal + alti)*np.cos(lat)*np.sin(lon)
    zCoor = (latitudeNormal*(1-eccenSquare) + alti)*np.sin(lat)

    return (xCoor,yCoor,zCoor)
#end function: geodetic2ecef

#Function: ecef2enu
#function converts data from ecef to enu
#inputs:
#-x current
#-y current
#-z current
#-x reference
#-y reference
#-z reference
#-latitude reference
#-longitude reference
def ecef2enu(xc,yc,zc,xr,yr,zr, latr,lonr):
    latr *= np.pi/180
    lonr *= np.pi/180
    #-subtract reference
    deltaPos = np.array([[xc-xr],[yc-yr],[zc-zr]])
    rotation = np.array([[-np.sin(lonr),np.cos(lonr),0],[-np.sin(latr)*np.cos(lonr),-np.sin(latr)*np.sin(lonr),np.cos(latr)]\
                          ,[np.cos(latr)*np.cos(lonr),np.cos(latr)*np.sin(lonr),np.sin(latr)]])
    return [np.dot(rotation,deltaPos[:,:,i]) for i in xrange(len(xc))]
#end function: ecef2enu
