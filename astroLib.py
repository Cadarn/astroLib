"""file astroLib.py

@brief A repository for useful astronomy tools: See individual tools for more details

@author Adam Hill <A.Hill@soton.ac.uk>

$Id: astroLib.py,v 1.0 2012/16/02 07:55:00 .
"""

import numpy as np
import urllib2
import calendar
import time
import math
import pyfits
import re
from bs4 import BeautifulSoup
import requests
from pprint import pprint

def _loopList(func):
    '''Internal decorator to be used with sphdist function to allow input of single values or lists of coordinates \n'''

    def _wrapper(long1,lat1,long2,lat2):
        '''Wrapper function to be used with sphdist function to allow input of single values or lists of coordinates \n'''
        #Attempt to convert values to floats in case non-floats are input
        try:
            long1 = float(long1)
            lat1 = float(lat1)
            long2 = float(long2)
            lat2 = float(lat2)
        except:
            pass

        #Simplest case of single values for all coords
        if type(long1) == float and type(lat1) == float and type(long2) == float and type(lat2) == float:
            return func(long1,lat1,long2,lat2)

        #Pair of cases when one pair of single coords and a vector to be compared with
        if type(long1) == float and type(lat1) == float:
            if len(long2) == len(lat2):
                dist = np.array(long2)*0.0
                for idx, val in enumerate(long2):
                    dist[idx] = func(long1, lat1, long2[idx], lat2[idx])
                return dist
            elif len(long2) != len(lat2):
                print 'Longitude and latitude vectors of unequal length'
                return None

        if type(long2) == float and type(lat2) == float:
            if len(long1) == len(lat1):
                dist = np.array(long1)*0.0
                for idx, val in enumerate(long1):
                    dist[idx] = func(long1[idx], lat1[idx], long2, lat2)
                return dist
            elif len(long1) != len(lat1):
                print 'Longitude and latitude vectors of unequal length'
                return None

        #Case where two vector pairs of coords are being used
        if len(long1) == len(lat1) and len(long1) == len(long2) and len(long2) == len(lat2):
            dist = np.array(long1)*0.0
            for idx, val in enumerate(long1):
                dist[idx] = func(long1[idx],lat1[idx], long2[idx], lat2[idx])
            return dist

    _wrapper.__doc__ = _wrapper.__doc__ + '\n' + func.__doc__
    return _wrapper            
                                                    

def rec2pol(x, y, radians=False):
    '''Function to return the 2D polar coordinates from a pair of cartesian input arguments.

    Description:
        Calculate r and theta 2D polar coordinates assuming standard Euclidean geometry
        with x,y input coordinates.  Will always return a positive value of theta between
        0 and 2pi.  Can return values in degrees (default) or radians.

    Useage:
        rec2pol(x,  y, radians=False)

    Inputs:
        x -- x cartesian coordinate.
        y -- y cartesian coordiante.

        [optional]
        radians -- Boolean value:   True to return theta in radians.
                                    False to return theta in degrees.

    Outputs: (r, theta)
        r  -- the distance from the origin of the point.
        theta -- angle of point relative to the x-axis.
    '''
    
    #Checking coords to prevent errors from (0,0) coords and keeping angle positive
    if (x,y) == (0,0):
        #If at origin then set r and theta to zero for the return
        return (0,0)

    #Calculate theta angle
    theta = math.atan(y/x)

    #If theta is negative add 2pi
    if theta < 0:
        theta += 2.0*math.pi
        
    #Convert to degrees if required
    if radians == False:
        theta = math.degrees(theta)

    #Calculate r
    r = math.sqrt(x*x + y*y)

    #Return final values
    return (r,theta)
    
@_loopList
def sphdist(long1, lat1, long2, lat2, radians=False):
    '''Function to return the spherical distance between two sets of longitude and latitiude points. \n

    Description:
        Calculate the spherical distance between two points for which longitude and
        latitude are known.  Can use values in degrees (default) or radians.

    Useage:
        sphdist(long1, lat1, long2, lat2, radians=False)

    Inputs:
        long1 -- Longitude of position 1.
        lat1 -- Latitude of position 1.
        long2 -- Longitude of position 2.
        lat2 -- Latitude of position 2.
        
        [optional]
        radians -- Boolean value:   True if inputs are in radians.
                                    False if inputs are in degrees.

    Outputs: theta
        theta -- the spherical distance between the input points in degrees (default)
        or in radians.
    '''

    #Convert to radians if given in degrees
    if radians == False:
        long1 = math.radians(long1)
        lat1 = math.radians(lat1)
        long2 = math.radians(long2)
        lat2 = math.radians(lat2)

    #Convert from spherical polar coords on unit sphere to cartesian
    vec1 = np.array([math.cos(lat1)*math.cos(long1),
                     math.cos(lat1)*math.sin(long1),
                     math.sin(lat1)])

    vec2 = np.array([math.cos(lat2)*math.cos(long2),
                     math.cos(lat2)*math.sin(long2),
                     math.sin(lat2)])

    #Calc the dot product of the vectors
    vecDot = np.dot(vec1, vec2)

    #Calc the magnitude of the cross product
    vecCross = np.cross(vec1, vec2)
    vecCrossM = np.sqrt(np.dot(vecCross, vecCross))

    #Convert to polar coords
    (r,theta) = rec2pol(vecDot, vecCrossM, radians=True)
                        
    #Convert distance to degrees if needed
    if radians == False:    
        theta = math.degrees(theta)

    return theta

def timeCon(t, inTime='MJD', outTime='MET'):
    '''Function to convert an input value from one time type to another.

    Description:
        Convert times from one type to another. Options include:
            JD -- Julian Date 
            MJD -- Modified Julian Date 
            IJD -- INTEGRAL Julian Date
            MET -- Fermi Mission Elapsed Time
            UTC -- Universal Time as year-month-dayThours:minutes:sec
                   must be given as a string e.g. '2001-10-01T12:30:34'
                   
    Useage:
        timeCon(t, inTime='MJD', outTime='MET')

    Inputs:
        t -- Input time
        inTime -- Format of input time.
        outTime -- Format of desired output time.

    Outputs: (conTime)
        conTime -- the requested time.
    '''


    def _sec2day(t):
        '''Internal function to convert seconds into days.'''
        return t/86400.0

    def _day2sec(t):
        '''Internal function to convert days into seconds.'''
        return t*86400.0
        
    #Dictionary to hold the time conversion reference points
    refTime = {'MET': 0.,
               'MJD': 51910.0,
               'JD': 2451910.5,
               'IJD': 366.0,
               'UTC': _sec2day(calendar.timegm((2001,01,01,0,0,0)))}

    #Check that the input time are catered for
    if inTime not in refTime or outTime not in refTime:
        print 'WARNING: One or both of the input time types is incorrectly entered.'
        return None

    #If input is MET need to convert from seconds to days as all other ref points
    #are in days.
    if inTime == 'MET':
        t  = _sec2day(t)

    #If input is UTC format need to convert it into seconds then to days since epoch
    if inTime == 'UTC':
        tempT = calendar.timegm(time.strptime(t.replace("-",""), "%Y%m%dT%H:%M:%S"))
        t = _sec2day(tempT)

    #Can now convert between values
    conTime = (t - refTime[inTime]) + refTime[outTime]

    #If output is MET need to convert days to seconds before returning
    if outTime == 'MET':
        conTime = _day2sec(conTime)

    #If output is UTC need to convert to days and then into UTC format
    if outTime == 'UTC':
        tempT = _day2sec(conTime)
        conTime = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(tempT))

    return conTime

def fixTSFITS(file):
    '''This function fixes the CTYPE1 header keyword in the pointLike output TS maps.
        
    Inputs:
        file -- path to desired TS map fits input file.

    Outputs: Saves back to original file with amended keyword    '''
    
    #Open file and update the CTYPE1 keyword
    hdulist = pyfits.open(file, mode='update')
    hdulist[0].header.update('CTYPE1', 'RA--ZEA')
    
    #Save the amended file
    hdulist.flush()
    hdulist.close()

    return

def HMS2dec(ra, delimiter):
    '''Function to convert a delimited sexagesimal RA frm HH:MM:SS to decimal degrees.
    
    Input:
        ra -- string in form HH:MM:SS
        delimiter -- string indicating the delimiter e.g. :
            
    Output:
        returns the RA in decimal degrees    
    '''

    #Check if delimiter is blank
    if delimiter == '':
        coord = ra.split()
    else:
        coord = ra.split(delimiter)
    #Check to see how many coords are given
    if len(coord) == 2:
        coord.append('0')
    if len(coord) == 1:
        coord.append('0')
        coord.append('0')
 
    coord[0] = float(coord[0])
    coord[1] = float(coord[1])/60.0
    coord[2] = float(coord[2])/3600.0
    
    raDecimal = np.sum(coord)*360.0/24.0
    
    return round(raDecimal,6)
    
def DMS2dec(dec, delimiter):
    '''Function to convert a delimited sexagesimal Dec frm DD:MM:SS to decimal degrees.
    
    Input:
        dec -- string in form DD:MM:SS
        delimiter -- string indicating the delimiter e.g. :
            
    Output:
        returns the Dec in decimal degrees    
    '''

    #Check if delimiter is blank
    if delimiter == '':
        coord = dec.split()
    else:
        coord = dec.split(delimiter)
        
    #Check to see if the dec is negative
    scale = 1.0
    if coord[0].strip()[0] == '-':
        scale = -1.0
  
    #Check to see how many coords are given
    if len(coord) == 2:
        coord.append('0')
    if len(coord) == 1:
        coord.append('0')
        coord.append('0')
     
    coord[0] = np.abs(float(coord[0]))
    coord[1] = float(coord[1])/60.0
    coord[2] = float(coord[2])/3600.0
    
    decDecimal = np.sum(coord)*scale
    
    return round(decDecimal,6)
  
        
class source(object):
    '''Class that holds a basic astronomical source data structure'''
        
    def __init__(self, srcID, ra, dec):
        self.id = srcID
        self.ra = ra
        self.dec = dec
        
    def dump(self):
        '''Function to dump out all of the attributes of the source to screen'''
        pprint(self.__dict__)
    

class SIMBADsrc(source):
    '''Source class based around the data that can be pulled from a SIMBAD query
    
    Initialise:
        SIMBADsrc(source)
    Input:
        source -- string containing the source name you want to query from SIMBAD.'''
    
    def __init__(self, srcID):
        self.id = srcID
        self.querySIMBAD()
        
    def querySIMBAD(self):
        '''Function that searches the SIMBAD service for basic info regarding a given source name.
    
        Inputs:
            Operates on the SIMBADsrc class instance which require a source name to be initialised
        
        Outputs:
            Updates the internal attributes of the SIMBADsrc class instance with the SIMBAD data
        '''

        #This is the general page to send the ID query to
        path = 'http://simbad.u-strasbg.fr/simbad/sim-script?submit=submit+script&script=query+id+'
        #Complete the source quer path by adding source name with spaces removed and replaced by '+' signs
        id = re.sub('\+', '%2B', self.id)
        path2 = path+'+'.join(id.split())
    
        #Use urllib2 library to quer the webpage
        request = urllib2.Request(path2)
        response = urllib2.urlopen(request)
        data = response.read()

        #Check that the source was identified
        loc = data.find('Unrecogniezd identifier')
        locb = data.find('Identifier not found in the database')
        if (loc != -1 or locb != -1):
            print 'ERROR: SIMBAD failed to recognise the source ID: %s' %self.id
            self.ra = 999
            self.dec = 999
            return
                        
    
        #Now to extract out the ra, dec
        loc1 = data.find('coord')
        loc2 = data.find('\n', loc1)
        coordData = data[loc1:loc2]
        
        #extract object type
        loc1 = data.find('main otype: ')
        loc2 = data.find('\n', loc1)
        if loc1 > -1:
            otype = data[loc1:loc2]
            self.otype = otype[12:-2]
        else:
            self.otype = 'Not found'
            
        #Check to see if we need to split on + or - depending on the dec.   
        temp = coordData.split(': ')[1].split('(')[0]    
        if temp.find('+') == -1:
            coords = temp.split('-')
            coords[1] = '-'+coords[1]
        else:
            coords = temp.split('+')
            coords[1] = '+'+coords[1]
                
        #Update internal values relating to coords
        self.ra  = HMS2dec(coords[0],'')
        self.dec = DMS2dec(coords[1],'')
                
        self.coordRef = coordData[-19:]
        
        #Hunting through the list of identifications for the variable star name if it exists
        loc1 = data.find('liste')
        loc2 = data.find('parents')
        
        #Check to see if there s a V* entry
        loc3 = data.find('V*',loc1,loc2)
        if loc3 > -1:
            loc4 = data.find('\n',loc3,loc2)
            self.vstar = data[loc3+3:loc4]
        else:
            self.vstar = 'NULL'
            
        #Hunt for the V magnitude and reference
        term = 'flux: V (Vega)'
        loc1 = data.find(term)
        if loc1 > -1:
            loc2 = data.find('\n',loc1)
            temp = data[loc1+len(term):loc2].split()
            self.Vmag = float(temp[0])
            self.VmagRef = temp[-1]
        else:
            self.Vmag = -99
            self.VmagRef = 'Not found'

        #Hunt for the Spectral type and reference
        term = 'Spectral type:'
        loc1 = data.find(term)
        if loc1 > -1:
            loc2 = data.find('\n',loc1)
            temp = data[loc1+len(term):loc2].split()
            self.specType = temp[0]
            self.specTypeRef = temp[-1]
        else:
            self.specType = 'Not listed'
            self.specTypeRef = 'Not found'        
                        
        #Store the whole SIMBAD data stream
        self.SIMBAD = data
    
def SIMBADsearch(ra,dec,rad=5.):
    '''Coordinate search around radius using SIMBAD
    
    Input:
        ra -- region central right ascension.
        dec -- ragion central declination.
        rad -- radius of region to search in arcminutes: Default = 5
        
    Output:
        Returns a list of srcs in the region from nearest to central coords to furthest
        Each row contains, source name, dist in arcsecs from centre, object type, ra, dec.'''
    #Construct the webaddress from which to search for sources    
    baseurl = 'http://simbad.u-strasbg.fr/simbad/sim-coo?Coord='
    dec = str(dec)
    if dec[0] == '+':
        dec = re.sub('+', '%2B', dec)
    elif dec[0] == '-':
        dec = re.sub('-', '%2D', dec)
    else:
        pass
    
    dataurl = baseurl+str(ra)+"+"+dec+"&Radius="+str(rad)+"&Radius.unit=arcmin"
    r = requests.get(dataurl)
    data = r.text
    r.close()

    soup = BeautifulSoup(data, "lxml")
        
    field = []
    #Search through the html and find the appropriate cols in the source table
    for tag in soup.find_all("table", {"class": "sortable"}):
        for tag2 in tag.find_all("tr"):
            row = []
            count = 0
            for td in tag2.find_all("td"):
                if (count > 0) & (count < 6):
                    row.append(re.sub('\n', '',td.get_text().strip()))
                count += 1
            if len(row) > 0:
                field.append(row)
    #Return all the sources in the field
    return field
