AstroLib Project
================

Small library of useful astronomy routines
------------------------------------------

#Introduction
Useful set of tools that I find I need quite regularly in some of my astronomy analysis scripts.
---
####rec2pol
Function to return the 2D polar coordinates from a pair of cartesian input arguments.
```Description:
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
        theta -- angle of point relative to the x-axis.```
---
####sphdist
Function to return the spherical distance between two sets of longitude and latitiude points.
```    Description:
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
        or in radians.```
---
####timeCon
Function to convert an input value from one time type to another.
```    Description:
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
```
---
####HMS2dec
Function to convert a delimited sexagesimal RA frm HH:MM:SS to decimal degrees.
```    Input:
        ra -- string in form HH:MM:SS
        delimiter -- string indicating the delimiter e.g. :
            
    Output:
        returns the RA in decimal degrees    
```
---
####DMS2dec
Function to convert a delimited sexagesimal Dec frm DD:MM:SS to decimal degrees.
```    Input:
        dec -- string in form DD:MM:SS
        delimiter -- string indicating the delimiter e.g. :
            
    Output:
        returns the Dec in decimal degrees ```
---
####source
Class that holds a basic astronomical source data structure
####dump
Function to dump out all of the attributes of the source to screen

---
####SIMBADsrc
Source class based around the data that can be pulled from a SIMBAD query
```    Initialise:
        SIMBADsrc(source)
    Input:
        source -- string containing the source name you want to query from SIMBAD.```
#####querySIMBAD
Function that searches the SIMBAD service for basic info regarding a given source name.
```        Inputs:
            Operates on the SIMBADsrc class instance which require a source name to be initialised
        
        Outputs:
            Updates the internal attributes of the SIMBADsrc class instance with the SIMBAD data
```
---
####SIMBADsearch
Function to perform a coordinate search around ra and dec with given radius using SIMBAD
```    Input:
        ra -- region central right ascension.
        dec -- ragion central declination.
        rad -- radius of region to search in arcminutes: Default = 5
        
    Output:
        Returns a list of srcs in the region from nearest to central coords to furthest
        Each row contains, source name, dist in arcsecs from centre, object type, ra, dec.```
