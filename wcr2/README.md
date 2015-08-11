WCR2
====

Wyoming Cloud Radar 2

In this directory find examples of IDL routines that read and process WCR2 data.
Note that these are merely examples and there are embedded functions that are not provided here.

If you would like more information on the routines or to find out about obtaining the
full library (WCR2TOOLS) please contact:

Dr. Samuel Haimov (haimov@uwyo.edu)

Below is information from Sam:

WCR2TOOLS library provides idl routines to work with the WCR (WCR2) raw data 
files, UW KingAir data files, NCAR C130 aircraft, and WCR data archived in 
netCDF files from August 2009 - 2014. Before August 2009 an obsolete WCRTOOLS
library must be used for handling raw and processed data. It is associated 
with the previous(retired) WCR hardware/software (WCR1).A few of the routines
in WCRTOOLS library are either used or useful for WCR2 and are described below.

It is expected that this library and the authors of the routines used
will be acknowledged by the users if the results of its use are published.

Generally, the use of the library requires Samuel Haimov's idltools idl library 
and Ray Sterner's jhuapl idl library. You need to set the IDL path to include
the IDL directories for wcr2tools, idltools and jhuapl.


Following is a list of the available idl routines with short descriptions.
The same information is accessable within idl by executing wcrtools.pro routine
(IDL> wcrtools).  Every routine includes a keyword HELP providing additional
information.  More general information is available in the Readme file in
idltools library.
________________________________________________________________________________

##Reading variables from AC, WCR or WCR/AC NetCDF files

ac2loadnc        - returns KingAir or NCAR C130 aircraft data 
                   variable from a netCDF file. There is only one 
                   time dimension in the returned variable. 
                    
wcr2loadnc       - load variables from a WCR netcdf file

##Miscelaneous routines related to WCR and AC NetCDF data

wcr2calthreshnc - return thresholded reflectivity in dBZ
                  Optionally it can also return thresholded Doppler
                  velocity field (expnded version of wcr2threshnc)
                  
wcr2plotnc      - display reflectivity and Doppler velocity for 
                  selected WCR2  netcdf files

wcrud2altnc -   Returns WCR up(beam1), down(beam1), or down&up(beam1&beam2) 
                reflectivity(or velocity) resampled into the vertical plane and 
                adjusted for the radar altitude. Deviations of the beam(s) in the 
                horizontal plane are ignored (horizontal homogeneity is assumed). 
                If beam(s) pointing angles are not given assume no deviation from 
                vertical (straight leg, small attitude angles)
                (WCRTOOLS library).
                
wcrwdvnc    -   Adjust the airborne radar radial Doppler velocity from the 
                Up/Down-pointing antenna for horizontal wind contribution.
                The routine uses Doppler velocity already corrected for the 
                aircraft motion contribution.
                (WCRTOOLS library).
             
wcrwunfold - Returns shifted/unfolded WCR Doppler corrected for AC motion
             velocity around given wind value component into the WCR beam
             (WCRTOOLS library).
________________________________________________________________________________

##DISCLAIMER:

  The above described routines may be used, copied, or redistributed as long 
  as they are not sold and the copyright notice (included in the routine) is 
  reproduced on each copy made. If results using these routines are published
  appropriate acknowledgment is expected.

  The software routines are provided as is and without any expressed
  or implied warranties whatsoever.  All warranties including, but not
  limited to, performance, merchantability or fitness for a particular
  purpose are hereby disclaimed.  You assume the entire risk and liability
  of using the software to include use in compliance with any third party
  rights.  You are advised to test the software thoroughly before relying
  on it.  In no event shall Samuel Haimov or UWyo/DAS be held responsible 
  for any direct, indirect, consequential or inconsequential damages or  
  monetary losses arising out of the use or inability to use the software.