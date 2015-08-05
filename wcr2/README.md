WCR2
====

Wyoming Cloud Radar 2

In this directory find examples of IDL routines that read and process WCR2 data.
Note that these are merely examples and there are embedded functions that are not provided here.

If you would like more information on the routines or to find out about obtaining the
full library (WCR2TOOLS) please contact:

Dr. Samuel Haimov (haimov@uwyo.edu)

Below is information from Sam:

Wcr2tools library provides idl routines to work with the WCR (WCR2) raw data 
files, UW KingAir data files, NCAR C130 aircraft, and WCR data archived in 
netCDF files after August 2009. Before August 2009 an obsolete WCRTOOLS
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

##Routines for handling WCR raw data

closewcr2.pro         - closes open WCR2 raw file

initwcr2.pro          - returns a WCR2 parameters structure for a 
                        WCR2 raw data file
                        
loadwcr2conf.pro      - loads WCR2 configuration structure

openwcr2.pro          - opens WCR2 raw file, reads file header(s), and 
                        returns a structure withfile information
                        
plotwcr2crcal.pro     - processing and plotting Corner Reflector 
                        calibration
                        
plotwcr2qam.pro       - make ps, pdf and txt files with data quality 
                        plots and radar parameters statistics for 
                        multiple WCR2 raw files
                        
plotwcr2qa.pro        - makes plots of radar noise and other data 
                        quality checks
                        
printwcr2info.pro     - print a WCR2 raw file parameters and statistics

readwcr2data.pro      - returns a WCR2 data structure of a WCR2 data file

readwcr2raw.pro       - reads WCR2 raw file records (radar profiles)

reopenwcr2.pro        - re-opens WCR2 raw file for which wcr2finfo 
                        structure exists
                        
wcr2beamcal.pro       - analyze WCR2 antenna beam angles in aircraft 
                        reference frame and returns the optimized beam 
                        unit vector for that antenna
                        
wcr2beamfind.pro      - returns beam unit vector (pointing angles in 
                        aircraft(AC) coordinate system) of a WCR2 antenna
                        beam optimized using ground returns from SPP or 
                        CPP WCR2 raw data files
                        
wcr2calcoeff.pro      - calculate hh-pol calibration coeff. from corner 
                        reflector measurements
                        
wcr2conf__define.pro  - defines WCR2 configuration named structure

wcr2conf.pro          - returns experiment specific WCR2 configuration 
                        structure with information about the radar 
                        installation and calibration parameters
                        
wcr2crdata.pro        - returns WCR2 Corner Reflector(CR) measurement data 
                        structure
                        
wcr2fh__define.pro    - defines WCR2 raw data first, fixed size file 
                        header named structure (for all radar algorithms)
                        
wcr2crproc.pro        - extracts WCR2 corner reflector power calibration 
                        data
                        
wcr2pwrcalco.pro      - calibrate WCR2 received co-polarized power in 
                        equivalent Z [mm^6/m^3]. Mean (moving average 
                        with Nswin window) noise power is subtracted (no 
                        thresholding) and range correction is applied. 
                        Power values below the mean noise (negative) are
                        left in the output (to allow proper power 
                        averaging later on if desired).
                        
wcr2pwrcalcr.pro      - calibrate WCR2 received cross-polarized power in
                        equivalent Z [mm^6/m^3].
                        
wcr2pwrprod.pro       - returns reflectivity structure for a given WCR2 
                        raw file
                        
wcr2rh__define.pro    - defines WCR2 raw data record (radar profile)
                        header named structure
                        
wcr2rxnfg.pro         - calculates Noise Figure and receiver gain for WCR
                        from hot load and warm load power measurements. 
                        Sufficient smoothing (SM keyword) must be used to
                        bring the uncertainty in the power measurements
                        well below the expected noise power increase for
                        the hot load (around 0.2 dB)
                        
wcr2timeoff.pro       - calculates time offset between WCR2 time stamps 
                        and platform data system time stamps.  It also 
                        returns the corresponding correlation coefficient.
                        
writewcr2raw.pro      - re-writes(copy) a fraction of a  WCR2 raw file to a 
                        new file.
                        

##Routines for creating WCR NetCDF files from WCR raw data

wcr2cdl2nc.       - Generate WCR netcdf file from a prototype cdl 
                    file (ncgen utility)
                    
wcr2createnc      - generate netcdf file(s) from a given cdl file, 
                    WCR2 raw data file(s), and optional WCR2 platform
                    motion(aircraft) navigational data file
                    
wcr2threshnc      - Return range adjusted reflectivity threshold,
                    fltarr(rangegates,profile), given radar range gates in 
                    meters and reflectivity noise standard deviation in 
                    equivalent dBZ@1km
                    
wcr2write1nc      - calculate radar received power related products
                    from a WCR2 raw file and write them to an already
                    generated (by wcr2cdl2nc) netcdf file.
                    Other ancillary variables are also recorded
                    
wcr2write2nc      - calculate radar received Doppler velocity related
                    products from a WCR2 raw file and write them to 
                    an already generated (by wcr2cdl2nc) netcdf file 
                    Other ancillary variables are also recorded,

##Routines associated with WCR Doppler correction for AC motion

ac2wcr     - Interpolates an aircraft variable time dimension 
             to match WCR time dimension (profiles) or WCR variable
             profile dimension to match AC time (WCRTOOLS library)
             
ac2gt      - Returns single precision AC (aircraft) to EARTH or 
             EARTH to AC coordinate transformation matrix (WCRTOOLS library).
             
wcrdopcor  - Returns WCR radial Doppler velocity corrected for AC motion.
             Positive velocity is toward the radar (WCRTOOLS library).
             
wcrwunfold - Returns shifted/unfolded WCR Doppler corrected for AC motion
             velocity around given wind value component into the WCR beam
             (WCRTOOLS library).

acdata     - returns UW_KingAir/NCAR_C130 selected INS and 
             other AC data variables from a netCDF file. 
             There is only one time dimension in the loaded
             variables. See the notes.  Time variable is 
             resampled to match data sampling rate (sps).
             
acwcr2data - returns an aircraft(AC) data structure associated
             with a WCR2 configuration
             
wcr2dopcor - returns WCR2 radial Doppler velocity w.r.t. 
             ground corrected for aircraft (AC) motion. 
             Positive velocity is toward the radar. 
             
wcr2dopprod- returns corrected Doppler velocity (and misc.
             velocity related variables) structure for given
             WCR2 and corsponding aircraft data files
             If aircraft data is not available provides the 
             uncorrected Doppler. 

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
                
wcrud2udnc  -   Returns the combined reflectivity or velocity field (nrg, nrpofs).  
                from the 'up' and the 'down' WCR beams. The ranges added in the
                gap between the first valid 'up' and 'down' gates are filled with 
                missing values.  The missing values of the output are always NaN.
                (WCRTOOLS library).
                
wcrwdvnc    -   Adjust the airborne radar radial Doppler velocity from the 
                Up/Down-pointing antenna for horizontal wind contribution.
                The routine uses Doppler velocity already corrected for the 
                aircraft motion contribution.
                (WCRTOOLS library).
                
h2omwir     -   Returns complex refraction index of pure water/ice for a given
                microwave frequency in GHz and temperature in degree Celsius
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