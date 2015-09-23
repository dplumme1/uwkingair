#! /usr/bin/env python
#******************************
# make_uwka_l2_flight.py
#******************************
"""
Create an L2 user file using L1 proccessed UWKA flight data.
Command line operation, use -h option for help.
Optional mapping dictionary can be pulled in to customize
L2 file.

Author
------
Nick Guy  6 Feb 2015.

Usage
-----
python make_uwka_l2_flight -h

TODO
----
Improve error handling.
File check for zipped files.
"""

# Load the needed packages
from netCDF4 import Dataset, stringtochar
import datetime
import numpy as np
import argparse
import os
#import sys
import imp
import traceback
#===============================================================
# Initialization defaults for variables
VERSION = '0.0.2'

# Attribute dictionaries
GLOB_ATTS = ['Source', 'Address', 'Phone', \
             'Project', 'ProjectLocation', 'ProjectPI', 'ProjectNumber',\
             'Aircraft', 'FlightNumber', 'Platform', \
             'time_coverage_start', 'time_coverage_end', \
             'FlightDate', 'TimeInterval', \
             'Conventions', 'ConventionsURL', 'ConventionsVersion',\
#             'geospatial_lat_min', 'geospatial_lat_max',\
#             'geospatial_lon_min', 'geospatial_lon_max',\
             'CenterCoordLon0','CenterCoordLat0','CenterCoordName']
             
STD_ATTS = ['long_name', 'standard_name', 'units', \
            'instrument']
            
CLD_PHYS_ATTS = STD_ATTS + ['CellSizes', 'CellSizesComment', \
                'CellSizeUnits', 'SerialNumber',\
                'ComputationMethod']

CLD_PHYS_PROBES = ['CIP', '2D-P']
#===============================================================
# BEGIN FUNCTIONS
#===============================================================

class FlightDictMap(object):
    """Class method to convert UWKA flight data to L2 user file."""
    def __init__(self, setupfilepath=None):
        '''
        Initialize class object.

        Parameters [optional]
        ----------
        setupfilepath : str
            The longpath name of the python script that contains
            a dictionary mapping of variable names to be used.
            If supplied the setup file mapping is used. 
            If not supplied then a default is provided below.
            
        '''
        self.setupfilepath = setupfilepath

        # Create dictionaries to use for mapping
        self.variables = {}
        self.varattr = {}
        self.globattr = {}

        if setupfilepath is None:
            # Use the simple default here
            self.variables, self.varattr, self.globatt = self.default_setup()
        else:
            # Use the specified setup dictionary file
            #sys.path.append(setupfilepath)
            #import l2_setup
            RunSetup = imp.load_source('run', setupfilepath)
            self.variables, self.varattr, self.globatt = RunSetup.run()

    def default_setup(self):
        variables = {
                    'time': 'time',
                    'LONC': 'longitude',
                    'LATC': 'latitude',
                    'GALT': 'altitude',
                    'PALT': 'alt_pressure',
                    'tas': 'tas',
                    'AVewvel': 'grnd_spd_u',
                    'AVnsvel': 'grnd_spd_v',
                    'AVzvel': 'grnd_spd_z',
                    'AVthead': 'heading',
                    'AVpitch': 'pitch_angle',
                    'AVroll': 'roll_angle',
                    'AVuwind': 'u_wind',
                    'AVvwind': 'v_wind',
                    'AVwwind': 'w_wind',
                    'wind_dir': 'wind_dir',
                    'wind_spd': 'wind_spd',
                    'topo': 'topo',
                    'C2DPsz2_OBL': 'conc_pms_2d',
                    }

        varattr = {
                  'longitude': (['standard_name', 'instrument'], 
                          ['Longitude', 'Applanix AV-410']),

                  'u_wind': (['long_name', 'standard_name', 'instrument'],
                             ['Horizontal wind, E-W component', 'x_wind', 'Applanix AV-410']),

                  'topo': (['instrument', 'standard_name'],
                           ['N/A', 'Topography height']),

                  'conc_pms_2d': (['ComputationMethod'], 
                                  ['All-in Method 2']),
                  }

        globattr = {
                    'Revision_status': 'Level2 transferred',
                   }
        return  variables, varattr, globattr


class L2_Process(object):
    """Class to hold the L2 file processing."""
    def __init__(self, MapDict, filePath=None, use_gattr=False):
        '''
        Initialize the class.

        Parameters [optional]
        ----------
        MapDict: FlightDictMap class
            The mapping dictionary to be used.
        filepath : str
            Long pathname to the input NetCDF file to process.
        use_gattr : bool
            True changes the MapDict mapping to suggested variables
            contained in the global attributes.
            False does nothing.    
        '''
        # Set some parameters
        self.SetupMap = MapDict
        self.infilepath = filePath

        self.outfilepath = os.path.join(os.path.dirname(filePath), 
                                       'L2_' + os.path.basename(filePath))

        # Read in the variable mapping
        #self.SetupMap = self.read_setup_mapping_file()

        # Read in the processed Flight Data NetCDF file
        self.read_uwka_flight()

        # Check whether to use global attributes for some variable mapping
        if use_gattr:
            print("Trying global attributes for mapping...")
            coords = self.ncFile.coordinates.split()
            winds = self.ncFile.wind_field.split()
            print(coords + winds)
            self._replace_dictionary_keys(coords + winds)

#            # Try to pull out coordinate information from global attributes
#            try:
#                self._test_coordinates()
#            except:
#                import sys
#                e = sys.exc_info()[0]
#                print("Error: %s" % e)
#                print("Coordinate mapping not in global attributes")
                    
            # Try pulling out wind information from global attributes
#            try:
#                self._test_winds()
#            except:
#                print("Wind mapping not in global attributes")

        # Create the output L2 User Flight NetCDF file
        self.print_msg("\nCreating %s"%self.outfilepath)
        self.outfile = Dataset(self.outfilepath, 'w', format='NETCDF4')

        # Write the dimension information to file
        self.write_dimensions()

        # Write the global attributes
        self.write_glob_attr()

        # Add variables
        self.print_msg("\nAdding variables to L2 file")
        for var in self.SetupMap.variables.keys():
            if var.lower() == 'time':
                self.write_time_var(var, self.SetupMap.variables[var])
            else:
                # Check to see if variable is found
                try:
                    #self.ncFile.variables[var]
                    self.write_var(var, self.SetupMap.variables[var])
                except:
                    self.print_msg("%s NOT FOUND IN INPUT FILE - IGNORING"%var)
                    traceback.print_exc()
                    #del self.SetupMap.variables[var]

        # Add variable attributes
        self.print_msg("\nAdding variable attributes to L2 file")
        for var in self.SetupMap.varattr.keys():
            # Check that lists have same length
            if len(self.SetupMap.varattr[var][0]) != len(self.SetupMap.varattr[var][1]):
                self.print_msg("Check %s attributes. Lists must have same length!"%var)
                return
            else:
                for index, attr in enumerate(self.SetupMap.varattr[var][0]):
                    # Check to see if variable exists
                    try:
                        self.add_attr(self.outfile.variables[var], 
                                      attr, self.SetupMap.varattr[var][1][index])
                    except:
                        pass

        # Add global attributes
        self.print_msg("\nAdding global attributes to L2 file")
        for gatt in self.SetupMap.globattr.keys():
            self.add_glob_attr(gatt, self.SetupMap.globattr[gatt])

        # Close the input file
        self.ncFile.close()

        # Close the output file
        self.outfile.close()

    def _test_coordinates(self):
        coords = self.ncFile.coordinates.split()
        self._replace_dictionary_keys(coords)

    def _test_winds(self):
        winds = self.ncFile.wind_field.split()
        self._replace_dictionary_keys(winds)

    def _replace_dictionary_keys(self, str_list):
        for key, val in self.SetupMap.variables.items():
            if val == 'longitude':
                self.SetupMap.variables[str_list[0]] = self.SetupMap.variables.pop(key)
                print("substituting %s for %s variable"%(str_list[0], key))
            elif val == 'latitude':
                self.SetupMap.variables[str_list[1]] = self.SetupMap.variables.pop(key)
                print("substituting %s for %s variable"%(str_list[1], key))
            elif val == 'altitude':
                self.SetupMap.variables[str_list[2]] = self.SetupMap.variables.pop(key)
                print("substituting %s for %s variable"%(str_list[2], key))
            elif val == 'time':
                self.SetupMap.variables[str_list[3]] = self.SetupMap.variables.pop(key)
                print("substituting %s for %s variable"%(str_list[3], key))
            elif val == 'wind_spd':
                self.SetupMap.variables[str_list[4]] = self.SetupMap.variables.pop(key)
                print("substituting %s for %s variable"%(str_list[4], key))
            elif val == 'wind_dir':
                self.SetupMap.variables[str_list[5]] = self.SetupMap.variables.pop(key)
                print("substituting %s for %s variable"%(str_list[5], key))
            elif val == 'w_wind':
                self.SetupMap.variables[str_list[6]] = self.SetupMap.variables.pop(key)
                print("substituting %s for %s variable"%(str_list[6], key))

    def print_msg(self, str):
        print(str)

    #########################
    ##  READ FILE METHODS  ##
    #########################

    def read_setup_mapping_file(self):
        '''Read in the setup variable mapping file'''
        map = np.genfromtxt(self.setupfilepath, dtype=None, \
                               names=['col1', 'col2', 'col3'], \
                               delimiter="	", skip_header=15) 
        return map

    def read_uwka_flight(self):
        '''Read in processed NetCDF file.'''

        # Read the NetCDF
        self.ncFile = Dataset(self.infilepath,'r')

        # Read in Dimensions
        self.nc_dims = self.ncFile.dimensions
        self.nc_dimnames = [dim for dim in self.ncFile.dimensions]

        # Read in Variable names
        self.ncvars = self.ncFile.variables

        # Grab the metadata stored in global attributes as a dictionary
        self.metadata = self.ncFile.__dict__  # Gets entire attribute info
        self.globattr = self.ncFile.ncattrs() # Gets only attribute names

    ############################
    ##  WRITE GLOBAL METHODS  ##
    ############################

    def write_dimensions(self):
        # Write the dimension information
        for dimname in self.ncFile.dimensions:
            self.outfile.createDimension(dimname, len(self.ncFile.dimensions[dimname]))

    def write_glob_attr(self):
        '''Write the global attributes.'''
        self.global_atts = {
                            'description': "UWKA Level 2 Data",
                            'documentation': "http://flights.uwyo.edu/n2uw/users/",
                            'created_UTC': datetime.datetime.utcnow().isoformat(),
                            'coordinates': "longitude latitude altitude time",
                            'latitude_coordinate': "latitude",
                            'longitude_coordinate': "longitude",
                            'zaxis_coordinate': "altitude",
                            'time_coordinate': "time",
                            'wind_field': "wind_dir wind_spd w_wind",
                            }
        for attname in GLOB_ATTS:
            self.global_atts[attname] = str(getattr(self.ncFile, attname))
        self.outfile.setncatts(self.global_atts)

    def add_glob_attr(self, attname, attval):
        '''Add a global attribute from setup file.'''
        self.outfile.setncattr(attname, attval)

    ##########################
    ##  WRITE FILE METHODS  ##
    ##########################       
    def write_time_var(self, tIn, tOut):
        '''Create a time variable for user file'''
        self.outfile.createVariable(tOut, \
                                    self.ncvars[tIn].datatype, \
                                    self.ncvars[tIn].dimensions)
        self.outfile.variables[tOut].long_name = "time of measurement"
        self.outfile.variables[tOut].standard_name = "time"
        self.outfile.variables[tOut].units = self.ncvars['time'].units
        self.outfile.variables[tOut].strptime_format = "seconds since %F %T %z"
        self.outfile.variables[tOut][:] = self.ncvars[tIn][:]

    def write_var(self, VarIn, VarOut):
        '''Create output variable'''
        # Get the attributes associated with a variable
        InputVaratts = self.ncFile.variables[VarIn].ncattrs()

        # Check for _FillValue
        if hasattr(self.ncvars[VarIn], '_FillValue'):
            fillval = self.ncvars[VarIn]._FillValue
        elif hasattr(self.ncvars[VarIn], 'missing_value'):
            fillval = self.ncvars[VarIn].missing_value
        else:
            fillval = None

        self.outfile.createVariable(VarOut, \
                                       self.ncvars[VarIn].datatype, \
                                       self.ncvars[VarIn].dimensions,\
                                       fill_value=fillval)

        # Copy the variable attributes
        self.outfile.variables[VarOut].setncatts({k: self.ncvars[VarIn].getncattr(k) for k in self.ncvars[VarIn].ncattrs()})
        
        # Fill the data values
        self.outfile.variables[VarOut][:] = self.ncvars[VarIn][:]

        # Add the category attribute
#        self.outfile.variables[VarOut].Category = self.CAT

        # Set which attribute list to use
#        if self.CAT == 'Cloud_physics':
#            attribute_list = CLD_PHYS_ATTS

#        else:
#            attribute_list =  STD_ATTS

#        for name in attribute_list:
#            # If the attribute is there then assign
#            if hasattr(self.ncFile.variables[VarIn], name):
#                attval = self.ncFile.variables[VarIn].__getattribute__(name)
#                self.outfile.variables[VarOut].__setattr__(name, attval)
#            else:
#                print("'%s' attribute not found - need to assign"%name)


    def add_attr(self, Var, Attname, Attvalue):
        Var.__setattr__(Attname, Attvalue)
        self.print_msg("   '%s' added"%Attname)

#######################
##  PROGRAM STARTUP  ##
#######################
if __name__ == '__main__':

    # Check for directory
    parser = argparse.ArgumentParser(\
             description='Start UWKA Flight Data User File Creation .')
    parser.add_argument('inputfile', type=str, help=' Full path of input file to process')

    igroup = parser.add_argument_group(
    title='Set input directory, optional',
    description=(''
    'An input file may be specified. If not a default will be created. '
    'The path for the input mapping directory can be set. '
    'If not chosen, the current working directory is assumed. '
    ' '
    'The following flags are availble.'
    ' '))

    igroup.add_argument('-s', '--setupfilepath',
    help='Full path of the setup mapping file',
    default=None)
    
    parser.add_argument('-b', action='store_true',
    help='Use built-in coordinate information in global attributes')

    parser.add_argument('-v', '--version', action='version',
    version='Version %s' % (VERSION))

    # Parse the args
    args = parser.parse_args()

    # If there is a mapping file then process
    SetupMap = FlightDictMap(setupfilepath=args.setupfilepath)
    
    L2_Process(SetupMap, filePath=args.inputfile, use_gattr=args.b)