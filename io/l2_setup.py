"""
l2_setup.py

DESCRIPTION
-----------
This script allows the user to easily create variable mapping of
data from a University of Wyoming King Air Research Aircraft
flight data file to a Level 2 product output file.

The accompanying program 'make_uwka_l2_flight.py' expects to 
find this script in the 'setupfilepath' variable.
If 'setupfilepath' is not provided,a simple default mapping
is performed and this script is not called.

NOTES
----- 
Python functions are used to define a dictionary for data
transfer. For further information on Python Dictionaries visit:
https://docs.python.org/2/tutorial/datastructures.html#dictionaries

For other projects (past or future):
Create a copy of this file.
Modify the dictionaries below to process.
Point 'make_uwka_l2_flight.py' to the correct directory/file.
The file may have any name, but the 'run' function must be present!

This procedure was created following the 2015 PECAN project. 
A setup script is included below for PECAN.

RETURNED
--------
variables : dictionary
    Create a variable.
    
    Structure: {InputVarName: OutputVarName}
    Example: {'GALT', 'altitude'}
    
varattr : dictionary
    Add and attribute to new variable using a list.
    The first position of the list is the attribute title, while
    the second position is the attribute information.
    
    Structure: { OutputVarName, [AttributeName, AttributeValue]}
    Example: {'longitude': ['instrument', 'Applanix AV-410']}
    
    This example will add instrument = "Applanix AV-410" to the
    'longitude' variable.
    
globattr : dictionary
    Add a global attribute.
    
    Structure: {GlobalAttributeName: GlobalAttribute}
    Example: {'Revision': 'Quality-Controlled'}
"""     
PROJECTS = ['pecan',]

#Infile = set(project.lower()).intersection(PROJECTS)
#if len(Infile)>0:
#    proj = Infile.pop()

def run():
    '''
    Convert King Air data file to a simplified default
    Level 2 output.
    '''
    variables = {
                'time': 'time',
                # Aircraft State
                'LONC': 'longitude',
                'LATC': 'latitude',
                'ztrue': 'altitude',
                'PALT': 'alt_pressure',
                'tas': 'tas',
                'aias': 'ias',
                'AVewvel': 'grnd_spd_u',
                'AVnsvel': 'grnd_spd_v',
                'AVzvel': 'grnd_spd_z',
                'AVthead': 'heading',
                'AVpitch': 'pitch_angle',
                'AVroll': 'roll_angle',
                'alpha': 'attack_angle',
                'beta': 'sideslip_angle',
                #'AVxdist': 'ew_position',
                #'AVydist': 'ns_position',
                
                # Atmospheric State
                'pmb': 'pressure',
                'ps_hads_a': 'pressure2',
                'trf': 'temperature',
                'tdplicor': 'dewpoint',
                'thetad': 'thetad',
                'thetae': 'thetae',
                'rh': 'rh',
                'mr': 'mixing_ratio',
                'lwc100': 'lwc',
                'irtc': 'irt',
                'irbc': 'irb',
                'swt': 'swt',
                'swb': 'swb',
                
                # Wind derivations
                'AVuwind': 'u_wind',
                'AVvwind': 'v_wind',
                'AVwwind': 'w_wind',
                'AVux': 'long_wind',
                'AVvy': 'lat_wind',
                'AVwdir': 'wind_dir',
                'AVwmag': 'wind_spd',
                
                # Licor Concentrations
                'co21s': 'co2_conc',
                'h2o1s': 'h2o_conc',
                
                # PCASP
                'AS200_OBR': 'pcasp_num',
                'CS200_OBR': 'pcasp_conc',
                'DBARP_OBR': 'pcasp_mean_diam',
                'PSFCP_OBR': 'pcasp_surf_area_conc',
                'PVOLP_OBR': 'pcasp_vol_conc',
                
                # CPC
                'conc_cpc': 'cpc_conc',
                
                # Miscellaneous
                'boom_pcor': 'boom_pcor',
                'turb': 'turb',
                'topo': 'topo',
                }

    varattr = {
              'lon': (['standard_name', 'instrument'], 
                      ['Longitude', 'Applanix AV-410']),

              'u_wind': (['long_name', 'standard_name', 'instrument'],
                         ['Horizontal wind, E-W component', 'x_wind', 'Applanix AV-410']),

              'topo': (['instrument', 'standard_name'],
                       ['N/A', 'Topography height']),

              #'conc_pms_2d': (['ComputationMethod'], 
              #                ['All-in Method 2']),
              }

    globattr = {
                'Revision_status': 'PECAN L2 Test File',
                }

    return variables, varattr, globattr
        
        
        
