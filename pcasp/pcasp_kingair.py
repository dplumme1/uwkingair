#!/usr/bin/env python
#####################
# pcasp_kingair - Quick quality check of PCASP
#####################
"""
 NOTES::
 To run:
 $ pcasp_kingair.py
 
 or alternatively
 $ python pcasp_kingair.py
 
 If this program is not executable.
 
 View the output using 'display' on bat or connecting to server and opening file in local OS.
 
 Measures particles  from 0.1 - 3 micron in diameter.
 There are 31 channels, this is due to the way the instrument is processed by NCAR software.
 Basically channel 0 is a throw away.
 In the processed file, variables with 'OBR' pertain to PCASP.  
   'W' is the regular var
   'WS' is the sheath var
   'WC' is density corrected.  Though this is the best to use post-processed look at the 
   other two for check.
 
 Generally decreases with altitude.
 Statistics are poor at the larger particle sizes.  This can be seen by the error
 bars on the second page spectra plots.
 
 On the PCASP/CPC size plot, the PCASP should be lower.
 Higher altitudes should show fewer particles and vice versa.
 The beginning will be noisy b/c of dust stirred up by props.
 The hotwire LWC shows clouds, spikes form PCASP here are due to particle shattering.
 The sample and ambient flow are 10:1 ratio and therefore this plot should line up
 closely b/c of the factor of 10 division.
 Ch0 plot: Smalest particles  goes to >100 if large vibrations are apparent (e.g. takeoff)
 Ch1 is the same.
 The correlation b/w LWC and CCN is an artifact due to shattering.  Known caveat by Jeff
 and is here b/c of the way it is processed.
 
 MODIFICATIONS::
 29 May 2015 - Nick Guy, ported from original IDL code by Jeff Snyder.
               A small fix was applied to dNdlog calculation.
"""
# import required packages
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib
#import xray
import netCDF4
import pandas as pd
import numpy as np

########################
### GLOBAL VARIABLES ###
########################
VERSION = "2015.07.15"

###########################
### USER INPUT BEGIN ###
###########################
# Set up path to file
PROJECT = "pecan15"

# Change the min/max of each plot by adjusting below
LIMS = {"tas" : (0., 150.),
        "alt" : (1000., 10000.),
        "flow" : (0., 2.),
        "ch0" : (-100., 900.),
        "ch1" : (-100., 900.),
        "chsum" : (-100., 900.),
        "cpc_comp" : (0., 5000.),
        "lwc" : (0., 1.5),
}
##fDir = "/netdata/kingair_data/" + project + "/work/"
#fDir = "/kingair_data/" + project + "/work/"
#nc_filename = '20150506.c1.nc'
#nc_path_filename  = fDir + nc_filename
#Imagepath = "/home/bguy/images/pcasp"
#createImage = True # True to save, False to display on screen


# Create a list of date-time pairs to use for subsetting
# To do so add any number of dtpairs.append() after dtpairs is defined
dts = []
#dts.append([(2015, 5, 27, 22, 5, 0), (2015, 5, 27, 22, 10, 0)])
#dts.append([(2015, 5, 27, 22, 35, 0), (2015, 5, 27, 22, 40, 0)])

#dts.append([(2015, 5, 6, 16, 50, 0), (2015, 5, 6, 17, 0, 0)])
#dts.append([(2015, 5, 6, 17, 31, 0), (2015, 5, 6, 17, 33, 0)])
#dts.append([(2015, 5, 6, 17, 40, 0), (2015, 5, 6, 17, 45, 0)])
#dts.append([(2015, 5, 6, 17, 59, 0), (2015, 5, 6, 18, 3, 0)])

###########################
### USER INPUT END ###
###########################

class CheckPCASP(object):
    """
    A class mentod to check PCASP measurements from UWyo King Air
    flight file.
    """
    def __init__(self, filepathIn, dts, flow3_range=(0., 5.)):
        '''Initialize the class.
        
        Parameters::
        ----------
        filepathIn = string 
            Full path to input file.
        dts = list
            A list of (start, stop) dates to be used to generate subsets.
            The format should be (year, month, day, hour, minute, second)
            for both start and stop and each value should an integer.
            An empty list results in no subset being processed or plotted.
        flow3_range = fltarr arr (2)
            (min, max) array to used to set flow values to allow to pass.
        '''
        self.filepathIn =  filepathIn
        self.filename = os.path.basename(self.filepathIn)
        self.dts = dts
        self.flow3_range = flow3_range
        self._grab_vars()

    def _grab_vars(self):
        '''Grab the variables from the netcdf file and create Pandas time series'''
        # Open the NetCDF file
        self.nc = netCDF4.Dataset(self.filepathIn, 'r')
#data = xray.open_dataset(nc_path_filename, decode_cf=False)
#df = data.to_dataframe()

        # Grab the variables from the file
        tas = self.nc.variables['tas']
        alt = self.nc.variables['ztrue']
        flow1 = self.nc.variables['PFLW_OBR']
        flow2 = self.nc.variables['PFLWS_OBR']
        flow3 = self.nc.variables['PFLWC_OBR']
        self.pcasp = self.nc.variables['AS200_OBR']
        cpc = self.nc.variables['conc_cpc']
        lwc = self.nc.variables['lwc100']
        
        # Remove the negative PCASP values
        self.pcasp[:][self.pcasp[:] < 0.] = np.nan
        
        # Calculate the PCASP mid-bin location
        self.pcasp_bin_mid = (self.pcasp.CellSizes[1:] + self.pcasp.CellSizes[0:-1]) /2.

        # Pull in time 
        timesec = self.nc.variables['time']
        self.time = netCDF4.num2date(timesec[:], timesec.units)
        
        # Create Pandas Series for easy workflow
        self.tass = pd.Series(tas, index=self.time, name='TAS')
        self.alts = pd.Series(alt, index=self.time, name='Alt')
        self.flow1s = pd.Series(flow1, index=self.time, name='Flow 1')
        self.flow2s = pd.Series(flow2, index=self.time, name='Flow 2')
        self.flow3s = pd.Series(flow3, index=self.time)
        self.cpcs = pd.Series(cpc, index=self.time, name='CPC')
        self.lwcs = pd.Series(lwc, index=self.time, name='LWC')
        # Calculate summations of PCASP channels, channel 0 is a bunk
        self.ch0s = pd.Series(self.pcasp[:,0,1], index=self.time, name='Channel 0')
        self.ch1s = pd.Series(self.pcasp[:,0,2], index=self.time, name='Channel 1')
        self.chsums = pd.Series(self.pcasp[:,0,3:29].sum(1), index=self.time, name='Channels 2-29')
        
        # Remove the flows outside of a determined range
        self.flow3s[(self.flow3s <= self.flow3_range[0])] = np.nan
        self.flow3s[(self.flow3s >= self.flow3_range[1])] = np.nan

        # Calculate d(log10{D}) for PCASP instrument
        self.dlog10d = np.log10(self.pcasp.CellSizes[1:] / self.pcasp.CellSizes[0:-1])        
        
        # Create Pandas DataFrame for PCASP data
        self.pcaspf = pd.DataFrame(data=self.pcasp[:,0,1:], index=self.time, \
                                   columns=self.pcasp_bin_mid)

        # Close the NetCDF file
        self.nc.close()
    
    def _draw_time_region(self, ax):
        '''A function to highlight time regions for further analysis.
        
        Parameters::
        ----------
        ax = Matplotlib axis instance 
            Used for drawing polygons on that axis.
        '''
        # Get the y-limits to match the polycollection
        ylim = ax.get_ylim()
        # Loop through the times and create
        for fills in self.dts:
            stime = netCDF4.datetime(fills[0][0], fills[0][1], fills[0][2],\
                        fills[0][3], fills[0][4], fills[0][5])
            etime = netCDF4.datetime(fills[1][0], fills[1][1], fills[1][2],\
                        fills[1][3], fills[1][4], fills[1][5])
            ax.fill_between(self.time, ylim[0], ylim[1], \
                            where=(self.time>=stime)&(self.time<=etime), alpha=0.3)
        
    def plot_timeseries(self, project, save=False, save_dir=None, fancy=True):
        '''Create an timeseries plot of specific variables.
        
        Parameters::
        ----------
        save = boolean 
            True to save an output plot, False displays on screen.
        fancy = boolean
            If set to true uses ggplot style, otherwise default is used.
        '''
        if fancy:
            matplotlib.style.use('ggplot')
        
        # Decrease the legend font size
        matplotlib.rcParams.update({'legend.fontsize': 'x-small',
                                    'legend.frameon':False,
                                    'legend.framealpha':0.5})
        
        # Make title
        figtitle = project + " NOAA PCASP " + self.filename
        
        # Create the figure/axes instances
        fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(nrows=8, ncols=1, figsize=(8, 11))
        
        # Plot True air speed
        self.tass.plot(style='k-', ax=ax1, sharex=True, legend=True)
        ax1.set_ylabel('m/s')
        ax1.set_ylim(bottom=LIMS['tas'][0], top=LIMS['tas'][1])
        ax1.grid(which='minor')
        ax1.set_title(figtitle)
        self._draw_time_region(ax1)

        # Plot Altitude
        self.alts.plot(style='k-', ax=ax2, legend=True)
        ax2.set_ylabel('m')
        ax2.set_ylim(bottom=LIMS['alt'][0], top=LIMS['alt'][1])
        ax2.grid(which='minor')
        self._draw_time_region(ax2)

        # Plot flow rates
        self.flow1s.plot(style='k-', ax=ax3, legend=True)
        self.flow2s.divide(10.).plot(style='r:', ax=ax3, legend=True) # Divide by 10
        ax3.set_ylabel('Flows')
        ax3.set_ylim(bottom=LIMS['flow'][0], top=LIMS['flow'][1])
        ax3.grid(which='minor')
        self._draw_time_region(ax3)

        # Plot Channel # 0 counts
        self.ch0s.plot(style='k-', ax=ax4, legend=True)
        ax4.set_ylabel('Counts')
        ax4.set_ylim(bottom=LIMS['ch0'][0], top=LIMS['ch0'][1])
        ax4.grid(which='minor')
        self._draw_time_region(ax4)

        # Plot Channel # 1 counts
        self.ch1s.plot(style='k-', ax=ax5, legend=True)
        ax5.set_ylabel('Counts')
        ax5.set_ylim(bottom=LIMS['ch1'][0], top=LIMS['ch1'][1])
        ax5.grid(which='minor')
        self._draw_time_region(ax5)

        # Plot Channels 2-29 Sum
        self.chsums.plot(style='k-', ax=ax6, legend=True)
        ax6.set_ylabel('Counts')
        ax6.set_ylim(bottom=LIMS['chsum'][0], top=LIMS['chsum'][1])
        ax6.grid(which='minor')
        self._draw_time_region(ax6)

        # Plot PCASP & CPC concentrations
        self.chsums.divide(self.flow3s).plot(style='k-', ax=ax7, legend=True)
        self.cpcs.plot(style='r:', ax=ax7, legend=True)
        ax7.set_ylabel(r"cm$^{-3}$")
        ax7.set_ylim(bottom=LIMS['cpc_comp'][0], top=LIMS['cpc_comp'][1])
        ax7.grid(which='minor')
        self._draw_time_region(ax7)

        # Plot LWC
        self.lwcs.plot(style='k-', ax=ax8, legend=True)
        ax8.set_ylabel(r"g cm$^{-3}$")
        ax8.set_ylim(bottom=LIMS['lwc'][0], top=LIMS['lwc'][1])
        ax8.grid(which='minor')
        self._draw_time_region(ax8)

        if save:
            # Set the output file name
            fname = "pcasp_timeseries_" + os.path.splitext(self.filename)[0] + ".png"
            outfile = os.sep.join([save_dir, fname])
            print "Creating image: " + outfile
            plt.savefig(outfile, format='png')
        else:
            plt.show()

        
    def plot_histograms(self, save=False, save_dir=None, fancy=True):
        '''Create an timeseries plot of specific variables.
        
        Parameters::
        ----------
        save = boolean 
            True to save an output plot, False displays on screen.
        fancy = boolean
            If set to true uses ggplot style, otherwise default is used.
        '''
        if fancy:
            matplotlib.style.use('ggplot')
            
        # Set the time string format
        fmt = '%H:%M:%S'
        
        # Create figure/axes instance        
        fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(8, 11), \
                                sharey="row")
        for ax, fills in zip(axs.flat, self.dts):
            # Calculate start and end times
            stime = netCDF4.datetime(fills[0][0], fills[0][1], fills[0][2],\
                        fills[0][3], fills[0][4], fills[0][5])
            etime = netCDF4.datetime(fills[1][0], fills[1][1], fills[1][2],\
                        fills[1][3], fills[1][4], fills[1][5])
            # Make title from these times
            figtitle = stime.strftime(fmt) + ' - ' + etime.strftime(fmt)
            
            ##if self.pcasp[ 
            # Calculate Poisson error from Cai et al. AMT (2013), Equation A4
            yerr = self.pcaspf.between_time(stime, etime).sum(axis=0).divide(self.flow3s.between_time(stime, etime).sum(axis=0)**2).pow(0.5).mul(1./self.dlog10d)

            # Calculate the dNd(log10{D}) value for this time period, 
            # sum at each diameter (bin), and replace zero values with NaN
            dndlogp = self.pcaspf.between_time(stime, etime).sum(axis=0).divide(\
                           (self.flow3s.between_time(stime, etime).sum(axis=0) * self.dlog10d))
            # Create a pandas Series for plotting
            dndlogpf = pd.Series(data=dndlogp.replace(to_replace=0., value=np.nan), \
                                 index=self.pcasp_bin_mid, name=figtitle)
            xlabels = self.pcasp_bin_mid.astype('|S5')
            dndlogpf.plot(ax=ax, logx=True, logy=True, yerr=yerr, legend=False, x=xlabels)

            ax.set_title(figtitle)
            ax.set_xlim(left=0.1, right=10.)
            ax.set_ylim(bottom=0.01, top=10000.)
            ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.set_xlabel(r'D ($\mu$m)')
            ax.set_ylabel(r"dN/dlog$_{10}$D (cm$^{-3}$)")
            
        fig.tight_layout()

        if save:
            # Set the output file name
            fname = "pcasp_timehist_" + os.path.splitext(self.filename)[0] + ".png"
            outfile = os.sep.join([save_dir, fname])
            plt.savefig(outfile, format='png')
            print "Creating image: " + outfile
        else:
            plt.show()
            
def parse():
    '''
    Parse the input command line.
    Parameters::
    ----------
    argv - string
        Input command line string.
    Notes::
    -----
    Returns directory and field for initialization.
    '''
    parser = argparse.ArgumentParser(
              description="Start PCASP Check Software.")

    parser.add_argument('-v', '--version', action='version',
                        version='ARTview version %s' % (VERSION))

    # Directory argument now optional
    parser.add_argument('-f', '--file', type=str,
                        help="File to process",
                        default=None)
    parser.add_argument('-d', '--directory', type=str,
                        help="Directory where file is located", 
                        default=None)
    parser.add_argument('-p', '--project', type=str,
                        help="Project (e.g. pecan15)",
                        default=None)
    parser.add_argument('-s', '--save', type=str,
                        help="True to save False",
                        default=None)
    parser.add_argument('-i', '--imagepath', type=str,
                        help="Path to save output image",
                        default=None)

    # Parse the args
    args = parser.parse_args()#argv[1::])
    print args
    
    if args.file is None:
        print "ERROR: Please Choose file to process"
        print "For useage help: pcaspy_kingair.py -h"
        print "-------------------------------------"
        return
    if args.project is None:
        args.project = PROJECT
    if args.directory is None:
        args.directory = os.sep.join(["/kingair_data", args.project, "work"])
    # This next line saves an image if '-s' or '--save' has anything
    if args.save is None:
        args.save = False
    else:
        args.save = True
    # If an image path is provided
    if args.imagepath is None:
        args.imagepath = os.path.dirname(os.path.realpath(__file__))
    else:
        if (os.path.exists(args.imagepath)):
            args.save = True
        else:
            args.save = False
            print "Directory does not exist, cannot save file!"

    return args

#####################
## RUN THE PROGRAM ##
#####################
if __name__ == '__main__':
    # Call the parser get command line input
    cmdargs = parse()
    nc_filepath = os.sep.join([cmdargs.directory, cmdargs.file])

    # Start the check class
    an = CheckPCASP(nc_filepath, dts)
    an.plot_timeseries(cmdargs.project, save=cmdargs.save, save_dir=cmdargs.imagepath)
    if len(dts) > 0:
        an.plot_histograms(save=cmdargs.save, save_dir=cmdargs.imagepath)
    else:
        print "No subset times selected, therefore no histogram plot created"
