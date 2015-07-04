#! /usr/bin/env python
#####################
# flthours.py - Create plot of flight hoursList of each flight
#####################

# import required packages
import os
import argparse
from datetime import datetime
import pandas as pd
import numpy as np

from bokeh.plotting import show, save, figure, output_file, ColumnDataSource
from bokeh.models import HoverTool
from collections import OrderedDict

##################################
# DEFINE GLOBAL VARIABLE #
##################################
# Project information: name, total hours, and dates

PROJ = 'pecan15'
PROJHRS = 120
FIRSTDATE = '2015-06-01'
LASTDATE  = '2015-07-15'
PROJDAYS = 45
##################################

class FltHrs(object):
    """
    A class mentod to retrieve and plot metar data.
    
    The metar retrieval is a modification of code found at:
    https://github.com/akrherz/iem/blob/master/scripts/asos/iem_scraper_example.py
    
    It is dependent upon the Iowa State Mesonet database.
    """
    def __init__(self, args):
        '''Initialize the class'''

        # Set date formats to be used with datetime
        self.d_fmt = "%Y-%m-%d %H:%M"

        # Use passed arguments
        self.project = args.project
        self.hours = args.hours

        # Set the URL to pull data
        self.url = "http://flights.uwyo.edu/projects/" + self.project + "/"

        # Create an date-time array for project hours
        self.start = datetime.strptime(args.start, "%Y-%m-%d")
        self.end = datetime.strptime(args.end, "%Y-%m-%d")
        #Calculate the number of days in project
        self.projdays = (self.end - self.start).days

        rng = pd.date_range(args.start, periods=self.projdays, freq='D')
        self.hrsrng = pd.Series(np.linspace(0., self.hours, num=self.projdays, endpoint=True), \
                                  index=rng)

    def get_web_data(self):
        '''Function to return flight data from web'''
        # Load in the URL and parse
        self.df = pd.read_html(self.url, header=1, parse_dates=True, index_col=0)

        # Sum the research, calibration, and test hours
        self.reshrs = self.df[1]['Hours'].sum()
        self.calhrs = self.df[2]['Hours'].sum()
        self.testhrs = self.df[3]['Hours'].sum()

        # Calculate the total hours in the field
        self.fieldhrs = self.reshrs + self.calhrs

        # Create arrays for research and calibration flights
        self.resflt = pd.concat([self.df[1],self.df[2]])
        self.calflt = self.df[2]
        self.resflt.sort(ascending=True, inplace=True)
        self.calflt.sort(ascending=True, inplace=True)

    def write_json(self, loc=None):
        '''Create a JSON text output file with total flight hours to date'''
        # Full path of JSON file
        if loc is None:
            webpath = os.sep.join(['/home/webadmin','flights','projects',self.project,'flthours.js'])
        else:
            webpath = os.sep.join([loc,'flthours.js'])
        nowdate = datetime.now()

        # Create string
        strout = "var flthours = 'As of %s, %.1f out of %d research hours were flown, "\
                 "%.1f remain.';\n"%\
                 (nowdate.strftime("%b %d, %Y"), self.fieldhrs, \
                  self.hours, self.hours-self.fieldhrs)

        # Open and write json file
        js = open(webpath, 'w')
        js.write(strout)
        js.write("var testhours = 'Test: %.1f';"%(self.testhrs))
        js.close()

    def plot_flight_hours(self, width=500, height=500, \
                          ln_wd=3, sz=10, fDir=".", fName='uwka_flthrs'):
        '''
        Create an output plot of research project flight hours
        
        Parameters::
        ----------
        width : int
            Height of plot in pixels.
        height : int
            Width of plot in pixels.
        ln_wd : int
            Line width to use for plots.
        sz : int
            Size of markers.
        fDir : str
            Output directory name to save file [default to current directory]
        fName : str
            Output filename [default 'meteogram']
        '''

        filepath = os.sep.join([fDir, fName + "_" + self.project +".html"])
        pagetitle = "UWKA Flight Hours"

        output_file(filepath, title=pagetitle)
        TOOLS = "pan,wheel_zoom,box_zoom,reset,save,hover"
        TITLE = "UWKA %s: Flight Hours"%(self.project)

        # Create a data source for 
        sourceres = ColumnDataSource(
                  data=dict(
                  Date=self.resflt.index.format(),
                  Flight=self.resflt['Flight#(*.kml)'].values,
                  Hours=self.resflt['Hours'].values,
                       )
                 )

        sourcecal = ColumnDataSource(
                  data=dict(
                  Date=self.calflt.index.format(),
                  Flight=self.calflt['Flight#(*.kml)'].values,
                  Hours=self.calflt['Hours'].values,
                       )
                 )

        self.hrsrng.plot()
        p = figure(x_axis_type = "datetime", tools=TOOLS, width=width, height=height)
        p.line(self.hrsrng.index.values, self.hrsrng.values, color='black', \
                legend='Allocated', line_width=ln_wd)
        p.circle(self.resflt.index, self.resflt['Hours'].cumsum().values, color='red', \
                  legend='Cumulative', size=sz, source=sourceres)
#        p.circle(self.calflt.index, self.calflt['Hours'].cumsum().values, color='blue', \
#                  legend='Cumulative', size=sz, source=sourcecal)
        p.title = TITLE
        p.grid.grid_line_alpha=0.3
        p.xaxis.axis_label = 'Date'
        p.yaxis.axis_label = 'Hours'
        p.legend.orientation = "top_left"

        hover = p.select(dict(type=HoverTool))
        hover.tooltips = OrderedDict([
            ("Cum. Hours", "@y"),
            ("Date", "@Date"),
            ("Flight", "@Flight"),
            ("Hours","@Hours"),
            ("index", "$index"),
        ])

        show(p)
        #save(obj=p, filename=filepath, title=pagetitle)

#####################
## RUN THE PROGRAM ##
#####################
if __name__ == '__main__':
    """
    NOTES::
    -------
    To run:
        $ flthours.py
    """

    # Define the command line structure
    parser = argparse.ArgumentParser(description="Create a flight hours plot for project")

    igroup = parser.add_argument_group(
            title="Set the station, begin and end date",
            description=("Run this executable to produce "
                        "an html plot for project flight hours. "
                        "Global environmental variables can be"
                        "set in the script (at top)."
                        "Optional arguments below can be passed"
                        "to override the global variables."
                        " "))

    igroup.add_argument('-p', '--project', type=str, default=None, 
                              help='project (e.g. pecan15)')

    igroup.add_argument('-t', '--hours', type=str, default=None, 
                              help='project (e.g. 120)')

    igroup.add_argument('-s', '--start', type=str, default=None, 
                              help='Start date to search (e.g. "2015-6-1")')
                       
    igroup.add_argument('-e', '--end', type=str, default=None, 
                              help='Start date to search (e.g. "2015-7-15")')

    # Parse the args
    args = parser.parse_args()
    if args.project is None:
        args.project = PROJ
    if args.hours is None:
        args.hours = PROJHRS
    if args.start is None:
        args.start = FIRSTDATE
    if args.end is None:
        args.end = LASTDATE

    # Start the Metar class
    flight = FltHrs(args)

    # Load the data
    flight.get_web_data()

    # Create a json file
    flight.write_json(loc=".")

    # Plot a meteogram
    flight.plot_flight_hours(fDir=".")