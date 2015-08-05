Display
=======

This directory contains display software that can be used to view data collected
by the University of Wyoming King Air research aircraft. 

##wcr2l2display.sav 

GUI tool written in the IDL language to display WCR2L2 data.
An IDL license isnot needed to run this tool. 
The free [IDL Virtual Machine](http://www.exelisvis.com/docs/creating_save_files_of_p.html#Save_3439793193_858495) is all that needed. 

A brief description of the software:
```
Controls:
New - Open a dialog_pickfile to select a new WCR file.
Destroy - Terminates the program and destroys the widget.
Quit - Terminates the program without destroying the widget.
Export - Export current display to postscript file.
#Image - Number of image to display.
Help - Display help window.

Color:
Reflectivity - Allow change in the color table for reflectivity.
Velocity - Allow change in the color table for velocity.

Min/Max Refl. - Allow change in the min/max value of reflectivity to be scaled for the display.

Min/Max Vel. - Allow change in the min/max value of velocity to be scaled for the display.

Time - Allow change in the time range (UTC).

Alt. - Allow change in the altitude range (km).

Mask:
Off6 - Remove beam more than 6 degee off of vertical.
Subsurface - Highlight subsurface gates.
Out of Range - Highlight gates outside max range.
Surface Clutter - Remove gates affected by surface clutter.
Transmitter Leakage - Remove gates affected by transmitter leakage.
```

Contact:
Dr. Yonggang Wang (WYG@uwyo.edu)
Dr. Bart Geerts (Geerts@uwyo.edu)