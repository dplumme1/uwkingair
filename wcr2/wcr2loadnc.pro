;$Id: wcr2loadnc.pro,v 1.2 2015/03/14 00:08:48 haimov Exp haimov $
;+
; NAME:
;       WCR2LOADNC
;
; PURPOSE:
;       Load variables from a WCR netcdf file.
;
; CATEGORY:
;       wcr2tools
;
; CALLING SEQUENCE:
;       result=wcr2loadnc(ncid,varname [,filename],startt=startt,stopt=stopt)
;
; INPUTS:
;       ncid:     netCDF variable ID, returned from a previous call to 
;                 ncdf_open, ncdf_create, ncdf_varload or ncdf_fcopy
;       Varname:  string, netCDF variable name
;       Filename: string, optional netCDF file name; if given will open netCDF 
;                 file first. In this case ncid is an output var
;
; KEYWORD PARAMETERS:
;       STARTT:   input, long, start time in hhmmss; (default is base_time or
;                   time[0], see Notes); if beyond the time interval of the file
;                   starts from the beginning
;       STOPT:    input, long, stop time in hhmmss, inclusive the last second 
;                   (default is the end if the file time, see Notes); if beyond
;                   the time interval of the file, read to the end of the file
;       HELP:     show this text ( use print,wcrloadnc(/help) )
;
; OUTPUTS:
;       ncid:     netCDF variable ID, returned if Filename is given
;
; SIDE EFFECTS:
;       The function returns an empty string if the usage is incorrect or 
;       -1 if the variable is not found or
;       -3 if 'time' does not exist.
;       -4 if STARTT and STOPT are not valid for this file
;
; NOTES:
;       The use of STARTT and STOPT requires the existence of netcdf variables
;       'time', which is the time stamp vector in unix seconds (seconds after
;       00:00:00 Jan 1, 1970).  The time dimension for all variables having it
;       is called 'profile'.  If the variable doesn not use 'profile' dimension
;       STARTT and STOPT are ignored. 
;
; EXAMPLES:
;
;       reflectivity=wcr2loadnc(ncid,'reflectivity') ;
;-
; MODIFICATION HISTORY:
;       Written by:  Samuel Haimov, November 2009
;       (this is a modified version of wcrloadnc, which takes care of
;        the variables format change between WCR1 and WCR2)
;       S.H., Dec 2011 - fixed a bug in count and offset for 2-D var
;       S.H., Jan 2013 - modified to work with L2 files when startt/stopt
;                        are defined.
;       S.H., Apr 2014 - fixed a bug in loading multi-dim vars when startt
;                        and/or stopt keywords are defined
;
; Copyright (C) 1996-present, Samuel Haimov, Dept. of Atmos. Sci., University of 
; Wyoming.  This software may be used, copied, or redistributed as long as it 
; is not sold and this copyright notice is reproduced on each copy made. This
; routine is provided as is without any expressed or implied warranties
; whatsoever.  Other limitations apply as described in the Readme file.
;

function wcr2loadnc,ncid,varname,filename,startt=startt,stopt=stopt,help=help

;  on_error,2

  if keyword_set(help) then begin
    doc_library,'wcr2loadnc'
    return,''
  endif

  if n_params() lt 2 then begin
    message,"Usage: res=wcr2loadnc(ncid,varname [,filename])",/info
    return,''
  endif else if n_params() eq 3 then ncid=ncdf_open(filename)

; Inquire about the varname

  vari=ncdf_varinq1(ncid,varname)
  if vari.id eq -1 then return,-1              ; varname doesn't exist

; Check if startt and stopt are given and is possible to use them

  indprof=(where(vari.dimname eq 'profile'))[0]
  if indprof eq -1 then begin
    strt=-1
    stpt=-1
  endif else begin
    if n_elements(startt) eq 0 then strt=-1 else strt=startt
    if n_elements(stopt) eq 0 then stpt=-1 else stpt=stopt
  endelse

  if strt ne -1 or stpt ne -1 then begin
    junk1=ncdf_varinq1(ncid,'time')
    if junk1.id eq -1 then return,-3     ; time doesn't exist
    time=ncdf_varload(ncid,'time')

    n=n_elements(time)

    if strt eq -1 then strt=0L else begin
       junk1=us2unixtime(time[0],/date,/utc)         
       junk2=long(time2sec(strt,unix=junk1))
       ind=(where(junk2 le long(time)))[0]
       if ind le 0 then begin
         junk1=us2unixtime(time[n-1],/date,/utc)         
         junk2=long(time2sec(strt,unix=junk1))
         strt=(where(junk2 le long(time)))[0]
       endif else strt=ind  
       if strt le 0 then strt=0L
    endelse

    if stpt eq -1 then stpt=n-1L else begin
       junk1=us2unixtime(time[n-1],/date,/utc)         
       junk2=long(time2sec(stpt,unix=junk1))
       ind=(where(junk2 le long(time)))[0]
       if ind le 0 then begin
         junk1=us2unixtime(time[0],/date,/utc)         
         junk2=long(time2sec(stpt,unix=junk1))
         stpt=(where(junk2 le long(time)))[0]
       endif else stpt=ind  
       if stpt le 0 then strt=n-1L
    endelse

    if stpt le strt then return,-4   ; invalid start/stop times

    ; Define offset and count: up to 4 dimensions are allowed
    ; (currently L1 and L2 have maximum of 3 dimensions)

    offset1=strt
    count1=stpt-offset1+1
    if vari.ndims eq 1 then begin
      offset=offset1
      count=count1
    endif else if vari.ndims eq 2 then begin
      if indprof eq 0 then begin
        offset=[offset1,0]
        count=[count1,vari.dimsize[1]]
      endif else begin
        offset=[0,offset1]
        count=[vari.dimsize[0],count1]
      endelse
    endif else if vari.ndims eq 3 then begin
      if indprof eq 0 then begin
        offset=[offset1,0,0]
        count=[count1,vari.dimsize[1:2]]
      endif else if indprof eq 1 then begin
        offset=[0,offset1,0]
        count=[vari.dimsize[0],count1,vari.dimsize[2]]
      endif else begin
        offset=[0,0,offset1]
        count=[vari.dimsize[0:1],count1]
      endelse
    endif else begin
      if indprof eq 0 then begin
        offset=[offset1,0,0,0]
        count=[count1,vari.dimsize[1:3]]
      endif else if indprof eq 1 then begin
        offset=[0,offset1,0,0]
        count=[vari.dimsize[0],count1,vari.dimsize[2:3]]
      endif else if indprof eq 2 then begin
        offset=[0,0,offset1,0]
        count=[vari.dimsize[0:1],count1,vari.dimsize[3]]
      endif else begin
        offset=[0,0,0,offset1]
        count=[vari.dimsize[0:2],count1] 
      endelse 
    endelse

    ; Load varname

    ncdf_varget,ncid,vari.id,tmp,offset=offset,count=count   

  endif else begin
    ncdf_varget,ncid,vari.id,tmp   
  endelse

  if n_elements(tmp) gt 1 then tmp=reform(tmp) $
  else tmp=tmp[0]

  return,tmp 

end
