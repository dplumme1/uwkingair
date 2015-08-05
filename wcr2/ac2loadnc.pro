;$Id: ac2loadnc.pro,v 1.1 2010/11/15 19:40:32 haimov Exp $
;+
; NAME:
;       AC2LOADNC
;
; PURPOSE:
;       Returns KingAir or NCAR C130 aircraft data variable from a netCDF file.
;       There is only one time dimension in the returned variable.  
;       See the notes.  
;
; CATEGORY:
;       wcrtools
;
; CALLING SEQUENCE:
;       result=ac2loadnc(Cdfid,varname [,filename],close=close,sps=sps, $
;              startt=startt,stopt=stopt,stime=stime,help=help)
;
; INPUTS:
;       Cdfid:    netCDF file ID, returned from a previous call to 
;                 ncdf_open, ncdf_create, ncdf_varload or ncdf_fcopy
;       Varname:  string, netCDF variable name
;                 ATTN: See STIME keyword and NOTES/EXAMPLES on how to load
;                 the time vector
;       Filename: string, optional netCDF file name; if given will open netCDF 
;                 file first. In this case Cdfid is an output parameter
;
; KEYWORD PARAMETERS:
;       CLOSE:  when set AC netcdf file will be closed.
;       SPS:    output, integer; samples per second for the loaded variable
;               input, integer; samples per second for the time variable only,
;                      requires STIME keyword to be set (time variable does not
;                      have sps dimension and SPS as input can be used to
;                      resample time to include fractional seconds)
;       STARTT: input, long, start time in hhmmss
;       STOPT:  input, long, stop time in hhmmss including last second record
;       STIME:  input, when set identifies varname as a time variable; this is
;                      needed to generate time in unix seconds; otherwise the
;                      time is recorded in seconds from the epoch defined in
;                      time variable units attribute
;                      ATTENTION: attribute units holds information about the
;                      epoch in the following strictly enforced format:
;                      'seconds since 20YY-MM-DD HH:MM:SS +0000'
;
;       HELP:   shows this text
;
; OUTPUTS:
;       Cdfid:  netCDF file ID (if defined), if filename given or returned from
;               a previous call to ncdf_open, ncdf_create, ncdf_varload, 
;               ncdf_fcopy or other routine
;
; SIDE EFFECTS:
;       The function returns -1 if the variable is not found or an empty string
;       if the usage is incorrect.
;       Variables cannot have more than 4 dimensions.
;
; NOTES:
;       KingAir netCDF files follow the NCAR-RAF/nimbus convention.  This
;       convention requires that the time dimension, 'time' (unlimited or not), 
;       is the first (last for idl, which uses column/fortran order) and the 
;       next dimension, if exists, is always samples per second, 'sps'.
;
;       This function returns the requested variable with the 'sps' dimension,
;       when exists, merged with the 'time' dimension.  Thus the output always
;       has only one time dimension and is sampled at the rate that can be 
;       extracted using the SPS keyword
;
;       If startt and stopt are given netcdf variable is extracted within the
;       specified time. If stopt is not defined the variable is read untill
;       the end of the file.  And if both, startt and stopt are not defined or
;       are invalid the complete data set for the variable is read.
;       
;       Loading a time variable:
;         Time variable is typically converted to UNIX seconds if it is not
;         recorded in UNIX seconds (e.g., NCAR Time is not in UNIX seconds)
;
;         The following basic time variables are available from the UWKA and
;         the NCAR C130 data files:
;           time        - long array of seconds from the epoch defined in
;                         units attribute (N2UW)
;           Time        - same as time (N130AR)
;         See the examples below on how to get time in unix seconds from 
;         time or Time
;        
;         For 'time' or 'Time' variable, 'units' attribute is expected in the 
;         the  following format: 'seconds since 20YY-MM-DD HH:MM:SS +0000'
;        
; EXAMPLES:
;
;      1) The following example assumes that a netcdf file was previously opened,
;         hroll variable exists and 121915-122252 is a legitimate time (hhmmss)
;         subinterval.
;
;         roll=ac2loadnc(ncid,'hroll',sps=sps,start=121915,stop=122252)
;
;         Note that sps is an output in this call and any input value for sps
;         would be ignored.
;      
;      2) Read the default time var for the whole data file.
;         Note that the output inlcudes whole seconds only.
;
;         time=ac2loadnc(ncid,'',/stime) 
;         time=ac2loadnc(ncid,'time')  ; the same as above but if 'time' does
;                                        not exist the routine will return -1
;
;      3) As 2) but resample time to include fractional seconds based on sps
;
;         time=ac2loadnc(ncid,'time',sps=25,/stime)  
;         time=ac2loadnc(ncid,'time',sps=25) ; no resampling, sps input ignored
;
;      4) Load 'Time' from NCAR C130 data set (nimbus revision 1.3 or later)
;
;         time=ac2loadnc(ncid,'',/stime)     ; load time in UNIX seconds starting
;                                            ; from the 'units' starting time
;         time=ac2loadnc(ncid,'Time')        ; load time in sec starting from 0.
;         time=ac2loadnc(ncid,'Time',/stime) ; load time in sec; Time is
;                                            ; converted to UNIX seconds using
;                                            ; 'units' attribute
; 
;-
; MODIFICATION HISTORY:
;       Written by:  Samuel Haimov, September 2009
;
; Copyright (C) 2009-present, Samuel Haimov, Dept. of Atmos. Sci., University of 
; Wyoming.  This software may be used, copied, or redistributed as long as it 
; is not sold and this copyright notice is reproduced on each copy made. This
; routine is provided as is without any expressed or implied warranties
; whatsoever.  Other limitations apply as described in the Readme file.
;

function ac2loadnc, Cdfid,varname,filename, close=close, $
                   sps=sps,startt=startt,stopt=stopt,stime=stime,help=help

  on_error,2
 
  if keyword_set(help) then begin
    doc_library,'ac2loadnc'
    return,''
  endif

  if n_params() lt 2 then begin
    message,"Usage: res=ac2loadnc(Cdfid,varname [,filename], $",/info
    print,"          close=close,sps=sps,startt=startt,stopt=stopt,help=help)"
    return,''
  endif else if n_params() eq 3 then begin
    if n_elements(Cdfid) eq 1 then ncdf_close,Cdfid
    Cdfid=ncdf_open(filename)
  endif  

; Inquire about varname var.  If does not exist return -1

  if varname ne '' then begin
    varinq=ncdf_varinq1(Cdfid,varname)
    if varinq.id eq -1 then return,-1

  ; If varname var is a scalar load it and exit

    if varinq.ndims eq 0 then begin
      ncdf_varget,Cdfid,varinq.id,vardata
      sps=0
      return,vardata
    endif
  endif

; Extract time and determine the time range for the variable if needed

  if n_elements(startt) gt 0 or n_elements(stopt) gt 0 or keyword_set(stime) $
  then begin
    platform=ncdf_attload(Cdfid,'Platform',/global,/silent) ; '-2' if none
    if keyword_set(stime) and varname ne '' then begin
      time=ncdf_varload(Cdfid,varname)
      if time[0] eq -1 then message,'Incorrect time var name.'
      junk=ncdf_attload(Cdfid,varname,'units')
      ; expected 'units' format: seconds since 20YY-MM-DD HH:MM:SS +0000'
      ind=strpos(junk,'20')
      junk=unixtime2us(ymd2date(strmid(junk,ind,4),strmid(junk,ind+5,2),$
        strmid(junk,ind+8),format='w$ n$ 0d$')+' '+strmid(junk,ind+11,8)+$
        ' '+strmid(junk,ind,4))
      time=double(time)+junk
    endif else begin  ; use default varname for time
      varnames=ncdf_varnames(Cdfid)
      if platform eq 'N130AR' then ind1=where(varnames eq 'Time') $
      else ind1=where(varnames eq 'time')
      if ind1[0] eq -1 then message,"'time' or 'Time' do not exist." $
      else time=ncdf_varload(Cdfid,varnames[ind1[0]])
      junk=ncdf_attload(Cdfid,varnames[ind1[0]],'units')
      ; expected 'units' format: seconds since 20YY-MM-DD HH:MM:SS +0000'
      ind=strpos(junk,'20')
      junk=unixtime2us(ymd2date(strmid(junk,ind,4),strmid(junk,ind+5,2),$
           strmid(junk,ind+8),format='w$ n$ 0d$')+' '+strmid(junk,ind+11,8)+$
           ' '+strmid(junk,ind,4))
      time=double(time)+junk
    endelse
    
    if n_elements(startt) gt 0 or n_elements(stopt) gt 0 then begin
      timeh=long(sec2time(time mod 86400))
      if n_elements(startt) eq 0 then ind1=0 else begin
         ind1=where(long(startt) eq timeh)
         if ind1[0] eq -1 then begin
            ind1=0
            message,'Start time not valid. Load from the beginning.',/info 
         endif else ind1=ind1[0]
      endelse

      if n_elements(stopt) gt 0 then begin
          ind2=where(long(stopt) eq timeh)
          if ind2[0] eq -1 then begin
            ind2=n_elements(timeh)-1 
            message,'Stop time not valid. Load to the end of the file.',/info 
          endif else ind2=ind2[n_elements(ind2)-1]
      endif else ind2=n_elements(timeh)-1
      n=ind2-ind1+1
      if n lt 1 then message,'Incorrect start/stop times'
    endif else begin
      ind1=0
      ind2=n_elements(time)-1
    endelse

    ; If time var is requested and sps > 1 rebin time vector to include 
    ; the fractional seconds end exit

    if keyword_set(stime) then begin
      if n_elements(sps) eq 1 && sps gt 1 then begin
        timeh=time[ind1:ind2]
        time=rebin(double(timeh),sps*n_elements(timeh))
      ; take care of the last second (it is not interpolated by rebin)
        time[n_elements(time)-sps:*]=time[n_elements(time)-sps:*]+findgen(sps)/sps
        return,time
      endif else return, time[ind1:ind2]
    endif

  endif else begin 
     ind1=0 
     n=varinq.dimsize[varinq.ndims-1]
     ind2=n-1
  endelse

; Determine offset and count

  offset1=ind1
  count1=n 
  sps=1
  if varinq.ndims eq 1 then begin
    offset=offset1
    count=count1
    sps=1
  endif else begin
    if strmid(varinq.dimname[varinq.ndims-2],0,3) ne 'sps' then begin
       message,'WARNING: '+varname+' 2nd dimension is not sps',/info
       print,'          Dimension name: ',varinq.dimname[varinq.ndims-2]
    endif
    sps=varinq.dimsize[varinq.ndims-2]
 
    if varinq.ndims eq 2 then begin
      offset=[0,offset1]
      count=[varinq.dimsize[0],count1]
    endif else if varinq.ndims eq 3 then begin
      offset=[0,0,offset1]
      count=[varinq.dimsize[0:1],count1]  
    endif else begin
      offset=[0,0,0,offset1]
      count=[varinq.dimsize[0:2],count1]  
    endelse
  endelse

; Load the variable

  ncdf_varget,Cdfid,varinq.id,vardata,offset=offset,count=count
        
; Reform

  if varinq.ndims gt 1 then begin
    if varinq.ndims eq 2 then  $                                      ; 2-D var
      vardata=reform(vardata,sps*n) $
    else if varinq.ndims eq 3 then  $                                 ; 3-D var
      vardata=reform(vardata,varinq.dimsize[0],sps*n) $
    else if varinq.ndims eq 4 then  $                                 ; 4-D var
      vardata=reform(vardata,varinq.dimsize[0],varinq.dimsize[1],sps*n)
  endif

  if n_elements(vardata) gt 1 then vardata=reform(vardata) $
  else vardata=vardata[0]

  if keyword_set(close) then begin 
    ncdf_close, Cdfid
    undefine, Cdfid
  endif
 
return,vardata

end  
       
