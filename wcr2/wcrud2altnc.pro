;$Id: wcrud2altnc.pro,v 1.1 2015/03/14 00:08:48 haimov Exp haimov $
;+
; NAME:
;       WCRUD2ALTNC
;
; PURPOSE:
;       Returns WCR up(beam1), down(beam1), or down&up(beam1&beam2) 
;       reflectivity(or velocity) resampled into the vertical plane and 
;       adjusted for the radar altitude. Deviations in the beam(s) data
;       in the horizontal plane are ignored (horizontal homogeneity is 
;       assumed). If beam(s) pointing angles are not given assume no 
;       deviation from vertical (straight leg, small attitude angles).
;
; CATEGORY:
;       wcrtools
;
; CALLING SEQUENCE:
;       result=wcrud2altnc(var,rng,alt,beam1,beam2)
;
; INPUTS:
;       Var:     fltarr(nrg,nprof); reflectivity or Doppler velocity for
;                 nrg range gates and nprof profiles; Note: var should not
;                 be rotated (input it as measured by the radar regardless
;                 of which beam it is, down, down-slant or up)
;       Rng:     fltarr(nrg), radar range in meters for Var; must be equally-
;                 spaced monotonically ascending vector (See SIDE EFFECTS)
;                for 'up' beam: Rng must be positive
;                for 'down' beam: Rng must be negative if BEAMDIR is not given
;                for 'down&up' combined: Rng must start with negative values 
;                for the 'down' beam followed by positive values for 'up' beam
;       Alt:     fltarr(nprof), radar altitude for every profile
;       Beam1:   fltarr(3,nporf), radar beam (up or down) unit vector
;                 in ground coordinates matching input data Var
;       Beam2:   fltarr(3,nprof), radar beam unit vector for the 'up' beam
;                 if Var is combined down&up beams; in this case Beam1 is the
;                 unit vector of the 'down' beam and Rng starts with the
;                 range for the down beam as negative and then the range for
;                 the 'up' beam as positive (e.g., -3075 to +3075 m)
;
; KEYWORD PARAMETERS:
;       ALTRNG:  output, fltarr(nalt), altitude vector in meters; 
;                        nalt - number of resampled altitudes
;       RELIEF:   input, when set (1) finds pixels set to Infinity and substitute
;                       them with 32767. This is used in conjunction with
;                       wcr2plotnc routine, where these pixels represent 
;                       subterrain and coloring of the terrain cross-section
;                       is requested; If relief=-1 - finds pixels that are
;                       outside the maximum radar range and substitute w/ 32767
;       BEAMDIR: input, string; beam direction: 'up','down','updown'
;                 if Beam1/2 input parameter(s) is/are given BEAMDIR is ignored.
;                 If Beam1/2 is/are not given and BEAMDIR is not define it is 
;                 determined by the sign of Rng
;       FLTLVL:  input, when set extends altitude by 250 m to include the  
;                 the aircraft flight level in the interpolated single beam
;                 data (flight level trace then can be added by a plotting
;                 routine, e.g., wcr2plotnc) or
;                input, long, altitude offset in meters (9 is not allowed)
;                   9 - does not remove leading/trailing rows of missing values)
;                
;       MISVAL:  input,float, missing value for Var (default -32767.)
;                      Note: -NaN and Infinity are not allowed as missing values
;       HELP:    shows this text
;
; OUTPUTS:
;       Returns array of (nalt,nprof) elements
;      
; SIDE EFFECTS:
;       1) Radar data is linearly interpolated to match the sampling altitudes.
;          If reflectivity is in dBZ no attempt to convert to Z is made before
;          interpolation (thus small mean eror/bias may be present). Missing 
;          values in the input and output are filled with NaNs. Thus the output
;          will have some interpolated values becoming NaNs when a missing
;          value participates in the interpolation.
;
;       2) When Beam1/2 are given:
;          a) The output is also projected to the vertical plane.  No attempt is 
;             made to discard data with "large" horizontal displacement (e.g.,
;             for 45 deg roll and 3 km range the radar data is displaced 2.1 km 
;             horizontally).
;          b) The sampling of the result preserves the range spacing in Rng.
;             Thus, if more than one data point falls into the vertical range
;             spacing the output is the mean value; if no value falls into
;             the output range missing value, NaN, is recorded.
;
;       3) There is one exception to the requirement for Rng to be equally 
;          spaced when Rng represents the combined 'up' and 'down' beam radar
;          ranges.  2 equal and larger range gate spacings are allowed for the 
;          first range gates above and below 0 m (they are in effect the WCR
;          range gate delay).  Example: ...,-145.,-105.,-75.,75.,105.,145.,...
;          or ...,-145.,-105.,-75.,0.,75.,105.,145.,... (where 0 m range is 
;          always missing value for all profiles in Var)
;
; NOTES:
;        
; EXAMPLES:
;
;-
; MODIFICATION HISTORY:
;       Written by:  Samuel Haimov, February 2008
;       S.H., Aug 09 - now handles when there are no valid points in the data
;       S.H., Oct 09 - added FLTLVL keyword
;       S.H., Jan 12 - added RELIEF keyword
;       S.H., May 12 - added FLTLVL=-9 option
;       S.H., Mar 14 - edited FLTLVL=-9 to FLTLVL=9 option (-9  was incorrect
;                      since only abs(FLTLVL) is used
;                      fixed a bug line 203 (beamdir='down') was after line 204
;                      (endelse)
;
; Copyright (C) 2007-present, Samuel Haimov, Dept. of Atmos. Sci., University of 
; Wyoming.  This software may be used, copied, or redistributed as long as it 
; is not sold and this copyright notice is reproduced on each copy made. This
; routine is provided as is without any expressed or implied warranties
; whatsoever.  Other limitations apply as described in the Readme file.
;

function wcrud2altnc, var,rng,alt,beam1,beam2,altrng=altrng,misval=misval,$
                      beamdir=beamdir,fltlvl=fltlvl,relief=relief,help=help

  on_error,2

; Check input <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
  if keyword_set(help) then begin
    doc_library,'wcrud2altnc'
    return,''
  endif

  nrg=n_elements(rng)
  rgs=abs(rng[nrg-1]-rng[nrg-2])
  if rgs lt 3. then begin
     message,'Invalid Rng vector',/info
     message,'Rng must be in meters with step greater than 3 m'
  endif

  if n_elements(fltlvl) eq 0 then flightlevel=0L $
  else if fltlvl eq 1 then flightlevel=250L $
  else flightlevel=long(abs(fltlvl))
  
  if n_elements(relief) eq 0 then relief=0

; ........................................................................... 
  if n_params() lt 3 then begin
    message,"Usage: res=wcrud2altnc(var,radrng,alt,beam1,beam2,$",/info
    message,"           altrng=altrng,beamdir=beamdir,misval=misval)",/info
    message,"Help:  print,wcrud2altnc(/help)",/info
    return,''
; ........................................................................... 
  endif else if n_params() eq 3 then begin
    if n_elements(beamdir) eq 0 then begin
      if rng[0] lt 0. and rng[nrg-1] gt 0. then begin
        message,"Beamdir not given. Assume 'updown' data", /info
        beamdir='updown'
      endif else if rng[nrg-1] lt 0. then begin    
        message,"Beamdir not given. Assume 'down' data", /info 
        beamdir='down'
      endif else if rng[nrg-1] gt 0. then begin    
        message,"Beamdir not given. Assume 'up' data", /info
        beamdir='up'
        flightlevel=-flightlevel
      endif
      radrng=rng
    endif else begin
      beamdir=strlowcase(beamdir)
      if beamdir eq 'up' then begin
        radrng=abs(rng)
        flightlevel=-flightlevel
      endif else if beamdir eq 'down' then begin
        if rng[nrg-1] lt 0 then begin
          if rng[nrg-1] ge rng[0] then $
            message,'Rng must be monotically ascending' $
          else radrng=rng
        endif else begin
          if rng[nrg-1] ge rng[0] then radrng=-rng $
          else message,'Rng must be monotically ascending'
        endelse     
      endif else if beamdir eq 'updown' then begin
        junk=0
        ind1=where(rng lt 0., compl=ind2)
        if ind1[0] eq -1 then begin
          message,'Rng must start with negative (maximum down) range.',/info
          junk=1
        endif
        if ind2[0] eq -1 then begin
          message,'Rng must have positive range values.',/info
        endif  else if junk then return,''
        radrng=rng
      endif
    endelse
; ........................................................................... 
  endif else if n_params() eq 4 then begin
    if beam1[2,0] ge 0 then begin
      radrng=abs(rng)
      flightlevel=-flightlevel
      beamdir='up'
    endif else begin
      if rng[nrg-1] lt 0 then begin
        if rng[nrg-1] ge rng[0] then $
          message,'Rng must be monotonically ascending' $
        else radrng=rng
      endif else begin
        if rng[nrg-1] ge rng[0] then radrng=-rng $
        else message,'Rng must be monotonically ascending'
      endelse     
      beamdir='down'        
    endelse
  endif else if n_params() eq 5 then begin
  ; Radrng must start with negative values and end with positive
    radrng=rng
    junk=0
    ind1=where(radrng lt 0., compl=ind2)
    if ind1[0] eq -1 then begin
      message,'Rng must start with negative (maximum down) range.',/info
      junk=1
    endif
    if ind2[0] eq -1 then begin
      message,'Rng must have positive range values for Beam2.',/info
      junk=1
    endif

  ; Beam1/2 must be down/up-pointing beams (beam1[2,k]<0./ beam2[2,k]>0.)

    ind1=where(beam1[2,*] gt 0.)
    ind2=where(beam2[2,*] lt 0.)
    if ind1[0] ne -1 then begin
      message,'Beam1 (down-pointing beam) contains positive z-components',/info
      junk=1
    endif  
    if ind2[0] ne -1 then begin
      message,'Beam2 (up-pointing beam) contains negative z-components'
    endif else if junk then return,''
    if total(beam1[2,*]) eq 0. or total(beam2[2,*]) eq 0. then $
      message,'Beam1 or Beam2  is 0',/info
    beamdir='updown'   
  endif

; Initialize <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  if n_elements(misval) eq 0 then misval=-32767.
  if finite(misval) eq 0 then nan=1 else nan=0

  if nan eq 0 then begin
    ind1=where(var eq misval)
    if ind1[0] ne -1 then var[ind1]=!values.f_nan
  endif

  nprof=n_elements(alt)
  minrng=min(radrng)
  maxrng=max(radrng)

  junk=size(var,/dim)
  if n_elements(junk) gt 2 then $
    message,'Input Var must be 2-D array' $
  else if junk[0] ne nrg or junk[1] ne nprof then $
    message,'Invalid input dimension size found.' 

; Calculate the output altitude vector
; (add/subtract the flight level for down/up images if requested)

  altmin=long(min(alt))+(flightlevel<0)
  junk=long(max(alt)+0.5)-altmin+abs(flightlevel)
  ind1=max((radrng-shift(radrng,1))[1:*])/rgs
  ind1=long(ind1)+(ind1 mod 1 gt 0)+long(junk/rgs)+((junk mod rgs) ne 0.)
  altrng=altmin+minrng+findgen(nrg+ind1)*rgs
  nalt=n_elements(altrng)

; Assign minimum deviation of the beam from vertical triggering correction

  minbdiv=0.97 ; 3% = 0.9 m (rgs=30m) ; ignore smaller diviations

; <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
; Interpolate data to account for beam-pointing angles and match altrng
; (no significant decimation is expected thus just interpolation, no filtering)  

  outvar=make_array(nalt,nprof,value=!values.f_nan)
  minalt=make_array(nprof,value=!values.f_nan)
  maxalt=minalt
  
  for k=0L,nprof-1 do begin  ; begin LOOP for all profiles
  ; ...........................................................................
  ; Interpolate data according to the beam(s) z-component projection;
  ; apply only if abs(beam[2,k]) lt 0.95 (>1.5 m for rgs=30m and tas=100m/s)

    if n_params() eq 5 then begin  
    ; up&down combined ........................................................

      out = var[*,k]

      if beam1[2,k] gt -minbdiv then begin ; Down beam
        ind2=where(radrng lt 0.)
        nind=n_elements(ind2)-1L
        radrnga=radrng[ind2]*beam1[2,k]
        ind1=(where(long(radrnga[0]) le long(reverse(-radrng[ind2]))))[0]
        if ind1 lt nind then begin
          junk=make_array(nind+1L,value=!values.f_nan)
          junk[nind-ind1:*]=congrid(var[ind2,k],ind1+1L, /interp)        
        endif else begin
          junk=interpol(var[ind2,k],radrnga,-radrng[ind2])
        endelse
        out[ind2]=junk
      endif

      if beam2[2,k] lt minbdiv then begin ; Up beam
        ind2=where(radrng gt 0.)
        nind=n_elements(ind2)-1L   
        radrnga=radrng[ind2]*beam2[2,k]
        ind1=(where(long(radrnga[nind]) le long(radrng[ind2])))[0]
        if ind1 lt nind then begin
          junk=make_array(nind+1L,value=!values.f_nan)
          junk[0:ind1]=congrid(var[ind2,k],ind1+1L, /interp)        
        endif else begin
          junk=interpol(var[ind2,k],radrnga,radrng[ind2])
        endelse       
        out[ind2]=junk
      endif

    endif else if n_params() eq 4 then begin
      if abs(beam1[2,k]) lt minbdiv then begin
      ; up or down beam only ....................................................

        radrnga=radrng*beam1[2,k] ; radrng projection to vertical

      ; resample to match range gate sampling (rgs)

        ind1=(where(long(radrnga[nrg-1]) le long(abs(radrng))))[0]

        if ind1 lt nrg-1L then begin
          out=make_array(nrg,value=!values.f_nan)
          out[0:ind1]=congrid(var[*,k],ind1+1L, /interp)        
        endif else begin
          out=interpol(var[*,k],radrnga,radrng)
        endelse    
    
      endif else out=var[*,k]
    endif else out=var[*,k]
 
  ; ...........................................................................
  ; Interpolate and shift data according to altitude

    minalt[k]=alt[k]+minrng
    maxalt[k]=alt[k]+maxrng
    ind1=where(altrng ge minalt[k] and altrng le maxalt[k])
    outvar[ind1,k]=interpol(out,alt[k]+radrng,altrng[ind1])

  endfor  ;  end LOOP for all profiles

; <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

; Remove gates with negative altitude if any (MSL of -100 m is used)

  ind1=where(altrng gt -100.,compl=ind2)
  if ind2[0] ne -1 then begin
    outvar=outvar[ind1,*]
    altrng=altrng[ind1]
    nalt=n_elements(ind1)
  endif

; Remove trailing columns (altitude gates) containing NaNs only if any      
; (will use while (slow); not a biggy)    

  if flightlevel le 0 and flightlevel ne -9 then begin
    break=1B
    k=nalt-1L
    while break do begin
      if total(finite(outvar[k,*])) ne 0 then break=0B $
      else k=k-1L
      if k lt 0 then begin
        message,'Check your input. No valid data in the output.',/info
        return,outvar
      endif
    endwhile
    outvar=outvar[0:k,*]
    altrng=altrng[0:k,*]   
  endif

; Remove leading columns (altitude gates) of NaNs only if any  

  if flightlevel ge 0 and flightlevel ne 9 then begin
    break=1B
    k=0L
    while break do $
      if total(finite(outvar[k,*])) ne 0 then break=0B $
    else k=k+1L
    outvar=outvar[k:*,*]
    altrng=altrng[k:*,*]
  endif

  if relief eq 1 then begin
  ; Find all pixels set to Infinity or -NAN after interpolation calculations)
  ; and replaces them with 32767. These are sub-surface pixels left in the image
  ; to be colored in a plotting routine.For this to work only sub-surface pixels 
  ; must be equal to Infinity (done in wcr2plotnc). 
    ind1=where(finite(outvar,/infinity) eq 1 or finite(outvar,/nan,sign=-1) eq 1)
    if ind1[0] ne -1 then outvar[ind1]=32767.
  endif else if relief eq -1 then begin
  ; Find all pixels that are out of radar range and mark as 32767.
    nalt=n_elements(altrng)
    if beamdir eq 'updown' or beamdir eq 'down' then begin
      ind1=value_locate(altrng,minalt)
      ind2=where(ind1 gt 0)
      if ind2[0] ne -1 then $ ; mark all gates between min(altrng) and minalt
        for k=0,n_elements(ind2)-1 do outvar[0:ind1[ind2[k]],ind2[k]]=32767.    
    endif
    if beamdir eq 'updown' or beamdir eq 'up' then begin         
      ind1=value_locate(altrng,maxalt)
      ind2=where(ind1 lt nalt-1)
      if ind2[0] ne -1 then $  ; mark all gates between maxalt and max(altrng)
        for k=0,n_elements(ind2)-1 do outvar[ind1[ind2[k]]:*,ind2[k]]=32767.
    endif    
  endif

; <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
; Hopefully, we are done

  return,outvar

end  
       
