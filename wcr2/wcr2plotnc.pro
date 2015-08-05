;$Id:
;+
; NAME:
;       WCR2PLOTNC
;
; PURPOSE:
;       Display reflectivity and Doppler velocity for WCR2 L1 netcdf file(s). 
;
; CATEGORY:
;       wcr2tools
;
; CALLING SEQUENCE:
;       wcr2plotnc,minmax,acwind=acwind,altbeam=altbeam,aspect=aspect,$
;                  asterdir=asterdir,breaki=breaki,colortb=colortb,fltlvl=fltlvl,$
;                  indir=indir,ldr=ldr,leak=leak,mark=mark,nav=nav,noerase=noerase,$
;                  ppz=ppz,profile=profile,relief=relief,satmark=satmark,$
;                  savedir=savedir,saveimg=saveimg,startt=startt,stopt=stopt,terrain=terrain,$
;                  thresh=thresh,zoomimg=zoomimg,zsrf=zsrf, help=help
;
; INPUTS:
;       Minmax:  fltarr(3),[dvminmax,zmax,zmin]
;                dvminmax - float, Doppler velocity min/max limit in m/s
;                zmax     - float, Reflectivity max limit in dBZ
;                Use -999 to use default values for parameters that you do not
;                want to set
;
; KEYWORD PARAMETERS:
;       ACWIND:  input, when set add the AC measured wind component into the 
;                radar beam to the velocity plot. For up-down merged beams, 
;                the blind area us used for plotting the AC wind; If a single 
;                beam is plotted FLTLVL should be set in order to plot the wind
;       ALTBEAM: input, integer, add radar altitude and optionally interpolate
;                to strictly vertical plane according to the radar beam-pointing
;                vectors (ignored for side beams)
;                0 - no action (default)
;                1 - interpolate to radar altitude
;                2 - interpolate to altitude and strictly vertical plane 
;                Notes: 1) If Netcdf file ALT variable is missing altbeam is 
;                          ignored. 
;                       2) The reflectivity interpolation is done in dBZ not Z
;       ASPECT:  float, profile-to-range ratio; use to average more profiles than
;                the actual aspect ; default use the actual aspect to display
;                1:1 flighttrack-to-range scale (this keyword is useful for very
;                large files that cannot be displayed with 1:1 scale)
;       ASTERDIR:input, string, path to ASTER DEM data directory
;                (ignored if TERRAIN keyword is not set)
;       BREAKI:  input, when set and the images are displayed, the routine will
;                 stop after every image plot; use .continue to continue
;                 plotting; if no more stops are needed assign 0 to BREAKI
;       COLORTB: intarr(2), [Zcolortb,DVcolortb], color tables for Z and DV
;                (default is [41,51] or [33,17] if 41/51 not available)
;       FLTLVL:  input, when set extends altitude by 250 m to include a trace  
;                 of the aircraft flight level (dashed line) in the image
;                 when ALTBEAM is set; or
;                input, long, desired altitude offset in meters
;       INDIR:   input, string, input directory for the WCR L1 netcdf files
;                 if given automatically loads all WCR L1 nc files
;                 (if not defined or '' opens dialog_pickfile)
;       LDR:     input, when set plot LDR instead of cross-pol reflectivity
;       LEAK:    input, when set(leak=1) and THRESH=9 removes pixels set in
;                 the reflectivity_mask as affected by leak
;       MARK:    input, string; text to be added at the top of the page
;                for ps/png/pdf (default is 'QUICK LOOK')
;       NAV:     input, average profiles and/or range gates before displaying
;                for increased sensitivity (requires THRESH keyword to be set)
;                input, npa, short integer, number of profiles to be avaraged
;                            (default is 1)
;                [npa,nrga], intarr(2), number of profiles and range gates 
;                            to average (default for nrga is 1)
;                ATTENTION: NAV is ignored if THRESH is not given or THRESH=9
;
;       NOERASE: when set will not erase the widgets or the graphic windows
;       PPZ:     when set display pulse-pair reflectivity instead of reflectivity
;       PROFILE: input, when set top axes is labeled profile/record number
;                instead of air distance
;       RELIEF:  input, when set, sub-terrain pixels for the reeflectivity 
;                       images are set to 32767 and plotted with the foreground
;                       color. Ignored if ALTBEAM=0 or THRESH ne 9. 
;       SATMARK: when set plot marks for partially or fully saturated pixels 
;                on the reflectivity images (ignored if NAV > 1 and/or zoom)
;       SAVEDIR: input, string, path/dir_name of directory to save images
;                       (default is <data_dir>/quicklook); ignored if
;                       SAVEIMG is 0 or not defined
;       SAVEIMG: input. [type,ori,res]; write a ps, pdf, or  png file   
;                type to a directory called quicklook:
;                 0 - does not write a ps/pdf/png files (default)
;                 1 - display the file and write ps file(s)
;                 2 - display the file and converts ps to png file(s)
;                 3 - display the file and converts ps to pdf file(s)
;                 If input is negative images are saved but not displayed
;                ori:
;                 0 - portrait
;                 1 - landscape(default) 
;                res: png resolution (number between 50 and 200, 
;                            default 200)
;       STARTT:  input, long, start time in hhmmss (see SIDE EFFECTS)
;       STOPT:   input, long, stop time in hhmmss, inclusive the last second 
;       TERRAIN: input, integer,remove/mark reflectivity pixels assoc. w/surface  
;                 0-no terrain processed(default)
;                 1-terrain/sub_terrain pixels are colored/removed
;                 2-down-slant terrain leak to up beam corrected 
;                 3-1+2
;                 If the input is negative terrain retrieval info is printed 
;       THRESH:  input, float, threshold level factor (must be positive or 0;
;                       0=no thresh,default); data below thresh*std(noise) is 
;                       not displayed; if set to 9 applies the reflectivity
;                       mask based on bit 2(keyword NAV is ignored);
;                       if set to 9X, where X is 0,1,2,3 uses reflectivity mask
;                       without thresholding (X=0) or with mask based on bit X-1  
;       ZOOMIMG: input, float, simple zoom in(>1) and zoom out(<1); default is
;                based on the size of the image
;       ZSRF:    input, when set, adds a strip of the surface reflectivity at 
;                 the bottom of the up-down image (ignored if THRESH is not 9
;                   and the mask bit for surface is not set(data before 2012))
;                 Note: in this case maximum surface reflectivity is limitted 
;                       to +40 dBZ
;                input, integer; maximum surface reflectivity (must be>35dBZ)
;       HELP:    show this text
;
; OUTPUTS:
;
; SIDE EFFECTS:
;       Input files for this routine are selected with the dialog_pickfile
;
;       If STARTT/STOPT is/are defined a single file selection only is allowed.
;       If STARTT and/or STOPT do/does not exist read from/to beginning/end of
;       file.
;
;       Opens one or more graphic windows and/or graphic widgets. Write image
;       file(s).
;
; NOTES:
;
; EXAMPLES:
;
;-
; MODIFICATION HISTORY:
;       Written by:  Samuel Haimov, November 2008
;       S.H., Oct 2009 - fixed several bugs; added SAVEDIR keyword;
;                        added flight level line for up-down images;
;                        added FLTLVL keyword
;       S.H., Nov 2009 - fixed missing altitude decimation when NAV > 1
;                        is given and ALTBEAM is set;  fixed missing
;                        zoom calc. when nav>1; fixed bugs when nav[1]>1;
;                        for thresh=9 now also applies acvcmask
;       S.H., Sep 2010 - added SATMARK keyword
;       S.H., Apr 2011 - added TERRAIN and ASTERDIR keywords
;       S.H., May 2011 - modified handling of SATMARK in lue of the expected
;                        change in the mask bits assignment
;       S.H., Jun 2011 - fixed bugs when applying threshold to velocity prods;
;                        fixed a bug when plotting saturation marks
;       S.H., Jul 2011 - fixed a bug marking terrain for multiple down and/or
;                        multiple down-slant reflectivities
;       S.H., Aug 2011 - added LDR keyword
;       S.H.,     2011 - modified thresh=9 (apply mask) to use the new mask bits
;       S.H., Dec 2011 - fixed a bug in the use of the mask for a single
;                        product file(e.g., one beam with one pwr and dv)
;       S.H., Jan 2012 - added RELIEF keyword
;       S.H., Mar 2012 - added flight level highlights for down and up beams
;                        more than 6 deg off of vertical - this happens
;                        when the aircraft is in steep pitch and/or is turning
;       S.H., Apr 2013 - added ACWIND and ZSRF keywords; fixed a bug for limited
;                        mask data(added line 678 ind3=-1); added keyword BREAKI
;       S.H., May 2013 - added keyword INDIR
;       S.H., Aug 2013 - fixed minor bug about max velocity displayed when the
;                        ground is plotted in black(using +32767 for this)
;       S.H., Oct 2013 - retrieve rx/txpol and rx/txpol2 from nvid instead of
;                        ppmag (ppmag is being replaced with ppreflectivity);
;                        added PPZ keyword; replaced wcr2calthreshnc with
;                        wcrthreshnc -in this case full noise must be available
;                        for NAV>1 to work correctly (this is always the case
;                        starting with ASCII13); modified THRESH keyword to
;                        allow 0 to 3 sigma thresholding when using the mask
;       S.H., Feb 2014 - added spectral width display plots if ppsw is available
;                        (combined up-down only and NO fancy masking is applied)
;       S.H., Mar 2014 - changed how xtop is calculated (use tas) to avoid
;                        its dependence on the aspect
;       S.H., Sep 2014 - changed velocity sign for side beam to match how I have
;                        in L2 quick looks.
;
; Copyright (C) 2008-present, Samuel Haimov, Dept. of Atmos. Sci., University of 
; Wyoming.  This software may be used, copied, or redistributed as long as it 
; is not sold and this copyright notice is reproduced on each copy made. This
; routine is provided as is without any expressed or implied warranties
; whatsoever.  Other limitations apply as described in the Readme file.
;

pro wcr2plotnc,minmax,acwind=acwind,breaki=breaki,indir=indir,$
               altbeam=altbeam,aspect=aspect,asterdir=asterdir,colortb=colortb,$
               fltlvl=fltlvl,leak=leak,mark=mark,nav=nav,noerase=noerase,$
               profile=profile,savedir=savedir,saveimg=saveimg,startt=startt,$
               stopt=stopt,thresh=thresh,zoomimg=zoomimg,satmark=satmark,$
               terrain=terrain,ldr=ldr,relief=relief,zsrf=zsrf,ppz=ppz,help=help

;  on_error,2

  if keyword_set(help) then begin
    doc_library,'wcr2plotnc'
    return
  endif

; Initilize

  resetps2x
  alt    =-999.    ; initialize altitude; -999=no altitude 
  minrefl=-55.     ; limit the minimum reflectivity displayed
  maxrefl= 35.     ; limit the max reflectivity displayed w/o ZSRF
  zmin=minrefl
  zmax=maxrefl
  ldrmin=-45.
  ldrmax=+5.
  zsrfc=-999.      ; will be used to hold surface dBZ (ZSRF set)
  maxvel=15.8      ; limit the min/max velocity displayed
  dvminmax=maxvel
  maxppsw=2.
  minppsw=1.
  polstr=['v','h'] ; polarization chars for pol Id 0 and 1
  dvstatus=['uncorrected','corrected'] ; for dvstatus Id 0 and 1

  if n_elements(asterdir) eq 0 then $
    asterdir='/netdata/kingair_data/aster/' ; '/net/tiger/aster/'
  if n_elements(terrain) eq 0 then terrain=0
  if terrain lt 0 then begin 
     verbose=1
     terrain=-terrain
  endif else verbose=0
  
  if keyword_set(ppz) then pps='pp' else pps=''
  if n_elements(altbeam) eq 0 then altbeam=0
  if n_elements(mark) gt 0 then comm=mark $
  else if n_elements(thresh) eq 0 then comm='QUICK LOOK' else comm=''
  if n_elements(thresh) eq 0 then begin
    sfactor=0 
  endif else begin
    sfactor=thresh>0.
    if sfactor eq 9 then sfactor=93
  endelse

  if n_elements(nav) eq 0 then nav=1 $
  else if sfactor eq 0 and total(nav) gt 1 then begin
    message,'No thresholding. NAV will be ignored !',/info
    nav=1
  endif
  if sfactor ge 90 then begin
    message,'Reflectivity mask is requested. NAV will be ignored !',/info
    message,'Thresholding set to bit'+num2str(sfactor-91,fr=0),/info
    nav=1
    sf=2^((sfactor-90)<3)-1
    sfactor=9  
  endif

  if n_elements(zsrf) eq 1 and sfactor eq 9 then $
    zsrf=5*(zsrf eq 1)+(zsrf-maxrefl)*(zsrf gt maxrefl) $
  else zsrf=0

; Set flags for reflectivity and velocity limits 

  if n_params() gt 0 then begin
    if minmax[0] ne -999 then dvminmax=minmax[0]
    if n_elements(minmax) gt 1 then $
      if minmax[1] ne -999 then zmax=minmax[1]
    if n_elements(minmax) eq 3 then $
      if minmax[2] ne -999 then zmin=minmax[2]
  endif

; Set image colors

  if n_elements(colortb) eq 0 then begin
    loadct,0,get_names=junk1
    if n_elements(junk1) lt 44 then begin
      zcol =33
      dvcol=17
    endif else begin
      zcol =41
      dvcol=51
    endelse
  endif else if n_elements(colortb) gt 1 then begin
    zcol =colortb[0]
    dvcol=colortb[1]
  endif else zcol=(dvcol=colortb[0])

  if dvcol eq 42 then setwdv=[-1,dvcol,1,1,1] else $
  if dvcol eq 43 or dvcol eq 51 then setwdv=[-1,dvcol,1,1,0] $
  else setwdv=[-1,dvcol,0,1,1]

; Select files to display

  if n_elements(indir) eq 0 || indir eq '' then begin
    if n_elements(startt)+n_elements(stopt) gt 0 then junk1=0 else junk1=1
    flst=dialog_pickfile(/read, filter='WCR.*.nc',path=datdir,multi=junk1)
    if flst[0] eq '' then return
  endif else begin
    flst=file_search(indir+'WCR.*.nc')
    if flst[0] eq '' then message,'No WCR nc files found in '+indir
  endelse

; Determine plot parameters 

  if n_elements(saveimg) eq 0 then saveimg=0 $
  else begin
    saveimg[0]=fix(saveimg[0])<3
    if n_elements(saveimg) eq 1 then begin res=200 & ori=1 & endif $
    else if n_elements(saveimg) eq 2 then begin ori=saveimg[1] & res=200
    endif else begin res=saveimg[2] & ori=saveimg[1] & endelse ; ori=1=landscape
    if ori then begin
      xysize  =[10.5,7.0]
      xyoff   =[0.4,10.75]
      xypgpos0=[0.0,1.0]   ; in inches
      xypgpos =[0.2,0.5]   ; in inches
      rot=1
    endif else begin
      xysize  =[7.2,8.5]
      xyoff   =[0.9,0.9]
      xypgpos0=[0.0,1.0]
      xypgpos =[0.5,0.5]
      rot=2
    endelse
  endelse

  nfiles=n_elements(flst)

; Loop over all selected files and display reflectivity and velocity

  for kk=0,nfiles-1 do begin

    wcrfnc=file_basename(flst[kk])
    print,'==================================================================='
    message,'Display '+wcrfnc+', '+num2str(kk+1)+' of '+ $
             num2str(nfiles),/info

  ; Determine ps/pdf/png base file name and create saveimg quicklook directory
  ; (start-stop times match the plot; later other identifiers are appended)

    if saveimg[0] ne 0 then begin
      if n_elements(startt) eq 0 then begin
        junk2=strpos(wcrfnc,'_')                  
        filen=strmid(wcrfnc,0,junk2+11)
      endif else begin
        junk1=time2timestr(startt,divider='')
        junk2=strpos(wcrfnc,'_')                  
        filen=strmid(wcrfnc,0,junk2-6)+junk1+strmid(wcrfnc,junk2,11)
      endelse
      if n_elements(stopt) gt 0 then begin
        junk2=strpos(filen,'_',/reverse_search)
        junk3=time2timestr(stopt,divider='')
        filen=strmid(filen,0,junk2+1)+junk3+strmid(filen,junk2+11)
      endif
      if n_elements(savedir) eq 0 then $
        savedir=file_dirname(flst[kk])+'/quicklook/' $
      else if strmid(savedir,strlen(savedir)-1L) ne path_sep() then $
        savedir=savedir+path_sep()
      file_mkdir,savedir
      filen=savedir+filen
    endif

  ; ---------------------------------------------------------------------------- 
  ; Read the appropriate data from a netcdf file

    ncid    =ncdf_open(flst[kk])
    daqalg  =strmid(ncdf_attload(ncid,'WCRdaqalg',/global),0,2)
    platform=ncdf_attload(ncid,'Platform',/global) ; platform id
    timeint =float(ncdf_attload(ncid,'WCRtimeint',/global)); between profiles
    radrng  =wcr2loadnc(ncid,'range') ; range gates in meters ;modified later
    rgs     =float(ncdf_attload(ncid,'WCRrangesampling',/global))
    nrg     =n_elements(radrng) ; this may change later if nav[1]>1
    beamid  =ncdf_attload(ncid,'WCR_BeamID',/global)
    beamstr =ncdf_attload(ncid,'WCR_BeamName',/global)
    beamstr =strtrim(strsplit(beamstr,',',/extract),2) ; ordered by beamid

    ; read reflectivity, reflectivity mask, and the needed attributes
    mom0=wcr2loadnc(ncid,pps+'reflectivity',startt=startt,stopt=stopt); in mm6/m3
    misval=float(ncdf_attload(ncid,'reflectivity','_FillValue'))
    if size(mom0,/n_dim) eq 3 then nz=(size(mom0,/dim))[2] else nz=1
    zant=ncdf_attload(ncid,'reflectivity','antenna')
    zbeam=long(ncdf_attload(ncid,'reflectivity','beamid'))
    zbeam1=zbeam  ; preserve the original; will use it for LDR plots
    zid=strcompress(ncdf_attload(ncid,'reflectivity','npid'),/remove_all)
    junk1=long(ncdf_attload(ncid,'reflectivity','rxpol'))
    junk2=long(ncdf_attload(ncid,'reflectivity','txpol'))
    zpol=polstr[junk1]+polstr[junk2]
    zmask=fix(wcr2loadnc(ncid,pps+'reflectivity_mask',startt=startt,stopt=stopt))
    zmasksat=8  ; 3rd bit; LATER use zmask 'valuename' attribute
    satval=-999 ; use this value to set all pixels affected by saturation

    ; read velocity, and the needed attributes
    mom1=wcr2loadnc(ncid,'velocity',startt=startt,stopt=stopt)
    if size(mom1,/n_dim) eq 3 then ndv=(size(mom1,/dim))[2] else ndv=1
    dvant=ncdf_attload(ncid,'velocity','antenna')
    dvbeam=long(ncdf_attload(ncid,'velocity','beamid'))
    dvid=strcompress(ncdf_attload(ncid,'velocity','nvid'),/remove_all)
    dvmax=float(ncdf_attload(ncid,'velocity','maxvel'))
    dvstatusid=float(ncdf_attload(ncid,'velocity','statusid'))
    if daqalg eq 'FF' then begin
      junk1=long(ncdf_attload(ncid,'psd','rxpol'))
      junk2=long(ncdf_attload(ncid,'psd','txpol'))
      dvpol2=(dvpol=polstr[junk1]+polstr[junk2])  ; RxTxRxTx polorization
    endif else begin      
      junk1=long(strmid(dvid,1,1))  ; rxpol
      junk2=long(strmid(dvid,2,1))  ; txpol
      dvpol=polstr[junk1]+polstr[junk2]  ; RxTx pol for 1st pulse in the pp
      junk1=long(strmid(dvid,3,1))  ; rxpol2
      junk2=long(strmid(dvid,4,1))  ; txpol2
      dvpol2=polstr[junk1]+polstr[junk2]  ; RxTxRxTx polorization
    endelse

    ; read pulse-pair spectral width var if available
      ppsw=wcr2loadnc(ncid,'ppsw',startt=startt,stopt=stopt) ; will return -1 if does not exist
      if ppsw[0] ne -1 then begin
        ind3=where(ppsw eq misval)
        if ind3[0] ne -1 then ppsw[ind3]=!values.f_nan
      endif
      
    ; read in situ measured wind contribution into the beam
    acwcb=wcr2loadnc(ncid,'acwcbeam',startt=startt,stopt=stopt)
    
    ; read air-relative speed
    tas=wcr2loadnc(ncid,'TAS',startt=startt,stopt=stopt)
    
    ; read time vector (in unix seconds) and convert to sec after midnight
    time  =wcr2loadnc(ncid,'time',startt=startt,stopt=stopt)
    time  =time-time[0]+(time[0] mod 86400d0) ; midnight left continuous
    timeax=1 ; set timeaxis to display in time format hh:mm:ss
    nprof =n_elements(time) ; this may change later if nav[0] > 1

  ; ----------------------------------------------------------------------------
  ; Determine range-to-profile aspect

    xtopid=1
    if n_elements(aspect) eq 1 then begin
      asp=aspect
      message,'Using '+num2str(asp,fr=1)+' aspect ratio',/info
    endif else if dvstatusid eq 0 then begin
      message,'Platform data not available.',/info
      if platform eq 'N2UW' or platform eq 'N130AR' then $
        asp=rgs/timeint/(90.*(platform eq 'N2UW')+115.*(platform eq 'N130AR')) $
      else begin
        message,'and the platform is not N2UW or N130AR.',/info
        asp=1.
        xtopid=0
      endelse
      message,'Use default profile/range aspect ratio: '+num2str(asp,fr=1),/info
    endif else begin
    ; read aspect ratio 
      asp=wcr2loadnc(ncid,'wcraspect',startt=startt,stopt=stopt)
    ; find if the AC was still on the ground (tas < 50 m/s)
      ind1=where(asp eq misval)
      if ind1[0] ne -1 then begin
        message,'Aircraft on the ground or wrong TAS',/info
        asp[ind1]=!values.f_nan
      endif
      asp=mean(asp,/nan)
      message,'Mean aspect based on average TAS: '+num2str(asp,fr=1),/info
    endelse

  ; ----------------------------------------------------------------------------
  ; Load altitude and check for dropped GPS signal (alt is normally GALT
  ; for C130. During VOCALS there were drops in GPS which affected
  ; aircraft gps corrected velocities leaving gaps in the - Chris Webster
  ; is working on a fix to fill the gaps)

    acvcmask=1. ; initialize if the variable is not available
    beamvecid=fix(ncdf_attload(ncid,'wcrbeamvector','beamid'))
    if dvstatusid ne 0 and altbeam ne 0 then begin
      alt=wcr2loadnc(ncid,'ALT',startt=startt,stopt=stopt)
      beamvec=wcr2loadnc(ncid,'wcrbeamvector',startt=startt,stopt=stopt)
      acvcmask=wcr2loadnc(ncid,'acvcmask',startt=startt,stopt=stopt)
      trng='Altitude [km]'
      ind1 =where(alt lt 0.)
      if ind1[0] ne -1 then begin
        message,'Found '+num2str(n_elements(ind1))+$
               ' GPS drops or negative PALT('+ $
                num2str(n_elements(ind1)/float(nprof),fr=1)+ $
                '%) between '+sec2timestr(time[ind1[0]])+$
                ' and '+sec2timestr(time[ind1[n_elements(ind1)-1L]]),/info
      endif
    endif else if altbeam ne 0 then begin
      message,'Altitude var ALT not available. Set ALTBEAM to 0',/info
      trng='Range [km]'
    endif else trng='Range [km]'

  ; ----------------------------------------------------------------------------
  ; Take care of terrain related radar reflectivity pixels if requested

    if terrain gt 0 and dvstatusid ne 0 then begin
      terr=abs(terrain)
      message,'Processing terrain pixels...',/info
      ind4=(where(beamstr eq 'down'))[0]
      ind2=where(zbeam eq beamid[ind4])   ; find down reflectivities
      ind4=(where(strmid(beamstr,0,5) eq 'down-'))[0]
      ind3=where(zbeam eq beamid[ind4])   ; find down-fore/aft reflectivities

      if ind2[0] ne -1 or ind3[0] ne -1 then begin
        if alt[0] eq -999. then begin     ; beamvec(3,nprof,nbeams)
          beamvec=wcr2loadnc(ncid,'wcrbeamvector',startt=startt,stopt=stopt)
          beamvecid=fix(ncdf_attload(ncid,'wcrbeamvector','beamid'))
          alt=wcr2loadnc(ncid,'ALT',startt=startt,stopt=stopt)
        endif
        lat=wcr2loadnc(ncid,'LAT',startt=startt,stopt=stopt)
        lon=wcr2loadnc(ncid,'LON',startt=startt,stopt=stopt)
        junk2=min(abs(lon))-0.2
        if junk2 lt 100. then junk2='0'+num2str(junk2,fr=1) $
        else junk2=num2str(junk2,fr=1)
        junk3=max(abs(lon))+0.2
        if junk3 lt 100. then junk3='0'+num2str(junk3,fr=1) $
        else junk3=num2str(junk3,fr=1)  
        junk1=['N'+num2str(min(lat)-0.2,fr=1),$
               'N'+num2str(max(lat)+0.2,fr=1),$
               'W'+junk2,'W'+junk3]

        message,'Retrieving ASTER DEM data for sector:',/info
        message,'('+strjoin(junk1,',')+')',/info
 
        aster2hll,asterdir,dema,lata,lona,sector=junk1
      endif

      junk2=nrg-1L
      junk3=lindgen(nprof)
      junk5=10.^(0.1*(maxrefl+1.))

      if ind2[0] ne -1 and terr ne 2 then begin ; find down beam terrain
        ind5=where(beamvecid eq (zbeam[ind2])[0])
        srrg=wcr2srr(alt,beamvec[*,*,ind5],dema,lata,lona,lat,lon,radrng,$
                     tol=15.,/gates,verbose=verbose); down beam

        ind6=where(finite(srrg))
        message,'Down beam: found '+num2str(n_elements(ind6))+ $
                 ' out of '+num2str(nprof)+' profiles with surface return.',/info
        if ind6[0] eq -1 then begin
          if verbose then $
            message,'For down beam terrain not found.',/info
        endif else begin
          ; color terrain pixels starting from 2 gates above
          dBZ=mom0[*,*,ind2]     
          if n_elements(ind2) eq 1 then begin     
            dBZ[(srrg[ind6]-2)>0,    junk3[ind6]]=junk5
            dBZ[(srrg[ind6]-1)>0,    junk3[ind6]]=junk5
            dBZ[ srrg[ind6],         junk3[ind6]]=junk5
            dBZ[(srrg[ind6]+1)<junk2,junk3[ind6]]=junk5
            dBZ[(srrg[ind6]+2)<junk2,junk3[ind6]]=junk5
            ; remove sub-terrain pixels; mark them as Infinity
              ; to allow coloring (keyword background)
            junk4=(srrg[ind6]+3)<junk2
            for k=0,n_elements(ind6)-1L do $
              dBZ[junk4[k]:junk2,ind6[k]]=!values.f_infinity
          endif else begin
            for k=0,n_elements(ind2)-1 do begin
              dvt=dBZ[*,*,ind2[k]]
              dvt[(srrg[ind6]-2)>0,    junk3[ind6]]=junk5
              dvt[(srrg[ind6]-1)>0,    junk3[ind6]]=junk5
              dvt[ srrg[ind6],         junk3[ind6]]=junk5
              dvt[(srrg[ind6]+1)<junk2,junk3[ind6]]=junk5
              dvt[(srrg[ind6]+2)<junk2,junk3[ind6]]=junk5
            ; remove sub-terrain pixels; mark them as Infinity
              ; to allow coloring (keyword background)
              junk4=(srrg[ind6]+3)<junk2
              for l=0,n_elements(ind6)-1L do $
                dvt[junk4[l]:junk2,ind6[l]]=!values.f_infinity
              dBZ[*,*,ind2[k]]=dvt            
            endfor
          endelse
          mom0[*,*,ind2]=dBZ
        endelse
      endif

      if ind3[0] ne -1 then begin  ; find down-slant beam terrain
        ind5=where(beamvecid eq (zbeam[ind3])[0])
        srrg=wcr2srr(alt,beamvec[*,*,ind5],dema,lata,lona,lat,lon,radrng,$
                     tol=15.,/gates,verbose=verbose); down-slant beam

        ; remove down slanted beam surface return mean leak to up beam
        ; assume 56 dB mean isolation and no more than 120 m spread in range

        ind6=where(finite(srrg))
        message,'Down-slant beam: found '+num2str(n_elements(ind6))+ $
                 ' out of '+num2str(nprof)+' profiles with surface return.',/info
        if ind6[0] eq -1 then begin
          if verbose then $
            message,'For down-slant beam terrain not found.',/info
        endif else begin 
        ; if requested and up or side beam is active remove the leak
          if terr ge 2 then begin
            ind4=(where(beamstr eq 'up'))[0]   ; find up or side reflectivities
            ind1=where(zbeam eq beamid[ind4]) 
            if ind1[0] eq -1 then ind4=(where(beamstr eq 'side'))[0]
            ind1=where(zbeam eq beamid[ind4]) 
            if ind1[0] ne -1 then begin
              dnsupiso=10.^(5.6) 
            ;  dnsupiso=10.^(4.2) ; 20110701 16:55, C130 ~650m alt; down pwr saturated 
              message,'Remove leak from down-slant to up beam using '+$
                       num2str(10.*alog10(dnsupiso),fr=0)+ ' dB isolation',/info  
              dBZ=mom0[*,*,ind1]
              dvt=mom0[*,*,ind3]
              junk1=round(100./rgs)
              for k=-round(30./rgs),round(90./rgs) do begin
                ind4=srrg+k
                ind5=where(ind4 ge 0 and ind4 lt nrg) ;this also eliminates NaNs
                dBZ[ind4[ind5],junk3[ind5]]=dBZ[ind4[ind5],junk3[ind5]]-$
                                            dvt[ind4[ind5],junk3[ind5]]/dnsupiso
              endfor
              mom0[*,*,ind1]=dBZ
            endif
          endif 

          if terr ne 2 then begin
            ; color terrain pixels
            dBZ=mom0[*,*,ind3]
            if n_elements(ind3) eq 1 then begin     
              dBZ[(srrg[ind6]-2)>0,    junk3[ind6]]=junk5
              dBZ[(srrg[ind6]-1)>0,    junk3[ind6]]=junk5
              dBZ[ srrg[ind6],         junk3[ind6]]=junk5
              dBZ[(srrg[ind6]+1)<junk2,junk3[ind6]]=junk5
              dBZ[(srrg[ind6]+2)<junk2,junk3[ind6]]=junk5
              ; remove sub-terrain pixels; mark them as Infinity
              ; to allow coloring (keyword background)
              junk4=(srrg[ind6]+3)<junk2
              for k=0,n_elements(ind6)-1L do $
                dBZ[junk4[k]:junk2,ind6[k]]=!values.f_infinity
            endif else begin
              for k=0,n_elements(ind3)-1 do begin
                dvt=dBZ[*,*,ind3[k]]
                dvt[(srrg[ind6]-2)>0,    junk3[ind6]]=junk5
                dvt[(srrg[ind6]-1)>0,    junk3[ind6]]=junk5
                dvt[ srrg[ind6],         junk3[ind6]]=junk5
                dvt[(srrg[ind6]+1)<junk2,junk3[ind6]]=junk5
                dvt[(srrg[ind6]+2)<junk2,junk3[ind6]]=junk5
              ; remove sub-terrain pixels; mark them as Infinity
              ; to allow coloring (keyword background)
                junk4=(srrg[ind6]+3)<junk2
                for l=0,n_elements(ind6)-1L do $
                  dvt[junk4[l]:junk2,ind6[l]]=!values.f_infinity
                dBZ[*,*,ind3[k]]=dvt            
              endfor
            endelse
            mom0[*,*,ind3]=dBZ
          endif
        endelse
      endif    
    endif
  ; ---------------------------------------------------------------------------- 
  ; Check the threshold input and threshold the data (dBZ, dvt)
 
    if sfactor eq 0 then begin 
      message,'No mask/detection applied/performed.',/info
      ind1=where(mom0 le 0.)
    ; the commented line assume that np=nv; so it cannot be used for modes
    ; that have different number of pwr and pp products
      if ind1[0] ne -1 then begin
        mom0[ind1]=!values.f_nan
    ;    mom1[ind1]=!values.f_nan
      endif
      dBZ=10.*alog10(mom0)
      dvt=mom1
    endif else if sfactor eq 9 then begin ; use reflectivity mask
      message,'Reflectivity mask applied...',/info
      dBZ=mom0

      if (where(zmask gt 9))[0] eq -1 and $
      (where(zmask gt 1 and zmask lt 9))[0] eq -1 then begin
        message,'Reflectivity mask points with mask=0 only',/info
        ind1=where(zmask eq 0, comp=ind2) ; find no target values
        if ind1[0] ne -1 then dBZ[ind1]=!values.f_nan
      endif else begin

      ; process the new mask (temporary procedure)
        if sf gt 0 then begin       
          ind1=where(zmask lt sf)               ; remove signal below sf
          if ind1[0] ne -1 then dBZ[ind1]=!values.f_nan
        endif else begin
          ind1=where(dBZ le 0.)
          if ind1[0] ne -1 then dBZ[ind1]=!values.f_nan
        endelse
        ind1=where(zmask ge 2048 and zmask lt 2055);remove leak below sf 
        if ind1[0] ne -1 then dBZ[ind1]=!values.f_nan
        ind1=where((zmask AND 256) eq 256)     ; remove surface clutter pixels
        if ind1[0] ne -1 then dBZ[ind1]=!values.f_nan ;dBZ[ind1]/10.
        ind1=where((zmask AND 512) eq 512)     ; highlight surface gates 
        if ind1[0] ne -1 then dBZ[ind1]=1.e10  ; (+100 dB)
        ind1=where((zmask AND 1024) eq 1024)   ; remove sub-surface gates
        if ind1[0] ne -1 then dBZ[ind1]=!values.f_infinity
        if n_elements(leak) eq 1 then begin
          message,'Remove leak pixels if any.',/info
          ind1=where((zmask AND 2048) eq 2048)            ; remove all cross-talk leak 
          if ind1[0] ne -1 then dBZ[ind1]=!values.f_nan
        endif

      ; to improve visualization of the surface in the quick looks
      ; thicken them by using 4 more subsurface gates

        ind1=where(dBZ eq 1.e10)   ; lonarr(m), m=n_elements(ind1)
        if ind1[0] ne -1 then begin
          ind2=array_indices(dBZ,ind1) ; (3,m), where 3=(rgind,profind,npind)
          if n_elements(ind2[*,0]) lt 3 then ind3=fltarr(n_elements(ind2[0,*])) $
          else ind3=ind2[2,*]
          for k=1,4 do dBZ[(ind2[0,*]+k)<(nrg-1),ind2[1,*],ind3]=1.e10
        endif

      endelse
 
      dBZ=10.*alog10(dBZ)
      dvt=mom1
   
    ; for C130 remove dvt missing ins data points   
      if dvstatusid ne 0 then begin
        ind1=where(acvcmask eq 0.)
        if ind1[0] ne -1 then begin
          message,'There are missing values in ACVCBEAM !!!!!',/info
          message,'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',/info
          if acvcmask[0] ne -1 then dvt[*,ind1,*]=!values.f_nan
        endif
      endif

    ; figure which mask to use with which velocity product
      ; I'll mask velocity products that have matching reflectivities only

      for ll=0,ndv-1L do begin
        ind2=where(zid eq dvid[ll]) ; find matching reflectivity product(s)
        if ind2[0] ne -1 then begin
          if n_elements(ind2) gt 1 then begin
            message,'For DV('+num2str(ll)+') ID '+dvid[ll]+$
             ' velocity there are '+num2str(n_elements(ind2))+$
             ' matching reflectivity products',/info
            message,'Will use first match. I know, there is an issue here',/info
            ind2=ind2[0]
          endif

          if (where(zmask gt 9))[0] eq -1 and $
          (where(zmask gt 1 and zmask lt 9))[0] eq -1 then begin
          ;  message,'Velocity mask points with mask=0 only',/info
            ind1=where(zmask[*,*,ind2] eq 0)
            ind3=-1
          endif else begin
          ; remove all data points based on masked dBZ
            ind1=where(finite(dBZ[*,*,ind2]) eq 0)
            ind3=where(finite(dBZ[*,*,ind2],/infinity) eq 1)
          endelse

          if ind1[0] ne -1 then begin
            junk3=dvt[*,*,ll]
            junk3[ind1]=!values.f_nan
            if ind3[0] ne -1 then junk3[ind3]=!values.f_infinity
            dvt[*,*,ll]=junk3
          endif

        endif else begin
          message,'For DV('+num2str(ll)+') ID '+dvid[ll]+$
            ' velocity there is NO '+$
            ' matching reflectivity product. Mask not applied.',/info
        endelse
      endfor

    endif else begin ; use detection/thresholding procedure

      message,'Reflectivity thresholding applied...',/info
      if n_elements(nav) gt 0 then begin ; Calculate nprof and nrg if nav
        nprof=nprof/nav[0]
        asp=asp/float(nav[0]) ; correct aspect ratio accordingly
        if n_elements(nav) gt 1 then begin
           nrg=nrg/nav[1]
           asp=asp*float(nav[1]) ; correct aspect ratio accordingly
        endif
        message,'NAV Corrected aspect: '+num2str(asp,fr=1),/info
      endif

      dBZ=fltarr(nrg,nprof,nz)  ; initilize thresholded reflectivity
      dvt=fltarr(nrg,nprof,ndv) ; initilize thresholded velocity
      if nav[0] gt 1 then begin
        noise=wcr2loadnc(ncid,pps+'noise',startt=startt,stopt=stopt)
        noisemn=wcr2loadnc(ncid,pps+'noisemn',startt=startt,stopt=stopt)
      endif else $
        noise=wcr2loadnc(ncid,pps+'noisestd',startt=startt,stopt=stopt)
 
      ; I'll threshold velocity products that have matching reflectivities only

      for ll=0,nz-1 do begin
        junk3=num2str(zbeam[ll])
        ind2=where(dvid eq zid[ll]) ; find matching velocity product(s)
        if ind2[0] eq -1 then begin  ; Z w/o matching velocity product
          message,'Thresholding '+zpol[ll]+' reflectivity for beam '+junk3,/info
;          dBZ[*,*,ll]=wcr2calthreshnc(mom0[*,*,ll],noise[*,ll],radrng,$
;                      sfactor,leak=leak,/nan,nav=nav,rradrng=rradrng)
          if nav[0] gt 1 then $ 
            dBZ[*,*,ll]=wcr2threshnc(mom0[*,*,ll],noise[*,ll],radrng,$
                                     mnoise=noisemn[*,ll],sfactor,/nan,$
                                     nav=nav,rradrng=rradrng) $
          else dBZ[*,*,ll]=wcr2threshnc(mom0[*,*,ll],noise[*,ll],radrng,$
                          sfactor,/nan,nav=nav,rradrng=rradrng)
        endif else begin ; Z w/ matching velocity product(s)
        ; if more than one velocity products apply the threshold to all
          message,'Thresholding '+zpol[ll]+' reflectivity and velocity(ies)'+$
                  ' for beam '+junk3,/info
;          dBZ[*,*,ll]=wcr2calthreshnc(mom0[*,*,ll],noise[*,ll],radrng,$
;                      sfactor,mom1[*,*,ind2[0]],leak=leak,/nan,nav=nav,$
;                      rradrng=rradrng,dvt=junk4)
          if nav[0] gt 1 then $ 
            dBZ[*,*,ll]=wcr2threshnc(mom0[*,*,ll],noise[*,ll],radrng,$
                  sfactor,dvc=mom1[*,*,ind2[0]],dvmax=dvmax[ind2[0]],$
                 /nan,nav=nav,rradrng=rradrng,dvt=junk4,mnoise=noisemn[*,ll]) $
          else dBZ[*,*,ll]=wcr2threshnc(mom0[*,*,ll],noise[*,ll],radrng,$
                      sfactor,dvc=mom1[*,*,ind2[0]],dvmax=dvmax[ind2[0]],$
                      /nan,nav=nav,rradrng=rradrng,dvt=junk4)
          dvt[*,*,ind2[0]]=junk4
          if n_elements(ind2) gt 1 and nav[0] le 1 then begin
            ind1=where(finite(dvt[*,*,ind2[0]]) eq 0)
            if ind1[0] ne -1 then begin
               for mm=1,n_elements(ind2)-1L do begin
                 junk4=mom1[*,*,ind2[mm]]
                 junk4[ind1]=!values.f_nan
                 dvt[*,*,ind2[mm]]=junk4
               endfor  
            endif
          endif
        endelse
      endfor
    endelse
    
    nprof = n_elements(dBZ[0,*,0])
    nrg = n_elements(dBZ[*,0,0])
    if nrg lt n_elements(radrng) then radrng=rradrng

  ; Figure out if image zoom is needed; if asp is less than 1 (many profiles 
  ; were averaged) expand range(y-axis)

    if n_elements(zoomimg) eq 1 then $
      imgz=float(abs(zoomimg)) $
    else if asp le 0.15 then imgz=6 $
    else if asp lt 0.30 then imgz=3 $
    else if asp lt 0.50 then imgz=2 $
    else if nprof/asp lt 3.*(n_elements(range)>50) then $
      imgz=2.+(nprof/asp le 100.) else imgz=1. 
    message,'Display image zoom set to: '+num2str(imgz,fr=1),/info
    if imgz ge 0.95 and imgz le 1.05 then imgz=1.

  ; Adjust altitude for NAV > 1 if needed

    if nav[0] gt 1 and altbeam gt 0 and alt[0] ne -999 then begin 
      alt0=alt ; save original alt for non-matching velocity products
    ; interpolate alt to match the decimated number of profiles
      alt=congrid(alt,n_elements(dBZ[0,*,0]),/interp)
    endif 

  ; ----------------------------------------------------------------------------- 
  ; Calculate flown air-distance margins to be displayed

    if keyword_set(profile) or xtopid eq 0 then begin
      xtop=[0L,nprof-1]
      xtt ='Profile'
    endif else begin
      xtop=[0,timeint*mean(tas)*nprof/1000.]
      xtt ='Air Relative Distance [km]'
    endelse

  ; Check if number of profiles to be displayed exceed 65k for pdf conversion

    message,'Displaying/saving image(s); '+num2str(nprof)+' profiles',/info

    if nprof gt 65536 and abs(saveimg[0]) eq 3 then begin
      message,'=========================================================',/info
      message,'Too many profiles. Cannot convert to pdf using Ghostscript',/info
      message,'For '+wcrfnc+' - will skip writing pdf file(s).',/info
      message,'=========================================================',/info
      pdfflg=0
    endif else pdfflg=1

 ; Initilize widget ID vector

    widgetid=lonarr(nz+ndv)-1L

 ; Prepare range axes for up-down combined beam(s)

   rgs=radrng[1]-radrng[0]   ; rangegate sampling
   if altbeam eq 0 or keyword_set(acwind) then begin
     blindrg=fix(2.*radrng[0]/rgs) ; # of gates between up and down with no data
     if (blindrg*rgs-2.*radrng[0]) lt -rgs/2 then blindrg=blindrg+1 $
     else if (blindrg*rgs-2.*radrng[0]) gt rgs/2 then blindrg=blindrg-1
     junk1=(rgs-radrng[0])*(blindrg gt 1)+findgen(blindrg)*rgs
   endif else begin blindrg=1 & junk1=1 & end
   udrng  = [-reverse(radrng),junk1,radrng]
 ; zrg=(where(udrng ge 0))[0]
 ; if udrng[zrg] gt abs(udrng[zrg-1]) then zrg=zrg-1

  ; ----------------------------------------------------------------------------
  ; If SATMARK is set and eligible determine satflag

    satflag=0  ; no saturation or no use of SATMARK
    if keyword_set(SATMARK) and nav[0] eq 1 and imgz eq 1 then begin
      ind1=where((zmask AND zmasksat) eq zmasksat)
      if ind1[0] ne -1 then begin
        satflag=1
        zmask[ind1]=satval
        message,'Pixels with saturation assigned -999.',/info
      endif  
    endif
 ; ----------------------------------------------------------------------------- 
 ; Display/save up-down beam pair images

   ; Process reflectivities

   ind3=(where(beamstr eq 'up'))[0]
   ind1=where(zbeam eq beamid[ind3])     ; find up beam reflectivities
   ind3=(where(beamstr eq 'down'))[0]
   ind2=where(zbeam eq beamid[ind3])     ; find down beam reflectivities

   if ind1[0] ne -1 and ind2[0] ne -1 then begin 

   ; before combining the avialble up-down pairs sort them by co-/cross-pol

     ind3=where(strmid(zpol[ind1],0,1) eq strmid(zpol[ind1],1,1),comp=ind4)
     if ind3[0] ne -1 and ind4[0] ne -1 then ind1=[ind1[ind3],ind1[ind4]]
     ind3=where(strmid(zpol[ind2],0,1) eq strmid(zpol[ind2],1,1),comp=ind4)
     if ind3[0] ne -1 and ind4[0] ne -1 then ind2=[ind2[ind3],ind2[ind4]]

   ; now combine the co-pol and then cross-pol for up and down
   ; (co-cross sorting will not work perfectly if there are different number
   ;  of co- and cross-pol products for the up and the down beams)

     nupdown=min([n_elements(ind1),n_elements(ind2)])
     if ppsw[0] ne -1 then $
      widgetid=lonarr(nz+ndv+nupdown)-1L
       
     for ll=0,nupdown-1L do begin

     ; if ZSRF is requested extract surface reflectivity
       if keyword_set(zsrf) then begin
         ind3=where(zmask[*,*,ind2[ll]] eq 512)
         if ind3[0] ne -1 then begin
           zsrfc=make_array(nprof,value=!values.f_nan)
           ind4=ind3 mod nrg ; find surface range gates
           ind3=ind3/nrg     ; find profiles with srfc return
           zsrfc[ind3]=10.*alog10((mom0[*,*,ind2[ll]])[ind4,ind3])
         endif else zsrfc=-999.      
       endif 

     ; prepare the combined down-up reflectivity data
       mom0=[rotate(dBZ[*,*,ind2[ll]],5),$
             make_array(blindrg,nprof,value=!values.f_nan),$
             dBZ[*,*,ind1[ll]]]
       if satflag then begin  ; prepare saturation mask for the up-down
         mask=[rotate(zmask[*,*,ind2[ll]],5),$
               make_array(blindrg,nprof,value=0),$
               zmask[*,*,ind1[ll]]]
         ind3=where(mask ne satval)
         if ind3[0] eq -1 then message, 'sat ind not found - cannot happen' $
         else mask[ind3]=0 
       endif

     ; and spectral width if available
       if ppsw[0] ne -1 then begin
       ; check if if I have to rebin it 
         if nprof ne n_elements(ppsw[0,*,0]) then begin
           ind3=n_elements(ppsw[0,*,0])
           ind3=ind3/nprof*nprof-1L
           ppsw=rebin(ppsw[*,0:ind3,*],size(dBZ,/dim))
         endif

         mom2=[rotate(ppsw[*,*,ind2[ll]],5),$
             make_array(blindrg,nprof,value=!values.f_nan),$
             ppsw[*,*,ind1[ll]]]
       endif 

     ; interpolate to vertical plane if requested
       if altbeam eq 1 and alt[0] ne -999. then begin 
       ; interpolate image to flight altitude
         mom0=wcrud2altnc(mom0,udrng,alt,altrng=altrng,relief=relief,misval=!values.f_nan)
         if ppsw[0] ne -1 then begin
           mom2=wcrud2altnc(mom2,udrng,alt,altrng=altppsw,relief=relief,misval=!values.f_nan)
           altppsw=altppsw/1000.
         endif 
         if satflag then begin
           mask=wcrud2altnc(float(mask),udrng,alt,altrng=altrngs)
           altrngs=altrngs/1000. ; change altitude to km
         endif
         altrng=altrng/1000.   ; change altitude to km
       endif else if altbeam eq 2 and alt[0] ne -999. then begin 
       ; interpolate image to flight altitude and strictly vertical plane
         ind3=(where(beamvecid eq zbeam[ind1[ll]]))[0]
         ind4=(where(beamvecid eq zbeam[ind2[ll]]))[0]
         mom0=wcrud2altnc(mom0,udrng,alt,beamvec[*,*,ind4],beamvec[*,*,ind3],$
                          altrng=altrng,relief=relief,misval=!values.f_nan)
         if ppsw[0] ne -1 then begin
           mom2=wcrud2altnc(mom2,udrng,alt,beamvec[*,*,ind4],beamvec[*,*,ind3],$
                          altrng=altppsw,relief=relief,misval=!values.f_nan)
           altppsw=altppsw/1000.
         endif 
         if satflag then begin
           mask=wcrud2altnc(float(mask),udrng,alt,beamvec[*,*,ind4],$
                            beamvec[*,*,ind3],altrng=altrngs)
           altrngs=altrngs/1000. ; change altitude to km
         endif
         altrng=altrng/1000. ; change altitude to km     
       endif else begin
         altrng=udrng/1000. ; change range to km for plotting
         altrngs=(altppsw=(altrng))
       endelse

     ; add a strip (3 gates) of surface reflectivity  
       if zsrfc[0] ne -999. then begin 
         mom0=[fltarr(8,nprof),mom0]
         mom0[0:1,*]=!values.f_nan
         mom0[2:5,*]=replicarr(zsrfc,4,/dim)
         mom0[6:7,*]=!values.f_nan
         altrng=[altrng[0:7]-7*(altrng[1]-altrng[0]),altrng]
       endif

     ; determine reflectivity min/max
       if zmin eq minrefl then junk2=min(mom0,/nan)>minrefl $
       else junk2=zmin>minrefl
       if zmax eq maxrefl then junk3=max(mom0,/nan)<(maxrefl+zsrf) $
       else junk3=zmax<(maxrefl+zsrf)
       if junk2 eq junk3 then begin
         junk2=0.9*junk2
         junk3=1.1*junk3
       endif

     ; if there are no valid pixels to display/save print a message and skip

       if total(finite(mom0)) eq 0 then begin
         message,'Up-Down reflectivity image:',/info
         message,'No valid/signal pixels detected. Skip the image.',/info
       endif else begin
       ; display and/or save the reflectivity plot
         tit='Up|Down '+pps+'Reflectivity ('+zpol[ind1[ll]]+'|'+$
             zpol[ind2[ll]]+', IDs: '+zid[ind1[ll]]+' | '+zid[ind2[ll]]+')'
         if saveimg[0] ge 0 then begin ; display the reflectivity image plot
           aximage,mom0,time,altrng,/tr,/bar,timeax=timeax,setw=[-1,zcol,0,1], $
             tit=wcrfnc+', '+tit,xtit=['UTC','dBZ',xtt],$
             ytit=trng,asp=asp,topx=xtop,imgz=imgz,widgetid=junk1,/silent, $
             min=junk2,max=junk3+(junk2 eq junk3)
           widgetid[ll]=junk1

           if ppsw[0] ne -1 then $
           aximage,mom2,time,altppsw,/tr,/bar,timeax=timeax,setw=[-1,33,0,1], $
             tit=wcrfnc+', Up|Down '+pps+' Spectral Width ('+zpol[ind1[ll]]+'|'+$
                 zpol[ind2[ll]]+', IDs: '+zid[ind1[ll]]+' | '+zid[ind2[ll]]+')',$
             xtit=['UTC','m/s',xtt],ytit=trng,asp=asp,topx=xtop,imgz=imgz,$
             widgetid=junk1,/silent,max=max(mom2,/nan)<maxppsw,min=minppsw
           widgetid[nupdown+ll]=junk1

           if altbeam gt 0 and alt[0] ne -999. and keyword_set(fltlvl) then $
           begin
             if imgz ne 1. then begin
               junk4=congrid(time,imgz*n_elements(mom1[0,*,0]),/interp)
               junk5=congrid(alt,imgz*n_elements(mom1[0,*,0]),/interp)/1000.         
               oplot,junk4,junk5,thick=2,lin=1

             ; overplot flight segments with down beam off vertical > 10 deg
               junk1=reform(90.+atan(beamvec[2,*,ind4]/ $
                     (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2)))*!radeg)
               junk1=congrid(reform(junk1),imgz*n_elements(mom1[0,*,0]),/interp)
               ind3=where(junk1 gt 10.)
               if ind3[0] ne -1 then $
                 oplot,junk4[ind3],junk5[ind3],symsize=0.5,psym=6
             endif else begin
               oplot,time,alt/1000.,thick=2,lin=1

             ; overplot flight segments with down beam off vertical > 10 deg
               junk1=reform(90.+atan(beamvec[2,*,ind4]/ $
                     (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2)))*!radeg)
               ind3=where(junk1 gt 10.)
               if ind3[0] ne -1 then $
                 oplot,time[ind3],alt[ind3]/1000.,symsize=0.5,psym=6
             endelse
           endif
           if satflag then begin
             ind3=where(mask lt 0.5*satval)
             if ind3[0] ne -1 then begin 
               ind5=array_indices(mask,ind3)
               oplot,time[ind5[1,*]],altrngs[ind5[0,*]],psym=4,thick=3
             endif
           endif
           if keyword_set(breaki) then stop          
         endif
         if saveimg[0] ne 0 then begin ; make ps and save in the requested format
           junk4=filen+'.dBZ'+zpol[ind2[ll]]+num2str(ll+1)+'.updown'
           setxpswin,99,1,landscape=ori,xysize=xysize,xyoff=xyoff,$
                     psfile=junk4+'.ps',/sil
           message,'Making ps file: '+junk4+'.ps',/info
           loadct1,zcol,rev=0,bf=2,/silent & !p.color=255B

           aximage,mom0,time,altrng,/tr,/bar,timeax=timeax, $
             tit=tit,xtit=['UTC','dBZ',xtt],$
             ytit=trng,topx=xtop,/silent,min=junk2,max=junk3+(junk2 eq junk3)

           if altbeam gt 0 and alt[0] ne -999. and keyword_set(fltlvl) then $
           begin $
             oplot,time,alt/1000.,thick=2,lin=1
             ; overplot flight segments with down beam off vertical > 10 deg
               junk1=reform(90.+atan(beamvec[2,*,ind4]/ $
                     (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2)))*!radeg)
             ind3=where(junk1 gt 10.)
             if ind3[0] ne -1 then $
               oplot,time[ind3],alt[ind3]/1000.,symsize=0.5,psym=6
           endif
           if satflag then begin
             ind3=where(mask lt 0.5*satval)
             if ind3[0] ne -1 then begin
               ind5=array_indices(mask,ind3)
               oplot,time[ind5[1,*]],altrngs[ind5[0,*]],psym=4,thick=1
             endif
           endif

           xyouts,xypgpos0[0],xypgpos0[1],comm,charsize=1.2,/normal,/noclip
           datestamp,comm=' '+file_basename(junk4),xypgpos=xypgpos,/inches,$
                     chars=0.9,rot=rot
           resetps2x

           if abs(saveimg[0]) eq 2 then begin
             message,'Converting to png: '+junk4+'.png',/info
             spawn,'idlps2png -r'+num2str(res)+' '+junk4+'.ps '+junk4+'.png'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'   
           endif else if pdfflg and abs(saveimg[0]) eq 3 then begin
             message,'Converting to pdf: '+junk4+'.pdf',/info
             spawn,'idlps2pdf '+junk4+'.ps '+junk4+'.pdf'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'     
           endif else if pdfflg eq 0 and abs(saveimg[0]) eq 3 then begin
             message,'Converting to png: '+junk4+'.png',/info
             spawn,'idlps2png -r'+num2str(res)+' '+junk4+'.ps '+junk4+'.png'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'   
           endif
         endif
         zbeam[ind1[ll]]=(zbeam[ind2[ll]]=-9)  ; use -9 to mark plotted beams 
       endelse
     endfor ; end loop up-down combined reflectivities; ll
   endif    ; end up-down combined reflectivities

   widgetcount=(where(widgetid eq -1))[0]

   ;  Process velocities

   ind3=(where(beamstr eq 'up'))[0]
   ind1=where(dvbeam eq beamid[ind3])     ; find up beam velocities
   ind3=(where(beamstr eq 'down'))[0]
   ind2=where(dvbeam eq beamid[ind3])     ; find down beam velocities

   if ind1[0] ne -1 and ind2[0] ne -1 then begin   

   ; before combining the avialble up-down pairs sort them by co-/cross-pol

     ppol=dvpol+dvpol2
     ind3=where(strmid(ppol[ind1],0,1) eq strmid(ppol[ind1],1,1) and $
                strmid(ppol[ind1],2,1) eq strmid(ppol[ind1],3,1) and $
                strmid(ppol[ind1],0,1) eq strmid(ppol[ind1],3,1),comp=ind4)
     if ind3[0] ne -1 and ind4[0] ne -1 then ind1=[ind1[ind3],ind1[ind4]]
     ind3=where(strmid(ppol[ind2],0,1) eq strmid(ppol[ind2],1,1) and $
                strmid(ppol[ind2],2,1) eq strmid(ppol[ind2],3,1) and $
                strmid(ppol[ind2],0,1) eq strmid(ppol[ind2],3,1),comp=ind4)
     if ind3[0] ne -1 and ind4[0] ne -1 then ind2=[ind2[ind3],ind2[ind4]]

   ; now combine the co-pol and then cross-pol for up and down
   ; (co-cross sorting will not work perfectly if there are different number
   ; of co- and cross-pol products for the up and the down beams; it is also
   ; more complicated with pulse pairs and I am not taking care of hhvv
   ; combination - two different co-pols are used for the pulse-pair)

     for ll=0,min([n_elements(ind1),n_elements(ind2)])-1L do begin
     ; prepare the combined down-up velocity data
     ; if AC wind requested add to the image else generate NaNs for blindrg
       if keyword_set(acwind) and dvstatusid ne 0 then begin  
          ind3=(blindrg-1)/2
          ind4=(blindrg-ind3-4)>1
          ind5=(where(beamvecid eq dvbeam[ind2[ll]]))[0]
          ind6=(where(beamvecid eq dvbeam[ind1[ll]]))[0]
          junk1=[replicate(!values.f_nan,2,nprof),$
           replicarr(-acwcb[*,ind5],ind3,/dim),$
           replicarr(acwcb[*,ind6],ind4,/dim),$
           replicate(!values.f_nan,2,nprof)]
       endif else junk1=make_array(blindrg,nprof,value=!values.f_nan)

       mom1=[rotate(dvt[*,*,ind2[ll]],5),junk1,-dvt[*,*,ind1[ll]]]

     ; make all NaNs positive (I use -NaNs to identify sub-surface pixels)
       ind3=where(finite(mom1) eq 0 and finite(mom1,/infinity) eq 0)
       mom1[ind3]=!values.f_nan

     ; interpolate to vertical plane if requested
       if altbeam eq 1 and alt[0] ne -999. then begin
         if nav[0] gt 1 and n_elements(mom1[0,*,0]) ne nprof then junk1=alt0 $
         else junk1=alt 
       ; interpolate image to flight altitude
         mom1=wcrud2altnc(mom1,udrng,junk1,altrng=altrng,relief=relief,misval=!values.f_nan)
         altrng=altrng/1000. ; change altitude to km
         ind4=(where(beamvecid eq dvbeam[ind2[ll]]))[0]
       endif else if altbeam eq 2 and alt[0] ne -999. then begin 
         if nav[0] gt 1 and n_elements(mom1[0,*,0]) ne nprof then junk1=alt0 $
         else junk1=alt 
       ; interpolate image to flight altitude and strictly vertical plane
         ind3=(where(beamvecid eq dvbeam[ind1[ll]]))[0]
         ind4=(where(beamvecid eq dvbeam[ind2[ll]]))[0]

         mom1=wcrud2altnc(mom1,udrng,junk1,beamvec[*,*,ind4],beamvec[*,*,ind3],$
                          altrng=altrng,relief=relief,misval=!values.f_nan)
         altrng=altrng/1000. ; change altitude to km     
       endif else altrng=udrng/1000. ; change range to km for plotting

     ; determine velocity minmax
       if dvminmax eq maxvel then begin
         ind3=where(mom1 ne 32767)
         junk2=max(abs(mom1[ind3]),/nan)<maxvel
       endif else junk2=dvminmax<maxvel

       if total(finite(mom1)) eq 0 then begin
         message,'Up-Down velocity image:',/info
         message,'No valid/signal pixels detected. Skip the image.',/info
       endif else begin
       ; display and/or save the velocity plot
         if saveimg[0] ge 0 then begin ; display the velocity image plot
           aximage,mom1,time,altrng,/tr,/bar,timeax=timeax,setw=setwdv, $
             tit='Up/Down '+dvstatus[dvstatusid]+' velocity (Nyquist:'+ $
             num2str(dvmax[ind1[ll]],fr=1)+'/'+num2str(dvmax[ind2[ll]],fr=1)+$
             ', IDs: '+dvid[ind1[ll]]+'/'+dvid[ind2[ll]]+')',$
             ytit=trng,asp=asp,topx=xtop,xtit=['UTC','m/s, positive is Up',xtt],$
             imgz=imgz,widgetid=junk1,/silent,min=-junk2,max=junk2+(junk2 eq 0.)
           widgetid[widgetcount+ll]=junk1
           if altbeam gt 0 and alt[0] ne -999. and keyword_set(fltlvl) then $
           begin
             if imgz ne 1 then begin          
               junk4=congrid(time,imgz*n_elements(mom1[0,*,0]),/interp)
               junk5=congrid(alt,imgz*n_elements(mom1[0,*,0]),/interp)/1000.         
               oplot,junk4,junk5,thick=2,lin=1

             ; overplot flight segments with down beam off vertical > 10 deg
               junk1=reform(90.+atan(beamvec[2,*,ind4]/ $
                     (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2)))*!radeg)
               junk1=congrid(reform(junk1),imgz*n_elements(mom1[0,*,0]),/interp)
               ind3=where(junk1 gt 10.)
               if ind3[0] ne -1 then $
                 oplot,junk4[ind3],junk5[ind3],symsize=0.5,psym=6
             endif else if nav[0] gt 1 and n_elements(mom1[0,*,0]) ne nprof then $
             begin
               oplot,time,alt0/1000.,thick=2,lin=1
             ; overplot flight segments with down beam off vertical > 10 deg
               junk1=reform(90.+atan(beamvec[2,*,ind4]/ $
                     (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2)))*!radeg)
               ind3=where(junk1 gt 10.)
               if ind3[0] ne -1 then $
                 oplot,time[ind3],alt0[ind3]/1000.,symsize=0.5,psym=6
             endif else begin
               oplot,time,alt/1000.,thick=2,lin=1
             ; overplot flight segments with down beam off vertical > 10 deg
               junk1=reform(90.+atan(beamvec[2,*,ind4]/ $
                     (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2)))*!radeg)
               ind3=where(junk1 gt 10.)
               if ind3[0] ne -1 then $
                 oplot,time[ind3],alt[ind3]/1000.,symsize=0.5,psym=6
             endelse
           endif
           if keyword_set(breaki) then stop          
         endif
         if saveimg[0] ne 0 then begin ; make ps and save in the requested format
           junk4=filen+'.DV'+num2str(ll+1)+'.updown'
           setxpswin,99,1,landscape=ori,xysize=xysize,xyoff=xyoff,$
                     psfile=junk4+'.ps',/sil
           message,'Making ps file: '+junk4+'.ps',/info
           loadct1,setwdv[1],rev=setwdv[2],bf=2,ccol=setwdv[4],/silent
           !p.color=255B
           aximage,mom1,time,altrng,/tr,/bar,timeax=timeax, $
             tit='Up/Down '+dvstatus[dvstatusid]+' velocity (Nyquist:'+ $
             num2str(dvmax[ind1[ll]],fr=1)+'/'+num2str(dvmax[ind2[ll]],fr=1)+$
             ', IDs: '+dvid[ind1[ll]]+'/'+dvid[ind2[ll]]+')',$
             ytit=trng,topx=xtop,/silent,xtit=['UTC','m/s, positive is Up',xtt],$
             min=-junk2,max=junk2+(junk2 eq 0.)
           if altbeam gt 0 and alt[0] ne -999. and keyword_set(fltlvl) then $
           begin $
             if nav[0] gt 1 and n_elements(mom1[0,*,0]) ne nprof then begin
               oplot,time,alt0/1000.,thick=2,lin=1
             ; overplot flight segments with down beam off vertical > 10 deg
               junk1=reform(90.+atan(beamvec[2,*,ind4]/ $
                     (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2)))*!radeg)
               ind3=where(junk1 gt 10.)
               if ind3[0] ne -1 then $
                 oplot,time[ind3],alt0[ind3]/1000.,symsize=0.5,psym=6
             endif else begin
               oplot,time,alt/1000.,thick=2,lin=1
             ; overplot flight segments with down beam off vertical > 10 deg
               junk1=reform(90.+atan(beamvec[2,*,ind4]/ $
                     (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2)))*!radeg)
               ind3=where(junk1 gt 10.)
               if ind3[0] ne -1 then $
                 oplot,time[ind3],alt[ind3]/1000.,symsize=0.5,psym=6
             endelse
           endif

           xyouts,xypgpos0[0],xypgpos0[1],comm,charsize=1.2,/normal,/noclip
           datestamp,comm=' '+file_basename(junk4),xypgpos=xypgpos,/inches,$
                     chars=0.9,rot=rot
           resetps2x

           if abs(saveimg[0]) eq 2 then begin
             message,'Converting to png: '+junk4+'.png',/info
             spawn,'idlps2png -r'+num2str(res)+' '+junk4+'.ps '+junk4+'.png'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'   
           endif else if pdfflg and abs(saveimg[0]) eq 3 then begin
             message,'Converting to pdf: '+junk4+'.pdf',/info
             spawn,'idlps2pdf '+junk4+'.ps '+junk4+'.pdf'
             message,'Removing '+junk4+'.ps',/info
            file_delete,junk4+'.ps'     
           endif else if pdfflg eq 0 and abs(saveimg[0]) eq 3 then begin
             message,'Converting to png: '+junk4+'.png',/info
             spawn,'idlps2png -r'+num2str(res)+' '+junk4+'.ps '+junk4+'.png'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'   
           endif
         endif
         dvbeam[ind1[ll]]=(dvbeam[ind2[ll]]=-9)  ; use -9 to mark plotted beams 
       endelse
     endfor ; end loop up-down velocities; ll
   endif    ; end up-down combined velocities

   widgetcount=(where(widgetid eq -1))[0]

 ; ----------------------------------------------------------------------------- 
 ; Display/save the rest of reflectivity and velocity images 

   ; Process the rest of the reflectivities

   ind1=where(zbeam ne -9)
   if ind1[0] ne -1 then begin
     for ll=0,n_elements(ind1)-1L do begin
       rrng=radrng
       mom0=dBZ[*,*,ind1[ll]]

     ; Process cross-pol products as LDR if requested; it is assumed that
     ; there is up-down dual-pol antenna combination - so this is not done
     ; for the up-down quick looks and cross-pol for up is present here

       ldrstr=''
       if keyword_set(LDR) then $ 
       ; for cross-pol reflectivity calculate and plot LDR
         if zpol[ind1[ll]] eq 'vh' or zpol[ind1[ll]] eq 'hv' then begin

         ; search for a co-pol signal if present calculate ldr
           if zpol[ind1[ll]] eq 'vh' then ind3=where(zpol eq 'hh') $
           else if zpol[ind1[ll]] eq 'hv' then  ind3=where(zpol eq 'vv')  

           if ind3[0] ne -1 then begin ; use the first matching co-pol
             ind2=(where(zbeam[ind1[ll]] eq zbeam1[ind3]))[0]
             if ind2 ne -1 then begin
               mom0=mom0-dBZ[*,*,ind3[ind2]]
               ldrstr=' LDR, ('+zpol[ind1[ll]]+'/'+zpol[ind3[ind2]]+$
                      ', ID:'+num2str(zid[ind1[ll]])+$
                      '/'+num2str(zid[ind3[ind2]])+')'
             endif
           endif
         endif

       if satflag then begin
         mask=zmask[*,*,ind1[ll]]
         ind3=where(mask ne satval)
         if ind3[0] ne -1 then mask[ind3]=0 
       endif
       ind2=(where(zbeam[ind1[ll]] eq beamid))[0]
       if strmid(beamstr[ind2],0,4) eq 'down' or $
          strmid(beamstr[ind2],0,2) eq 'up' then begin    
         altflg=1
       endif else begin
         altflg=0 ; do not interpolate for altitude if side beams
       endelse

       if alt[0] ne -999. and altbeam gt 0 then $
         ind4=(where(beamvecid eq zbeam[ind1[ll]]))[0]

     ; interpolate to vertical plane if requested
       trng='Range [km]'
       if altflg and altbeam eq 1 and alt[0] ne -999. then begin 
       ; interpolate image to flight altitude
         mom0=wcrud2altnc(mom0,rrng,alt,altrng=altrng,misval=!values.f_nan,$
              beamdir=strmid(beamstr[ind2],0,4),fltlvl=fltlvl,relief=relief)
         if satflag then begin
           mask=wcrud2altnc(float(mask),rrng,alt,altrng=altrngs,$
                            beamdir=strmid(beamstr[ind2],0,4),fltlvl=fltlvl)                          
           altrngs=altrngs/1000. ; change altitude to km
         endif
         altrng=altrng/1000.   ; change altitude to km
         trng='Altitude [km]'
       endif else if altflg and altbeam eq 2 and alt[0] ne -999. then begin 
       ; interpolate image to flight altitude and strictly vertical plane
         mom0=wcrud2altnc(mom0,rrng,alt,beamvec[*,*,ind4],fltlvl=fltlvl,$
                          altrng=altrng,relief=relief,misval=!values.f_nan)

         if satflag then begin
           mask=wcrud2altnc(float(mask),rrng,alt,altrng=altrngs,$
                            beamvec[*,*,ind4],fltlvl=fltlvl)
           altrngs=altrngs/1000. ; change altitude to km
         endif
         altrng=altrng/1000. ; change altitude to km     
         trng='Altitude [km]'
       endif else if strmid(beamstr[ind2],0,4) eq 'down' then begin
         mom0=rotate(mom0,5)
         if satflag then mask=rotate(mask,5)
         altrng=-reverse(rrng/1000.)
       endif else begin
         altrng=rrng/1000. ; change range to km for plotting
         altrngs=altrng
       endelse

     ; determine reflectivity min/max
       if ldrstr ne '' then begin
         junk2=ldrmin
         junk3=ldrmax
         dBZstr='dB'
       endif else begin
         if zmin eq minrefl then junk2=min(mom0,/nan)>minrefl $
         else junk2=zmin>minrefl
         if zmax eq maxrefl then junk3=max(mom0,/nan)<maxrefl $
         else junk3=zmax<maxrefl
         if junk2 eq junk3 then begin
           junk2=0.9*junk2
           junk3=1.1*junk3
         endif
         dBZstr='dBZ'
       endelse

       if total(finite(mom0)) eq 0 then begin
         message,beamstr[ind2]+' reflectivity image:',/info
         message,'No valid/signal pixels detected. Skip the image.',/info
       endif else begin
       ; display and/or save the reflectivity plot
         if ldrstr eq '' then $
           tit=' reflectivity ('+zpol[ind1[ll]]+', ID: '+zid[ind1[ll]]+')' $
         else tit=ldrstr
         if saveimg[0] ge 0 then begin ; display the reflectivity image plot
           aximage,mom0,time,altrng,/tr,/bar,timeax=timeax,setw=[-1,zcol,0,1], $
             tit=(wcrfnc+', '+beamstr[ind2]+tit)[0],$
             xtit=['UTC',dBZstr,xtt],ytit=trng,asp=asp,topx=xtop,imgz=imgz,$
             widgetid=junk1,/silent,min=junk2,max=junk3+(junk2 eq junk3)
           widgetid[widgetcount+ll]=junk1
           if altbeam gt 0 and alt[0] ne -999. and keyword_set(fltlvl) and $
           altflg then begin
             if imgz ne 1 then begin
               junk4=congrid(time,imgz*n_elements(mom1[0,*,0]),/interp)
               junk5=congrid(alt,imgz*n_elements(mom1[0,*,0]),/interp)/1000.         
               oplot,junk4,junk5,thick=2,lin=1
               if strmid(beamstr[ind2],0,5) ne 'down-' then begin
               ; overplot flight segments with down beam off vertical > 10 deg
                 junk1=reform(90.-abs(atan(beamvec[2,*,ind4]/ $
                       (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2))))*!radeg)
                 junk1=congrid(reform(junk1),imgz*n_elements(mom1[0,*,0]),/interp)
                 ind3=where(junk1 gt 10.)
                 if ind3[0] ne -1 then $
                   oplot,junk4[ind3],junk5[ind3],symsize=0.5,psym=6
               endif
             endif else begin
               oplot,time,alt/1000.,thick=2,lin=1
               if strmid(beamstr[ind2],0,5) ne 'down-' then begin
               ; overplot flight segments with down beam off vertical > 10 deg
                 junk1=reform(90.-abs(atan(beamvec[2,*,ind4]/ $
                       (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2))))*!radeg)
                 ind3=where(junk1 gt 10.)
                 if ind3[0] ne -1 then $
                   oplot,time[ind3],alt[ind3]/1000.,symsize=0.5,psym=6
               endif
             endelse
           endif
           if satflag then begin
             ind3=where(mask lt 0.5*satval)
             if ind3[0] ne -1 then begin 
               ind5=array_indices(mask,ind3)
               oplot,time[ind5[1,*]],altrngs[ind5[0,*]],psym=4,thick=3
             endif
           endif           
           if keyword_set(breaki) then stop          
         endif
         if saveimg[0] ne 0 then begin ; make ps and save in the requested format
           junk4=filen+'.dBZ'+zpol[ind1[ll]]+num2str(ll+1)+'.'+beamstr[ind2]
                 
           setxpswin,99,1,landscape=ori,xysize=xysize,xyoff=xyoff,$
                     psfile=junk4+'.ps',/sil
           message,'Making ps file: '+junk4+'.ps',/info
           loadct1,zcol,rev=0,bf=2,/silent & !p.color=255B

           aximage,mom0,time,altrng,/tr,/bar,timeax=timeax, $
             tit=(strupcase(+beamstr[ind2])+tit)[0],$
             xtit=['UTC',dBZstr,xtt],ytit=trng,topx=xtop,/sil,min=junk2,$
             max=junk3+(junk2 eq junk3)
           if altbeam gt 0 and alt[0] ne -999. and keyword_set(fltlvl) and $
           altflg then begin
             oplot,time,alt/1000.,thick=2,lin=1
             if strmid(beamstr[ind2],0,5) ne 'down-' then begin
             ; overplot flight segments with down beam off vertical > 10 deg
               junk1=reform(90.-abs(atan(beamvec[2,*,ind4]/ $
                     (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2))))*!radeg)
               ind3=where(junk1 gt 10.)
               if ind3[0] ne -1 then $
                 oplot,time[ind3],alt[ind3]/1000.,symsize=0.5,psym=6
             endif
           endif
           if satflag then begin
             ind3=where(mask lt 0.5*satval)
             if ind3[0] ne -1 then begin 
               ind5=array_indices(mask,ind3)
               oplot,time[ind5[1,*]],altrngs[ind5[0,*]],psym=4,thick=1
             endif
           endif           

           xyouts,xypgpos0[0],xypgpos0[1],comm,charsize=1.2,/normal,/noclip
           datestamp,comm=' '+file_basename(junk4),xypgpos=xypgpos,/inches,$
                     chars=0.9,rot=rot
           resetps2x

           if abs(saveimg[0]) eq 2 then begin
             message,'Converting to png: '+junk4+'.png',/info
             spawn,'idlps2png -r'+num2str(res)+' '+junk4+'.ps '+junk4+'.png'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'   
           endif else if pdfflg and abs(saveimg[0]) eq 3 then begin
             message,'Converting to pdf: '+junk4+'.pdf',/info
             spawn,'idlps2pdf '+junk4+'.ps '+junk4+'.pdf'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'     
           endif else if pdfflg eq 0 and abs(saveimg[0]) eq 3 then begin
             message,'Converting to png: '+junk4+'.png',/info
             spawn,'idlps2png -r'+num2str(res)+' '+junk4+'.ps '+junk4+'.png'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'   
           endif
         endif
       endelse
     endfor ; end loop for all other reflectivities
   endif ; end for all other reflectivities

   widgetcount=(where(widgetid eq -1))[0]

   ; Process the rest of the velocities

   ind1=where(dvbeam ne -9)
   if ind1[0] ne -1 then begin
     if keyword_set(acwind) and dvstatusid ne 0 then begin
       blindrg=fix(radrng[0]/rgs)+2 ; # of blind range gates
       junk1=reverse(radrng[0]-(findgen(blindrg)+1)*rgs)
       rrng=[junk1,radrng] 
     endif else rrng=radrng

     for ll=0,n_elements(ind1)-1L do begin
       ind4=(where(beamvecid eq dvbeam[ind1[ll]]))[0]
       ind2=(where(dvbeam[ind1[ll]] eq beamid))[0]
       if strmid(beamstr[ind2],0,4) eq 'down' then begin    
         mom1=dvt[*,*,ind1[ll]]
         altflg=1
         if keyword_set(acwind) and dvstatusid ne 0 then $
           mom1=[replicate(!values.f_nan,1,nprof),$  
                 replicarr(-acwcb[*,ind4],blindrg-2,/dim),$
                 replicate(!values.f_nan,1,nprof),$  
                 mom1]
         if beamstr[ind2] eq 'down' then $
           bartit='m/s, positive is Up' $
         else bartit='m/s, positive is toward WCR'
       endif else if beamstr[ind2] eq 'up' then begin
         altflg=1
         mom1=-dvt[*,*,ind1[ll]]
         bartit='m/s, positive is Up'
         if keyword_set(acwind) and dvstatusid ne 0 then $
           mom1=[replicate(!values.f_nan,1,nprof),$  
                 replicarr(acwcb[*,ind4],blindrg-2,/dim),$
                 replicate(!values.f_nan,1,nprof),$  
                 mom1]
       endif else begin
         altflg=0 ; do not interpolate for altitude if side beams
         mom1=dvt[*,*,ind1[ll]]
         bartit='m/s, positive is toward WCR'
         if keyword_set(acwind) and dvstatusid ne 0 then $
           mom1=[replicate(!values.f_nan,1,nprof),$  
                 replicarr(acwcb[*,ind4],blindrg-1,/dim),$
                 replicate(!values.f_nan,1,nprof),$  
                 mom1]
       endelse

     ; make all NaNs positive (I use -NaNs to identify sub-surface pixels)
       ind3=where(finite(mom1) eq 0 and finite(mom1,/infinity) eq 0)
       mom1[ind3]=!values.f_nan

     ; interpolate to vertical plane if requested
       trng='Range [km]'
       if altflg and altbeam eq 1 and alt[0] ne -999. then begin 
         if nav[0] gt 1 and n_elements(mom1[0,*,0]) ne nprof then junk1=alt0 $
         else junk1=alt 
       ; interpolate image to flight altitude
         mom1=wcrud2altnc(mom1,rrng,junk1,altrng=altrng,misval=!values.f_nan,$
              beamdir=strmid(beamstr[ind2],0,4),fltlvl=fltlvl,relief=relief)
         altrng=altrng/1000. ; change altitude to km
         trng='Altitude [km]'
       endif else if altflg and altbeam eq 2 and alt[0] ne -999. then begin 
         if nav[0] gt 1 and n_elements(mom1[0,*,0]) ne nprof then junk1=alt0 $
         else junk1=alt 
       ; interpolate image to flight altitude and strictly vertical plane
         mom1=wcrud2altnc(mom1,rrng,junk1,beamvec[*,*,ind4],fltlvl=fltlvl,$
                          altrng=altrng,relief=relief,misval=!values.f_nan)

         altrng=altrng/1000. ; change altitude to km     
         trng='Altitude [km]'
       endif else if strmid(beamstr[ind2],0,4) eq 'down' then begin
         mom1=rotate(mom1,5)
         altrng=-reverse(rrng/1000.)
       endif else altrng=rrng/1000. ; change range to km for plotting

     ; determine velocity minmax
       if dvminmax eq maxvel then begin
         ind3=where(mom1 ne 32767)
         junk2=max(abs(mom1[ind3]),/nan)<maxvel
       endif else junk2=dvminmax<maxvel

       if total(finite(mom1)) eq 0 then begin
         message,beamstr[ind2]+' velocity image:',/info
         message,'No valid/signal pixels detected. Skip the image.',/info
       endif else begin
       ; display and/or save the velocity plot
         if saveimg[0] ge 0 then begin ; display the velocity image plot
           aximage,mom1,time,altrng,/tr,/bar,timeax=timeax,setw=setwdv, $
             tit=(strupcase(beamstr[ind2])+' '+dvstatus[dvstatusid])[0]+$
             ' velocity (Nyquist:'+num2str(dvmax[ind1[ll]],fr=1)+$
             ', ID: '+dvid[ind1[ll]]+')',ytit=trng,asp=asp,$
             topx=xtop,/silent,min=-junk2,max=junk2+(junk2 eq 0.),$
             xtit=['UTC',bartit,xtt],imgz=imgz,widgetid=junk1           
           widgetid[widgetcount+ll]=junk1
           if altbeam gt 0 and alt[0] ne -999. and keyword_set(fltlvl) and $
           altflg then begin
             if imgz ne 1 then begin          
               junk4=congrid(time,imgz*n_elements(mom1[0,*,0]),/interp)
               junk5=congrid(alt,imgz*n_elements(mom1[0,*,0]),/interp)/1000.         
               oplot,junk4,junk5,thick=2,lin=1
               if strmid(beamstr[ind2],0,5) ne 'down-' then begin
               ; overplot flight segments with down beam off vertical > 10 deg
                 junk1=reform(90.-abs(atan(beamvec[2,*,ind4]/ $
                       (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2))))*!radeg)
                 junk1=congrid(reform(junk1),imgz*n_elements(mom1[0,*,0]),/interp)
                 ind3=where(junk1 gt 10.)
                 if ind3[0] ne -1 then $
                   oplot,junk4[ind3],junk5[ind3],symsize=0.5,psym=6
               endif
             endif else if nav[0] gt 1 and n_elements(mom1[0,*,0]) ne nprof then $
             begin
               oplot,time,alt0/1000.,thick=2,lin=1
               if strmid(beamstr[ind2],0,5) ne 'down-' then begin
               ; overplot flight segments with down beam off vertical > 10 deg
                 junk1=reform(90.-abs(atan(beamvec[2,*,ind4]/ $
                       (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2))))*!radeg)
                 ind3=where(junk1 gt 10.)
                 if ind3[0] ne -1 then $
                   oplot,time[ind3],alt0[ind3]/1000.,symsize=0.5,psym=6
               endif
             endif else begin
               oplot,time,alt/1000.,thick=2,lin=1
               if strmid(beamstr[ind2],0,5) ne 'down-' then begin
               ; overplot flight segments with down beam off vertical > 10 deg
                 junk1=reform(90.-abs(atan(beamvec[2,*,ind4]/ $
                       (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2))))*!radeg)
                 ind3=where(junk1 gt 10.)
                 if ind3[0] ne -1 then $
                   oplot,time[ind3],alt[ind3]/1000.,symsize=0.5,psym=6
               endif
             endelse
           endif
           if keyword_set(breaki) and ll lt n_elements(ind1)-1 then stop          
         endif
         if saveimg[0] ne 0 then begin ; make ps and save in the requested format
           junk4=filen+'.DV'+num2str(ll+1)+'.'+beamstr[ind2]
           setxpswin,99,1,landscape=ori,xysize=xysize,xyoff=xyoff,$
                     psfile=junk4+'.ps',/sil
           message,'Making ps file: '+junk4+'.ps',/info
           loadct1,setwdv[1],rev=setwdv[2],bf=2,ccol=setwdv[4],/silent
           !p.color=255B

           aximage,mom1,time,altrng,/tr,/bar,timeax=timeax, $
             tit=(strupcase(beamstr[ind2])+' '+dvstatus[dvstatusid])[0]+$
             ' velocity (Nyquist:'+num2str(dvmax[ind1[ll]],fr=1)+$
             ', ID: '+dvid[ind1[ll]]+')',ytit=trng,$
             topx=xtop,/silent,xtit=['UTC',bartit,xtt],$
             min=-junk2,max=junk2+(junk2 eq 0.)
           if altbeam gt 0 and alt[0] ne -999. and keyword_set(fltlvl) and $
           altflg then $
             if nav[0] gt 1 and n_elements(mom1[0,*,0]) ne nprof then begin
               oplot,time,alt0/1000.,thick=2,lin=1 
               if strmid(beamstr[ind2],0,5) ne 'down-' then begin
               ; overplot flight segments with down beam off vertical > 10 deg
                 junk1=reform(90.-abs(atan(beamvec[2,*,ind4]/ $
                       (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2))))*!radeg)
                 ind3=where(junk1 gt 10.)
                 if ind3[0] ne -1 then $
                   oplot,time[ind3],alt0[ind3]/1000.,symsize=0.5,psym=6
               endif
             endif else begin
               oplot,time,alt/1000.,thick=2,lin=1
               if strmid(beamstr[ind2],0,5) ne 'down-' then begin
               ; overplot flight segments with down beam off vertical > 10 deg
                 junk1=reform(90.-abs(atan(beamvec[2,*,ind4]/ $
                       (sqrt(beamvec[0,*,ind4]^2+beamvec[1,*,ind4]^2))))*!radeg)
                 ind3=where(junk1 gt 10.)
                 if ind3[0] ne -1 then $
                   oplot,time[ind3],alt[ind3]/1000.,symsize=0.5,psym=6
               endif
             endelse

           xyouts,xypgpos0[0],xypgpos0[1],comm,charsize=1.2,/normal,/noclip
           datestamp,comm=' '+file_basename(junk4),xypgpos=xypgpos,/inches,$
                     chars=0.9,rot=rot
           resetps2x

           if abs(saveimg[0]) eq 2 then begin
             message,'Converting to png: '+junk4+'.png',/info
             spawn,'idlps2png -r'+num2str(res)+' '+junk4+'.ps '+junk4+'.png'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'   
           endif else if pdfflg and abs(saveimg[0]) eq 3 then begin
             message,'Converting to pdf: '+junk4+'.pdf',/info
             spawn,'idlps2pdf '+junk4+'.ps '+junk4+'.pdf'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'     
           endif else if pdfflg eq 0 and abs(saveimg[0]) eq 3 then begin
             message,'Converting to png: '+junk4+'.png',/info
             spawn,'idlps2png -r'+num2str(res)+' '+junk4+'.ps '+junk4+'.png'
             message,'Removing '+junk4+'.ps',/info
             file_delete,junk4+'.ps'   
           endif
         endif
       endelse
     endfor ; end loop for all other velocities
   endif ; end for all other velocities

  ; ---------------------------------------------------------------------------- 
  ; Select how to continue

    if saveimg[0] ge 0 then begin
      junk1=''
      print,format='($,"Continue/Stop/Quit (c,s,q)")'
      read,junk1
 ;    read,prompt='Continue/Stop/Quit (c,s,q): ',junk1
      junk1=strmid(strlowcase(strtrim(junk1,2)),0,1)

      if junk1 eq 'q' then begin
        ncdf_close,ncid
        return
      endif else if junk1 eq 's' then stop
      ncdf_close,ncid

    ; Erase the plots before continuing if not requested otherwise

      if keyword_set(noerase) then noerase=1 else noerase=0
      ind1=where(widgetid ne -1)

      if noerase eq 0 and ind1[0] ne -1 then begin
        for ll=0,n_elements(ind1)-1 do $
          if widgetid[ind1[ll]] ge 10000 then wdelete,widgetid[ind1[ll]]-10000 $
          else widget_control,widgetid[ind1[ll]],/destroy
        widgetid[*]=-1 
      endif
    endif
  endfor ; end file loop; kk

end 
