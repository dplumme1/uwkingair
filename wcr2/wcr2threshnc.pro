;$Id: wcr2threshnc.pro,v 1.2 2015/03/14 00:08:48 haimov Exp haimov $
;+
; NAME:
;       WCR2THRESHNC
;
; PURPOSE:
;       Return thresholded reflectivity in dBZ and optionally the matching
;       Doppler velocity.
;       If reflectivity mask is given range gates with sub-surface data
;       points only are removed and the rest that are left are assigned no
;       target value
;
; CATEGORY:
;       wcr2tools
;
; CALLING SEQUENCE:
;       result=wcr2threshnc(rf,noise,radrng,sfactor)
;
; INPUTS:
;       Rf:      fltarr(ngate,nprof); reflectivity factor in mm^6/mm^3
;                  ngate-number of range gates, nprof-number of profiles
;                  ATTN: Rf should not have missing values
;       Noise:   fltarr(nprof), Standard deviation of the reflectivity
;                  noise in equivalent mm^6/m^3 at 1 km or
;                fltarr(nprof), the full noise in mm^6/m^3 at 1 km
;                 (requires MNOISE to be defined and NAV[0]>1)
;                  ATTN: Noise should not have missing values
;       Radrng:  fltarr(rangegates) radar range gates in meters
;       Sfactor: float, factor applied to StDev of the noise to form 
;                 the threshold (must be greater than 0., default is 3.)
;                
; KEYWORD PARAMETERS:
;       DVC:     input, fltarr(rg,np); Doppler velocity field corrected for 
;                  aircraft motion and matching Rf 
;                  (should be used in conjuction with keyword DVT)  
;                  ATTN: DVC should not have missing values
;       DVT:     output, fltarr(rg,np); tresholded Dv in m/s
;       DVMax:   input,float; Nyquist velocity in m/s for DVC
;                  Note: this is optional parameter when NAV>1 and DVC
;                        input is given; in this case the DVC averaging
;                        is done in I and Q phase rather than velocity
;                        to avoid averaging errors due to folding; 
;                        If DVMax is not given velocity values are averaged.
;       MASKR:   intarr(ngates,nprof), reflectivity mask; 
;                  (should be used in conjuction with keyword MASKT)
;                  MASK description is given in wcr2mask2.pro  
;       MASKT:   output, intarr((ngates/nrga)*nrga,(nprof/npa)*npa)
;                  returns MASK matching the thresholded reflectivity
;                  ATTN: bit 0 is set according to sfactor 
;                        bit 1 and 2 are set to 0
;       MISVAL:  input,float, missing value for the output (reflectivity)
;                  and velocity DVT (default -32767.)
;       MNOISE:  input,fltarr(nprof), mean noise in mm^6/m^3 at 1km, which is
;                  subtracted from reflectivity.  If this is given it is
;                  assumed that Noise is the full noise; ignored if NAV[0]=1
;                  ATTN: Do not input MNOISE if NAV[0]>1 and the full noise
;                        is not available. Leaving MNOISE while using Noise 
;                        St.Dev as input for Noise will assume it is the 
;                        full noise !
;       NAN:     input, when set replaces all missing values with NaNs if 
;                missing values are not NaNs
;       NAV:     input, [npa,nrga], number of profiles(npa, default 1) and 
;                  range gates(nrga, default is 1) to be averaged; or
;                  if MNOISE is given:
;                  [npa,nrga,nsw], where nsw is smoothing window size for
;                  calculating the moving-averaged mean noise and StDev
;                  (default for nsw is 101; it expected that npa<30 and nrga<10)
;       NFILTER: input, when set implements further noise reduction by removing
;                       isolated pixels; noise detection uses smooth(dBZ,[nrg,npr]),
;                       where nrg is # of gates and npr is # of prof to smooth 
;                       (default is no filtering):
;                  integer (odd number) to use for smooth(dBZ,[1,nfilter])
;                  intarr(2) - use smooth(dBZ,[nfilter[0],nfilter[1]])
;                  Note: this filtering causes some minor changes in the targets
;                        boundaries - use with caution; nfilter is not applied
;                        to MASKT. Typical values are from [1,3] to [5,5]
;       QUIET:   input, when set supress informational messages
;       RRADRNG: output,fltarr((ngates/nrga)*nrga), resampled range axis of 
;                  Rf when nrga > 1, otherwise returns Radrng; if reflectivity
;                  mask (MASKR) is given range gates with sub-surface pixels
;                  only are removed from RRADRNG.
;       RTIME:   output, dblarr((nprof/npa)*npa), resampled Time when NAV[0]>1 
;       TIME:    input,dblarr(nprof); time stamps for every profile in 
;                  Unix seconds(used in conjuction with RTIME)
;       ZMIN:    output, fltarr(nprof), minimum detectable signal at 1 km based
;                  on the requested thresholding           
;       HELP:    show this text (use print,wcr2threshnc(/help) )
;
; OUTPUTS:
;
; SIDE EFFECTS:
;       If Rf size is not exactly divisible by NAV Reflectivity size is 
;       truncated to the maximum divisible integer lower than the input size. 
;       Thus it is possible that not all input profiles are used.
;       
; NOTES:
;       If Noise is the full noise (and therefore NAV[0]>1 and MNOISE is given)
;       the thresholding is done by first restoring the removed mean noise in
;       the reflectrivity.  Then a new moving averaged noise as well as the
;       noise St.Dev. are calculated using "boxcar" window NAV[0]*nsw. And
;       finally, the detection threshold of StDev*Sfactor is applied.
;
;       ATTENTION: If nav[1]>1, average several range gates this routine
;                  deviates from what a radar with equivalent resolution in
;                  range will measure, especially when the range_resolution
;                  to range gate distance ratio is large. The radar averages
;                  the scatterer returns in the larger resolution volume 
;                  without power range correction. 
;                  In this routine the range gate averaging is done after the
;                  range correction is applied and therefore the new averaged
;                  reflectivity in the new larger range gate will more closely
;                  be represented with the true reflectivity associated with
;                  the geometric center of the gate. This metters for WCR since
;                  range gates of more than 200 m are significantly large when
;                  compared to the typical maximum range of 4-5 km.
;
; EXAMPLES:
;   ; Load the down beam from ASCII RF18 data file using Applanix and
;   ; processed after 28 March 2012 (full noise var was added) 
;     dnb=readwcr2beamnc('WCR.ASCII12av.20120303.175529_181439.CPP-H1H2V2.nc',2)
;
;   ; Threshold 3 St.Dev, no averaging, no nfilter
;     z2=wcr2threshnc(dnb.z,dnb.noisesd,dnb.radrng,3,dvc=dnb.dv,$
;                     maskr=dnb.zmask,dvt=dvt,maskt=maskt,/nan)
;
;   ; Average 10 profiles, threshold above 3 StDev and apply nfilter
;     z2a=wcr2threshnc(dnb.z,dnb.noise,dnb.radrng,3,nav=[10,1,101],nfilt=[5,5],$
;             dvc=dnb.dv,maskr=dnb.zmask,mnoise=dnb.noisemn,rradrng=rnga,$
;             time=dnb.time,rtime=rtimea,dvt=dvta,maskt=maskta,/nan)
;   ; Plot z2 and z2a
;     aximage,z2,/tr,dnb.time mod 86400,dnb.radrng,/timeax,/order,$
;             tit=dnb.fn,xtit='UTC',ytit='Range [m]',setw=[-1,41]
;     aximage,z2a,/tr,rtimea mod 86400,rnga,/timeax,/order,$
;             tit=dnb.fn,xtit='UTC',ytit='Range [m]',setw=[-1,41]
;-
; MODIFICATION HISTORY:
;       Written by:  Samuel Haimov, September 2007 
;                    (modified version of wcrcalthreshnc.pro)
;       S.H., Apr 2012 - major changes
;       S.H., Jan 2013 - changed nsw*npa to nsw when nav[0]>1 and mnoise is given             
;       S.H., Oct 2013 - added DVMax keyword and modified the DVC averaging
;                        when NAV>1
;       S.H., Apr 2014 - saved rngcor for use everywhere it is needed; added
;                        a note for nrga>1 ; May 14 found that while making
;                        the changes for rngcor line 320 (calv=rebin(calv,nr,np)
;                        was mistakenly moved after line 321 (if nrga gt 1 .. )-
;                        this is fixed now
;
; Copyright (C) 2007-present, Samuel Haimov, Dept. of Atmos. Sci., University 
; of Wyoming.  This software may be used, copied, or redistributed as long as 
; it is not sold and this copyright notice is reproduced on each copy made.
; This routine is provided as is without any expressed or implied warranties
; whatsoever. Other limitations apply as described in the WCR2TOOLS Readme file.

function wcr2threshnc,rf,noise,radrng,sfactor,dvc=dvc,maskr=maskr,dvmax=dvmax, $
                      dvt=dvt,maskt=maskt,misval=misval,mnoise=mnoise,$
                      nan=nan,nav=nav,nfilter=nfilter,rradrng=rradrng,$
                      rtime=rtime,time=time,zmin=zmin,quiet=quiet,help=help
;  on_error,2

  if keyword_set(help) then begin
    doc_library,'wcr2threshnc'
    return,''
  endif

  if n_params() lt 3 then begin
    message,"Usage: res=wcr2threshnc(rf,noise,radrng,sfactor)",/info
    return,''
  endif else if n_params() eq 3 then sf=3. else sf=sfactor

; Check the input and re-size the data if needed

  if sf le 0. then message,'Sfactor must be gerater than 0.'

  if n_elements(dvc) gt 0 then begin
    if total(size(dvc,/dim)) ne total(size(rf,/dim)) then $
      message,'Invalid input for DV; must be the same size as Rf.'
    dvflg=1
    if n_elements(dvmax) eq 0 then dvmax=0 else begin
      dvmax=abs(dvmax[0])
      if dvmax lt max(abs(dvc)) then $
        message,'DVMax is less than max(abs(DVC))'
    endelse
  endif else dvflg=0

  if n_elements(maskr) gt 0 then begin
    if total(size(maskr,/dim)) ne total(size(rf,/dim)) then $
      message,'Invalid input for MASKR; must be the same size as Rf.'
    maskflg=1
  endif else maskflg=0

  if n_elements(misval) eq 0 then misval=-32767.
  if n_elements(time) eq 0 then rtime=-32767. else rtime=time[0]

  if n_elements(nav) eq 0 then begin 
    nrga=1 & npa=1 & nsw=101
  endif else if n_elements(nav) eq 1 then begin
    nsw=101 & nrga=1 & npa=nav
  endif else if n_elements(nav) eq 2 then begin
    nsw=101 & nrga=nav[1] & npa=nav[0]
  endif else begin
    nsw=nav[2] & nrga=nav[1] & npa=nav[0]    
  endelse
  
  if n_elements(nfilter) eq 0 then nf=0 else nf=nfilter

; If reflectivity mask is given check for gates with sub-surface pixels 
; for all profiles and remove them

  if maskflg then begin
    ind1=where(maskr eq 1024) ; find sub-surface pixels
    if ind1[0] ne -1 then begin
      junk1=rf
      junk1[ind1]=0.
      junk1=total(junk1,2)
      ind2=where(junk1 ne 0.)
      if ind2[0] ne -1 then begin ; should never be -1
        calv=rf[ind2,*]
        if dvflg then dvt=dvc[ind2,*]
        maskt=maskr[ind2,*]
        rng=radrng[ind2]
      endif else begin
        calv=rf 
        if dvflg then dvt=dvc
        maskt=maskr
        rng=radrng
      endelse
    endif else begin
      calv=rf 
      if dvflg then dvt=dvc
      maskt=maskr
      rng=radrng
    endelse
  endif else begin
    calv=rf
    if dvflg then dvt=dvc
    rng=radrng
  endelse

  nprof =n_elements(calv[0,*])
  ngate =n_elements(calv[*,0])
  np    =nprof/long(npa) ; output profiles
  nr    =n_elements(rng)/long(nrga)           ; ouput range gates
  indp  =np*npa-1L
  indr  =nr*nrga-1L

; Resize

  anoise=noise[0:indp]                  ; noise cut to multiple of npa
  arng  =rng[0:indr]                    ; radrng cut to multiple of nrga
  rngcor=arng^2*1.e-6 
  if rtime[0] ne -32767. then rtime =time[0:indp]
  calv=calv[0:indr,0:indp]              ; rf cut to (np*npa)*(nr*nrga)
  if dvflg then dvt=dvt[0:indr,0:indp]

  if n_elements(mnoise) eq 0 and npa gt 1 then begin
    if keyword_set(quiet) eq 0 then begin
      message,'NAV[0]>1 and MNOISE is not given !',/info
      message,'Will assume NOISE is the StDev of the noise.',/info
    endif
  endif else if n_elements(mnoise) ne nprof and npa gt 1 then begin
    message,'MNOISE must have the same number of profiles as Rf'
  endif else if n_elements(mnoise) eq nprof and npa gt 1 then begin
  ; restore the removed mean noise in the reflectivity
  ; rngcor=replicarr(arng^2*1.e-6,indp+1L); range correction
  ; calv=calv/rngcor+replicarr(mnoise[0:indp],nindr+1L,dim=1)
  ;  rngcor=arng^2*1.e-6 ; move this line earlier, may use it several times
    calv=calv+rngcor#mnoise[0:indp]
  ; calculate the new mean noise and subtract it from reflectivity
  ;  nw=npa*nsw/2  ; changed it Jan 2013; too long smoothing, which may cause
  ;  thresh=smooth(anoise,nsw*npa) ; problem when there are not enough profs.
    nw=nsw/2                       ; I'll leave it to nsw, use larger value
    thresh=smooth(anoise,nsw)      ; given by the input if needed

    thresh=[reverse(thresh(nw:nw*2-1)),thresh(nw:nprof-nw-1),$
            reverse(thresh(nprof-nw*2:nprof-nw-1))]
  ; calv=(calv-replicarr(thresh,ngate,/dim))*rngcor
    calv=calv-rngcor#thresh
  endif

; Process ...

  if npa eq 1 and nrga eq 1 then begin
  ; -----------------------------------------------------------------------
  ; Process for non-averaging 

    ; Prepare the threshold
 
      if keyword_set(quiet) eq 0 then $
        message,'Assume NOISE input is the StDev of the noise.',/info

    ; thresh=replicarr(anoise*sf,nr,dim=1)*replicarr(rngcor,np)
      zmin=anoise*sf
      thresh=rngcor#zmin  
      zmin=10.*alog10(zmin)

    ; If mask is given remove the preset target bits (1,2,3 StDev)

    if maskflg then begin
      ind4= where((maskt AND 1) eq 1)
      if ind4[0] ne -1 then maskt[ind4]=maskt[ind4]-1
      ind4= where((maskt AND 2) eq 2)
      if ind4[0] ne -1 then maskt[ind4]=maskt[ind4]-2
      ind4= where((maskt AND 4) eq 4)
      if ind4[0] ne -1 then maskt[ind4]=maskt[ind4]-4
    endif

  endif else begin
  ; -----------------------------------------------------------------------
  ; Process for averaging 

    if npa gt 1 then begin
      anoise=rebin(anoise,np)                           ; adjust noise
      if rtime[0] ne -32767. then rtime=rebin(rtime,np) ; resample time
    endif
     
    ; Rebin cannot handle NaNs properly and it should not be used if
    ; there are missing values in Rf or Dv (see wcr2calthreshnc.pro)
  
    calv=rebin(calv,nr,np)
    if nrga gt 1 then begin
      arng=rebin(arng,nr)               ; adjust range
      rngcor=arng^2*1.e-6
    endif
    
    if dvflg then begin
      if dvmax eq 0 then begin
        message,'DVmax not given; will average velocities(bias due to folding can occur)',/info
        dvt=rebin(dvt,nr,np)
      endif else begin ; averaging in phase will avoid possible errors 
                 ; if DVC contains folded values
        junk3=!pi/dvmax
        junk1=rebin(cos(junk3*dvt),nr,np)
        junk2=rebin(sin(junk3*dvt),nr,np)
        dvt=atan(junk2,junk1)/junk3
      endelse
    endif
 
    ; Prepare the threshold from Noise

    if npa gt 1 then begin
      if n_elements(mnoise) eq nprof then begin
      ; anoise must be the non-smoothed noise
        if keyword_set(quiet) eq 0 then $
          message,'Assume NOISE input is the non-averaged full noise.',/info
        ind1   =np-nsw-1
        thresh=fltarr(np)
        for k =0L,ind1 do thresh(k)=std(anoise(k:k+nsw))
        thresh=thresh*sf
        thresh=[replicate(thresh[0],(nsw-1)/2),thresh[0:ind1],$
                replicate(thresh[ind1],(nsw-1)/2+1)]
      endif else begin ; noise must be the StDev of the noise
        thresh=anoise/sqrt(npa)*sf
      endelse
    endif else begin        ; Noise must be the StDev of the noise
    ; anoise must be the StDev of the noise
      if keyword_set(quiet) eq 0 then $
        message,'Assume NOISE input is the StDev of the noise.',/info
      thresh=anoise*sf
    endelse

    ; adjust the thresh for the additional averaging if nrga > 1
    ; (e.g., if nrga=2 reduce the StDev by 1.5dB (sqrt(nrga)))

    if nrga gt 1 then thresh=thresh/sqrt(nrga)

    ; thresh=replicarr(thresh,nr,dim=1)*replicarr(arng^2*1.e-6,np)
    zmin=thresh
    thresh=rngcor#zmin
    zmin=10.*alog10(zmin)

    ; Resample the mask

    if maskflg then begin
      maskt=maskt[0:indr,0:indp]         ; cut to (np*npa)*(nr*nrga)

    ; find surface gates resample them and restore in maskt after averaging

      ind1=where(maskt eq 512)
      if ind1[0] ne -1 then begin
        ind2=ind1/(indr+1L)     ; map ind1 to profiles
        ind3=ind1 mod (indr+1L) ; map ind2 to range gates

        junk1=max(ind3,junk2)  ; find how many gates are marked as surface clutter
        nsc=where(maskt(0:junk1,junk2) AND 256)
        if nsc[0] ne -1 then nsc=n_elements(nsc) else nsc=-1

      ; generate surface index array for every profile; if not present leave NaN

;       commented next line _ cannot use NaNs since there could be missing
;       surface points (in case of side beam NaNs would be the majority)
;       junk1=make_array(indp+1L,value=!values.f_nan)
        junk1=make_array(indp+1L,value=-32767.)
        junk1[ind2]=ind3  ; save surface gate for every profile it is present
        ind1=(round(rebin(junk1,np))) ; calc the averaged surface gate index
        ind4=where(ind1 le indr and ind1 gt 0)
        if ind4[0] ne -1 then ind1=ind1[ind4]
        ind2=lindgen(np) ; generate the prof. index array after averaging
;        ind3=where(finite(ind1))    ; find which are finite 
;        if ind3[0] ne -1 then begin ; and get the prof and rg indices for those
;          ind1=ind1[ind3]
;          ind2=ind2[ind3]
;        endif 
        ind2=ind2[ind4]     
      endif else begin ind3=-1 & nsc=-1 & endelse

    ; remove the surface clutter, surface, and sub-surface mask bits
  
      ind4=where((maskt AND 256) eq 256) 
      if ind4[0] ne -1 then maskt[ind4]=maskt[ind4]-256
      ind4= where(maskt eq 512)
      if ind4[0] ne -1 then maskt[ind4]=0
      ind4= where(maskt eq 1024)
      if ind4[0] ne -1 then maskt[ind4]=0

    ; resample maskt using nearest sample
  
      maskt=rebin(maskt,nr,np,/sample)           

      if ind3[0] ne -1 then begin

      ; determine range gate indices (ind1) for the surface if nrga > 1

        if nrga gt 1 then ind1=value_locate(arng,rng[ind1])

      ; re-mask surface clutter

        if nsc ne -1 then begin
          nsc=ceil(float(nsc)/float(nrga))
          for k=0,nsc-1L do $        
            maskt[(ind1-(nsc-k))>0,ind2]=maskt[(ind1-(nsc-k))>0,ind2]+256
        endif

      ; re-mask sub-surface pixels
        ind4=nr-1
        for k=0,n_elements(ind2)-1L do maskt[ind1[k]:ind4,ind2[k]]=1024

      ; re-mask surface pixels

        maskt[ind1,ind2]=512   ; save ind1 and ind2 for possible later use

      endif

      ; remove bits 0, 1, 2 that are originaly set in maskr

      ind4= where((maskt AND 1) eq 1)
      if ind4[0] ne -1 then maskt[ind4]=maskt[ind4]-1
      ind4= where((maskt AND 2) eq 2)
      if ind4[0] ne -1 then maskt[ind4]=maskt[ind4]-2
      ind4= where((maskt AND 4) eq 4)
      if ind4[0] ne -1 then maskt[ind4]=maskt[ind4]-4
    endif

  endelse ; end processing for averaging
  ; -----------------------------------------------------------------------

; Threshold reflectivity

  ind4=where(calv lt thresh,compl=ind3) 

  if ind4[0] ne -1 then begin
    calv[ind3]=10.*alog10(calv[ind3])
    calv[ind4]=-32767.
    if dvflg then dvt[ind4]=-32767.

  ; if mask is given set 0 bit to mark target pixels

    if maskflg then begin
      maskt[ind3]=maskt[ind3]+1
    ; repair surface and sub-surface mask pixels if any
      ind4=where(maskt eq 513 or maskt eq 1025)
      if ind4[0] ne -1 then maskt[ind4]=maskt[ind4]-1
    endif
  endif else calv=10.*alog10(calv)

; Assign misval to all subsurface pixels

  if maskflg then begin
    ind4=where(maskt eq 1024)
    if ind4[0] ne -1 then begin
      calv[ind4]=-32767.
      if dvflg then dvt[ind4]=-32767.
    endif
  endif

; Remove pixels identified as noise when nfilter is applied

  if nf[0] ne 0 then begin ; use -32767. assigned above
  ;  nf=3
  ;  junk1=calv
  ;  for k=1,np-2 do junk1[*,k]=total(junk1[*,k-1:k+1],2)/3.
  ;  ind3=where(junk1 lt -8000.)
  ;  if ind3[0] ne -1 then begin & $
  ;    junk1=calv & $
  ;    junk1[ind3]=!values.f_nan & $
  ;  endif

  ; The following should be the same as above but is not due to a bug in smooth;
  ; smooth has big rounding errors (4th significant digit) for large data range. 
  ; For this particular aplication as a detection filter smooth detects more
  ; noise pixels than the above moving-averaged filter. 

     
    if n_elements(nf) eq 1 then junk1=smooth(calv,[1,nf[0]]) $
    else begin
    ; junk1=smooth(calv,[1,nf[1]])
    ; junk1=smooth(junk1,[nf[0],1])
      junk1=smooth(calv,[nf[0],nf[1]]) ;does the same as above
    endelse    
    ind4=where(junk1 lt -8000.)  ; threshold at ~25% of -32767. 
    if ind4[0] ne -1 then begin
      junk1=calv
      junk1[ind4]=-32767.
    endif

  ; restore surface pixels if any (some may be lost due to noise filtering)

    if maskflg then begin
      ind3=where(maskt eq 512)
      if ind3[0] ne -1 then begin
        junk1[ind3]=calv[ind3]
        calv=junk1
      endif else calv=junk1
    endif else calv=junk1
  endif 

  ; Assign the proper misval

   ind3=where(calv eq -32767.)
   if ind3[0] ne -1 then begin
     if keyword_set(nan) and finite(misval) then begin
       calv[ind3]=!values.f_nan
       if dvflg then dvt[ind3]=!values.f_nan
     endif else begin
       calv[ind3]=misval
       if dvflg then dvt[ind3]=misval
     endelse
  endif

  rradrng=arng
  return,calv
end
