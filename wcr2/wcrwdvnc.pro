;$Id: wcrwdvnc.pro,v 1.1 2015/03/14 00:08:48 haimov Exp haimov $
;+
; NAME:
;	WCRWDVNC
;
; PURPOSE:
;
;    Adjust the airborne radar radial Doppler velocity from the Up/Down-pointing
;    antenna for horizontal wind contribution (if any). The routine uses Doppler 
;    velocity already corrected for the aircraft motion contribution. Velocity
;    sign is left with radar convention (positive velocity is toward the radar).
;
; CATEGORY:
;
;	wcrtools
;
; CALLING SEQUENCE:
;
;	wdv=wcrwdvnc(dvc,beamvector,hwind,hwindp=hwindp,help=help)
;
;
; INPUT:
;
;    Dvc:        fltarr(nrg,nprof), nrg range gates and nprof profiles of
;                corrected for aircraft motion WCR Doppler velocity from
;                Up or Down antennas
;    Beamvector: fltarr(3,nprof), East,North,Up WCR beam unit vector
;    Hwind:      Horizontal winds applicable to the Dvc field:
;                - fltarr(2), East(u) and North(v) components of the mean wind
;                  field applicable to all radar profiles;
;                - fltarr(nprof,2), East,North wind field along flight track
;
; KEYWORD PARAMETERS:
;
;    HWINDP:    fltarr(nrg,2), horizontal winds vertical profile normalized 
;               to the winds at 1st rangegate (i.e., hwindp[0,*]=1.; nrg range 
;               gates match the range gates of Dvc
;    HELP:      Show this text
;
; PROCEDURE:
;    DVC is corrected by the given HWIND assuming constant horizontal winds
;    in the vertical.  If this assumption does not hold a vertical profile
;    HWINDP can be used.  The estimated vertical Doppler velocity is calculated
;    from the following equation:
;
;    dvw = (dvc+[hwind_u,hwind_v].beamvector_u_v)/abs(beamvector_w) ;
;    where beamvector_w is the vertical component of beamvector,
;          hwind_u, hwind_v are the East and North wind components,
;          and . designates a dot product
;    
;    Note that if there are missing values in dvc and they are not NaNs wdv 
;    won't preserve them. Since the modifications in the missing values are
;    expected to be small (mssing values that are not NaNs are usually large
;    negative numbers) it should still be possible to identify them. One can
;    always replace the missing values with NaNs.
;
; EXAMPLES:
;    ; Adjust a WCR non truly vertical down beam Doppler velocity (corrected 
;    ; for AC motion) using AC measured horizontal winds; HiCu03 case:
;
;    acfile='/net/bat/R1/data/kingair_data/hicu03/processed/20030826.c25.nc'
;    wcrf='/net/fox/fox1/WCR/HICU03/aug26/W2003-08-26-17-16-52' ;v-beam is down
;
;    openwcrf,wcrf
;    readwcrdata,172600,172715,/timein,rawt=wcrt
;    dvv=wcrdopcorm(wcrt,acncf=acfile,beaminit=3,acdata=acdata,radbeam=wcrb)
;    beamvector=ac2wcr(wcrb,acdata.time,wcrt) ; resample radbeam onto wcrtime
;
;    dvvc=wcrwdvnc(dvv,beamvector,[mean(acdata.uw),mean(acdata.vw)])
;
;-
; MODIFICATION HISTORY:
;     Written by:  Samuel Haimov and Rick Damiani, May 2004
;
; Copyright (C) 2004-present, Samuel Haimov, Dept. of Atmos. Sci., University of 
; Wyoming.  This software may be used, copied, or redistributed as long as it 
; is not sold and this copyright notice is reproduced on each copy made. This
; routine is provided as is without any expressed or implied warranties
; whatsoever.  Other limitations apply as described in the Readme file.
;

function wcrwdvnc, dvc,beamvector,hwind,hwindp=hwindp,help=help

  on_error,2

  if keyword_set(help) then begin
    doc_library,'wcrwdvnc'
    return,''
  endif

  if (n_params() lt 3) then begin
     message,'Usage: wdv=wcrwdvnc(dvc,beamvector,hwind,hwindp=hwindp)',/info
     return,-1
  endif

; Find sizes and signs

  nrg  =n_elements(dvc[*,0])
  nprof=n_elements(dvc[0,*])
  bw   =replicarr(beamvector[2,*],nrg,/dim)

; Match the given winds to dvc

  if n_elements(hwind) eq 2 then begin
    hw=replicarr(replicarr(hwind,nprof,/dim),nrg,/dim)
  endif else if n_elements(hwind) eq 2L*nprof then begin
    hw=replicarr(hwind,nrg,/dim)
  endif else message,'Invalid Hwind input'

; Adjust for horizontal wind vertical profile

  if n_elements(hwindp) eq nrg then begin
    hw=hw*replicarr(replicarr(hwindp,nprof),2)
  endif else if n_elements(hwindp) eq 2L*nrg then begin 
    hw[*,*,0]=hw[*,*,0]*replicarr(hwindp[*,0],nprof)
    hw[*,*,1]=hw[*,*,1]*replicarr(hwindp[*,1],nprof)
  endif else if n_elements(hwindp) ne 0 then message,'Invalid Hwindp input' 

; Calculate horizontal wind contribution into the beam

  if n_elements(hwindp) eq 0 then begin
    wdv=fltarr(nprof)
    for k=0L,nprof-1 do $
      wdv[k]=beamvector[0:1,k]##reform(hw[0,k,*],1,2)
    wdv=replicarr(wdv,nrg,/dim)
  endif else begin
    wdv=fltarr(nrg,nprof)
    for l=0L,nrg-1 do $
      for k=0L,nprof-1 do $
        wdv[l,k]=beamvector[0:1,k]##reform(hw[l,k,*],1,2)
  endelse

; Calculate radar vertical velocity and quit

; wdv=(sign(bw)*dvc-wdv)/bw
  wdv=(dvc+wdv)/abs(bw)

  return,reform(wdv)

end
