;$Id: wcr2wunfoldnc.pro,v 1.1 2015/03/14 00:08:48 haimov Exp haimov $
;+
; NAME:
;      WCR2WUNFOLDNC
;
; PURPOSE:
;      Returns WCR Doppler velocity centered (and possibly unfolded) at 
;      given shift_velocity into the WCR beam. In some cases (e.g., 
;      Doppler measurements in horizontal planes) the aircraft measured
;      wind can be used successfully.
;
; CATEGORY:
;      wcr2tools
;
; CALLING SEQUENCE:
;      res=wcr2wunfoldnc(wcrdvf,vmax,wcrwind)
;    or
;      res=wcr2wunfoldnc(wcrdvf,vmax,wind,wcrbeam)
;
; INPUTS:
;      Wcrdvf:  fltarr(rg,prof),WCR measured radial Doppler velocity field
;               corrected for AC motion in m/s; rg: range gates, pro: radar 
;               profiles; positive value means toward the radar
;      Vmax:    float, Doppler Nyquist velocity in m/s
;      -------
;      Wcrwind: fltarr(prof),shift_velocity into the WCR beam
;               Positive values means away from the radar
;
;      -------  or
;
;      Wind:    fltarr(3),      Mean shift_velocity vector [EW,NS,Up] 
;                               to be applied for all vel. points or
;               fltarr(prof,3), Shift_velocity vector,[EW,NS,Up], for each 
;                               radar profile, in m/s 
;      Wcrbeam: fltarr(3,prof), WCR beam vectors in Earth corrdinate system
;      -------
;
; KEYWORD PARAMETERS:
;      WINDSCALE: fltarr(rg), if defined scales the wind given in wcrwind/wind 
;                 along the range gate axis
;      HELP:      show this text
;
; SIDE EFFECTS:
;      The function returns an empty string if the usage is incorrect
;
; NOTES:
;      All input variables are available from the WCR processed netDF files
;
; PROCEDURE:
;      This routine shifts the Doppler velocity data according to the shift
;      velocity(wcrwind or wind) given as an input.  If wind estimate is used 
;      it is assumed that its values are not changing along the range gate 
;      dimension of the measured Doppler field.  It is also assumed that any 
;      folding that takes place in the input Doppler field is caused by the
;      input shift velocity.  In this case this procedure will unfold any 
;      folded velocities even in cases of multiple folding as long as the
;      number of foldings is constant for the corrected region. The
;      unfolding may not be successful if the folded velocities are affected
;      by other phenomena, e.g., perturbations to the velocity field exceeding
;      the unambiguous range, strong shear, etc.
;
;      Windscale keyword is an attempt to deal with folding in the presence 
;      of a strong shear.  The combination of the wind defined at 0th range
;      gate by wcrwind or wind and windscale will generate a wind profile,
;      which then is used to offset/shift the Doppler data accordingly.
; 
;
; EXAMPLES:
;
;-
; MODIFICATION HISTORY:
;      Written by: Samuel Haimov and Rick Damiani, April 2011
;
; Copyright (C) 2011-  Samuel Haimov, Dept. of Atmos. Sci., University of 
; Wyoming.  This software may be used, copied, or redistributed as long as it 
; is not sold and this copyright notice is reproduced on each copy made. This
; routine is provided as is without any expressed or implied warranties
; whatsoever.  Other limitations apply as described in the Readme file.
;

function wcr2wunfoldnc,wcrdvf,vmax,wind,wcrbeam,windscale=windscale,help=help 

  on_error,2

  if keyword_set(help) then begin
    doc_library,'wcr2wunfoldnc'
    return,''
  endif

; Check input and calculate wind/shift_velocity component into wcr beam

  nrg=n_elements(wcrdvf[*,0])
  np=n_elements(wcrdvf[0,*])

  if n_params() lt 3 then begin
    message,"Usage: res=wcr2wunfoldnc(wcrdvf,vmax,wcrwind)",/info
    print,"            res=wcr2wunfoldnc(wcrdvf,vmax,wind,wcrbeam)"
    return,''
  endif else if n_params() eq 3 then begin
    wcracw=wind
  endif else if n_params() eq 4 then begin
  ; Check wind input and calculate wind component into the radar beam
    if n_elements(wind) eq 3 then begin ; Use mean wind for all profiles
      acwind=replicarr(wind,np,/dim)
      for k=0L,np-1 do wcracw[k]=wcrbeam[*,k]##acwind[k,*] 
    endif else begin
      for k=0L,np-1 do wcracw[k]=wcrbeam[*,k]##wind[k,*] 
    endelse
  endif
    
; Process the shift/unfolding

  ; Form wind/shift_velocity field to include all WCR range gates. 
  ; Take care of sign convention    
    w=replicarr(-wcracw,nrg,dim=1)

  ; Scale the wind along the radar range if given

  if n_elements(windscale) eq nrg then w=w*replicarr(windscale,np)

  ; Save the wind component directions (-1,0,1); 0 means zero shift
    sw=sign(w)

  ; Do some housework and calculate the foldings
    wcrdvu= wcrdvf
    vmax2 = 2.*vmax
    nvm   = fix(w/vmax2)*vmax2

  ; Calculate the folding border
    f=w-sw*(vmax+nvm)

  ; Find all velocity points left and right of the border and shift them
  ; accordingly
    w  = sw*wcrdvf
    f  = sw*f
    ind= where(w gt f)
    if ind[0] ne -1 then wcrdvu[ind]=wcrdvf[ind]+nvm 
    ind=where(w le f)
    if ind[0] ne -1 then wcrdvu[ind]=wcrdvf[ind]+nvm+sw*vmax2

; Done

  return,wcrdvu

end 
