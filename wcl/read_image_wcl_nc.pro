;*************************************************************************************
; Description:
;
; For reading and imaging the 2 channels Upward WCL netcdf fles

; Version: 1.0.0
; 
;Author: Min Deng
;*************************************************************************************
;
pro read_image_wcl_nc, path,fname

common wcl_data_p, p1r2,p2r2, pch,sch,dep
common wcl_data_data1, galt,trf,pmb, lon, lat, zenith,Time,h
common wcl_data_data2, pchbg,schbg,pchbgstd,schbgstd, pchsat, schsat

 

gif_dir=''
;path='/home/mdeng2/wcl_process/quicklooks/'
;fname='150522WCL.UP.pecan15.20150522.035000-035960.nc'
;fname='150522WCL.UP.pecan15.20150522.030000-030960.nc'
  cdfid=ncdf_open(path+fname)   ; open the netcdf file
  
   
   id_profile  = ncdf_dimid(cdfid,'profile' )
   id_Datarange  = ncdf_dimid(cdfid,'range' )
   id_vector3  = ncdf_dimid(cdfid,'vector3' )
  
  ncdf_diminq, cdfid, id_profile, char_strng, nray
  ncdf_diminq, cdfid,  id_Datarange, char_strng, nbin
  ncdf_diminq, cdfid,  id_vector3, char_strng, nvector

  ; get the id's of the variables to be read    
  Time_id=ncdf_varid(cdfid,'Time')
  timeSec_id=ncdf_varid(cdfid,'time')
  range_id=ncdf_varid(cdfid,'Range')  
  Pch_id=ncdf_varid(cdfid,'CopolPower')
  Sch_id=ncdf_varid(cdfid,'CrossPower')
  Pch_overlap_id=ncdf_varid(cdfid,'CopolOverlap')
  Sch_overlap_id=ncdf_varid(cdfid,'CrossOverlap')
  PchBG_id=ncdf_varid(cdfid,'CopolBG')
  SchBG_id=ncdf_varid(cdfid,'CrossBG')
  PchBGstd_id=ncdf_varid(cdfid,'CopolBGSTD')
  SchBGstd_id=ncdf_varid(cdfid,'CrossBGSTD')
  Pch_st_id=ncdf_varid(cdfid,'CopolSatur')
  Sch_st_id=ncdf_varid(cdfid,'CrossSatur')
  Badfl_id=ncdf_varid(cdfid,'Prof_qc_flag')  
  Galt_id=ncdf_varid(cdfid,'ALT')
  Ralt_id=ncdf_varid(cdfid,'Ralt')
  Pitch_id=ncdf_varid(cdfid,'Pitch')
  Roll_id=ncdf_varid(cdfid,'Roll')
  lon_id=ncdf_varid(cdfid,'LON')
  lat_id=ncdf_varid(cdfid,'LAT')
  trf_id=ncdf_varid(cdfid,'trf')
  pmb_id=ncdf_varid(cdfid,'pmb')
  zenith_id=ncdf_varid(cdfid,'Zenith')
  attitude_id=ncdf_varid(cdfid,'BeamVector')

  ncdf_varget, cdfid, Time_id,Time
  ncdf_varget, cdfid, timeSec_id,time_sec
  ncdf_varget, cdfid, range_id,h
  ncdf_varget, cdfid, pch_id, pch
  ncdf_varget, cdfid, sch_id, sch
  ncdf_varget, cdfid, pchbg_id, pchbg
  ncdf_varget, cdfid, schbg_id, schbg
  ncdf_varget, cdfid, pchbgstd_id, pchbgstd
  ncdf_varget, cdfid, schbgstd_id, schbgstd
  ncdf_varget, cdfid, pch_st_id, pchsat
  ncdf_varget, cdfid, sch_st_id, schsat
  ncdf_varget, cdfid, badfl_id, badfl
  ncdf_varget, cdfid, galt_id, galt
  ncdf_varget, cdfid, ralt_id, ralt
  ncdf_varget, cdfid, lon_id, lon
  ncdf_varget, cdfid, lat_id, lat
  ncdf_varget, cdfid, trf_id, trf
  ncdf_varget, cdfid, pmb_id, pmb
  ncdf_varget, cdfid, zenith_id, zenith
  ncdf_varget, cdfid, attitude_id,Beamvector
   
  ncdf_close, cdfid
  
   
   
   
   
   p1r2=fltarr(nray,nbin)
   p2r2=fltarr(nray,nbin)
   dep=fltarr(nray,nbin)
   for i=0, nray-1 do begin
   p1r2[i,*] = 1.*(pch[i,*])*((h/1000.)^2.)
   p2r2[i,*] = 1.*(sch[i,*])*((h/1000.)^2.)
   endfor
   
   sm=[11,33]
   p1r2_sm = smooth(p1r2,sm)
   p2r2_sm = smooth(p2r2,sm)   
   dep   = p2r2_sm/p1r2_sm/3.67589/4.5  ; calibration constant, needs to calibrated after each alignment
   help, p1r2,p2r2, pch,sch,dep
   
   
    
   
  name=FILE_BASENAME(fname )
  name=strmid(name,0,strlen(name)-3)
  plot_wcl_power_dep, p1r2_sm, p2r2_sm,galt,Time,h,gif_dir,name
   
  
  end
  
