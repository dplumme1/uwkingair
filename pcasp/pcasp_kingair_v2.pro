;
pro pcasp_kingair_v2
;
common plotting, xsize_inch, ysize_inch, width_inch, width_normal, height_inch, $
 height_normal, height_spacer_inch, height_spacer_normal, width_spacer_inch, width_spacer_normal, $
 yoffset_inch, yoffset_normal, xoffset_inch, xoffset_normal
;
common normal_coordinates, x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal
;
nc_filename       = '20150506.c1.nc'
nc_path_filename  = 'C:\jeff\data_kingair\' + nc_filename
;hhmmsss_size_dist = [-1.]
;hhmmsse_size_dist = [-1.]
hhmmsss_size_dist = [165000., 174000., 173100., 175900.]
hhmmsse_size_dist = [170000., 174500., 173300., 180300.]
tcnts_size_dist   = convert_hhmmss_to_seconds(hhmmsss_size_dist)
tcnte_size_dist   = convert_hhmmss_to_seconds(hhmmsse_size_dist)
;
ps_filename_frag = (strsplit(nc_filename,'.',/extract))[0]
ps_path_filename = 'C:\jeff\ccn_kapee\out\pcasp_kingair_v2' + ps_filename_frag + '_.ps'
start_ps_plot, ps_path_filename
;
iplot                 = 0
;
width_inch            =  7.00
width_normal          = width_inch / xsize_inch
height_inch           =  1.
height_normal         = height_inch / ysize_inch
height_spacer_inch    =  0.25
height_spacer_normal  = height_spacer_inch / ysize_inch
width_spacer_inch     =  1.5
width_spacer_normal   = width_spacer_inch / xsize_inch
yoffset_inch          =  0.3
yoffset_normal        = yoffset_inch / ysize_inch
xoffset_inch          =  1.0
xoffset_normal        = xoffset_inch / xsize_inch
;
make_tcnt_kingair, nc_path_filename, tcnt
;
get_nc_variable, 'tas', nc_path_filename, tas, nsec
get_nc_variable, 'ztrue', nc_path_filename, altitude, nsec
get_nc_variable, 'conc_cpc', nc_path_filename, cpc_conc, nsec
get_nc_variable, 'lwc100', nc_path_filename, lwc100, nsec
;
get_spectra, 'AS200_OBR', nc_path_filename, spectra, nsec, nchan ; raw counts
cnts         = spectra[*,1:*] ; see NetCDF header
result       = size(cnts)
nchan        = result[2]
index_true   = where(cnts ge 0, index_count_true, complement=index_false, ncomplement=index_count_false)
if index_count_false gt 0 then begin  ; creat NAN where count data is zero or smaller
 index2d = array_indices(cnts, index_false)
 cnts[index2d[0,*],index2d[1,*]] = !VALUES.F_NAN
endif
;
;..sizing
;
my_var_name = 'AS200_OBR'
my_attrib   = 'CellSizes'
get_attrib, nc_path_filename, my_var_name, my_attrib, adia
print, adia
adia_lower  = adia(0:29) ; see NetCDF header
adia_upper  = adia(1:30) ; see NetCDF header
adia_mid    = (adia_lower + adia_upper) / 2.
adia_dlog10d = alog10(adia_upper / adia_lower)
;
if nchan ne n_elements(adia_lower) then stop
;
;..flows
;
get_nc_variable, 'PFLWC_OBR', nc_path_filename, pflwc, nsec
index_true       = where(pflwc ge 0 and pflwc le 5, index_count_true, $
                   complement=index_false, ncomplement=index_count_false)
if index_count_false gt 0 then pflwc(index_false) = !VALUES.F_NAN
get_nc_variable, 'PFLW_OBR', nc_path_filename, pflw, nsec
get_nc_variable, 'PFLWS_OBR', nc_path_filename, pflws, nsec
;
;..total counts channels 2 to 29
;
sum_cnts_2_29    = total(cnts(*,2:29), [2])
;
;..total concentration channels 0 to 29
;
pcasp_conc_total = total(cnts(*,0:29), [2]) / pflwc
;
xrange = [tcnt[0], tcnt[-1]]
;
iplot = iplot + 1
set_p_position_yt, iplot
xyouts, x_coord_min_normal+0.01, y_coord_max_normal+0.005, 'NOAA PCASP  '+nc_path_filename, /normal
plot, /noerase, /nodata, fltarr(1), ytitle = 'TAS, m/s', yrange=[0., 150.], $
 xrange = xrange, xstyle=5, ystyle = 1, yticks = 3, $
 yminor = 1, yticklen = -0.01
oplot, tcnt, tas
timeaxis, tcnt, nticks=3, charsize=1, nolabels=1, ticklen = -10
;
iplot = iplot + 1
set_p_position_yt, iplot
plot, /noerase, /nodata, fltarr(1), ytitle = 'Altitude, m', yrange=[500.,6500.], $
 xrange = xrange, xstyle=5, ystyle = 1, yticks = 2, $
 yminor = 3, yticklen = -0.01
oplot, tcnt, altitude
timeaxis, tcnt, nticks=3, charsize=1, nolabels=1, ticklen = -10
;
iplot = iplot + 1
set_p_position_yt, iplot
plot, /noerase, /nodata, fltarr(1), ytitle = 'Flows', yrange=[0.0,2.0], $
 xrange = xrange, xstyle=5, ystyle = 1, yticks = 2, $
 yminor = 2, yticklen = -0.01
oplot, tcnt, pflw, linestyle = 0
oplot, tcnt, pflws / 10., linestyle = 1
timeaxis, tcnt, nticks=3, charsize=1, nolabels=1, ticklen = -10
;
iplot = iplot + 1
set_p_position_yt, iplot
plot, /noerase, /nodata, fltarr(1), ytitle = 'Chan 0, count', yrange=[-100., 900.], $
 xrange = xrange, xstyle=5, ystyle = 1, yticks = 5, $
 yminor = 2, yticklen = -0.01
oplot, tcnt, cnts(*,0)
timeaxis, tcnt, nticks=3, charsize=1, nolabels=1, ticklen = -10
;
iplot = iplot + 1
set_p_position_yt, iplot
plot, /noerase, /nodata, fltarr(1), ytitle = 'Chan 1, count', yrange=[-100., 900.], $
 xrange = xrange, xstyle=5, ystyle = 1, yticks = 5, $
 yminor = 2, yticklen = -0.01
oplot, tcnt, cnts(*,1)
timeaxis, tcnt, nticks=3, charsize=1, nolabels=1, ticklen = -10
;
iplot = iplot + 1
set_p_position_yt, iplot
plot, /noerase, /nodata, fltarr(1), ytitle = 'Chan 2-29, count', yrange=[-100., 900.], $
 xrange = xrange, xstyle=5, ystyle = 1, yticks = 5, $
 yminor = 2, yticklen = -0.01
oplot, tcnt, sum_cnts_2_29
timeaxis, tcnt, nticks=3, charsize=1, nolabels=1, ticklen = -10
;
iplot = iplot + 1
set_p_position_yt, iplot
plot, /noerase, /nodata, fltarr(1), ytitle = 'PCASP & CPC, cm!U-3!N', yrange=[0., 3000.], $
 xrange = xrange, xstyle=5, ystyle = 1, yticks = 3, $
 yminor = 1, yticklen = -0.01
oplot, tcnt, pcasp_conc_total, linestyle = 0
oplot, tcnt, cpc_conc, linestyle = 1
timeaxis, tcnt, nticks=3, charsize=1, nolabels=1, ticklen = -10
;
iplot = iplot + 1
set_p_position_yt, iplot
plot, /noerase, /nodata, fltarr(1), ytitle = 'lwc100, g m!U-3!N', yrange=[0., 1.5], $
 xrange = xrange, xstyle=5, ystyle = 1, yticks = 3, $
 yminor = 1, yticklen = -0.01
oplot, tcnt, lwc100, linestyle = 0
timeaxis, tcnt, nticks=3, charsize=1, nolabels=0, ticklen = -10, title = 'UTC'
;
;..draw in the time intervals for size distributions
;
y_coord_max_normal    = 1. - yoffset_normal
y_coord_min_normal    = y_coord_max_normal - 8.*height_normal - 7.*height_spacer_normal
!p.position           = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
;
yrange = [0., 1.]
;
plot, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = 5, xrange = xrange, xstyle = 5 
for i = 0, n_elements(tcnts_size_dist) - 1 do begin  
 oplot, [tcnts_size_dist[i], tcnts_size_dist[i]], [0., 1.], linestyle = 1
 oplot, [tcnte_size_dist[i], tcnte_size_dist[i]], [0., 1.], linestyle = 1
endfor
;
iplot = 1
;
xsize_inch            =  8.5
ysize_inch            = 11.0
width_inch            =  2.2
width_normal          = width_inch / xsize_inch
height_inch           =  2.2
height_normal         = height_inch / ysize_inch
height_spacer_inch    =  0.2
height_spacer_normal  = height_spacer_inch / ysize_inch
width_spacer_inch     =  0.2
width_spacer_normal   = width_spacer_inch / xsize_inch
yoffset_inch          =  0.75
yoffset_normal        = yoffset_inch / ysize_inch
xoffset_inch          =  1.00
xoffset_normal        = xoffset_inch / xsize_inch
;
xrange = [0.1, 10.]
xticks = 2
xminor = 9
xstyle = 1
xtitle = 'Diameter, !9m!7m'
;
yrange = [0.01, 1000.]
yticks = 5
yminor = 9
ystyle = 1
ytitle = 'dN/dlog!D10!ND, cm!U-3!N'
;
for i = 0, n_elements(tcnts_size_dist) - 1 do begin  
 if iplot            eq 1 then begin
  plot, /nodata, fltarr(1), xstyle=5, ystyle=5
  x_coord_min_normal    = xoffset_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = (ysize_inch - height_inch) / ysize_inch - 0.*height_spacer_normal - yoffset_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
   yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
   xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), $
   ytitle = ytitle
 endif else if iplot eq 2 then begin
  x_coord_min_normal    = x_coord_max_normal + width_spacer_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = y_coord_min_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), ytickname = replicate(' ', yticks + 1) 
 endif else if iplot eq 3 then begin
  x_coord_min_normal    = x_coord_max_normal + width_spacer_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = y_coord_min_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), ytickname = replicate(' ', yticks + 1) 
 endif else if iplot eq 4 then begin
  x_coord_min_normal    = xoffset_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = (ysize_inch - 2.*height_inch) / ysize_inch - 1.*height_spacer_normal - yoffset_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), $
  ytitle = ytitle
 endif else if iplot eq 5 then begin
  x_coord_min_normal    = x_coord_max_normal + width_spacer_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = y_coord_min_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), ytickname = replicate(' ', yticks + 1) 
 endif else if iplot eq 6 then begin
  x_coord_min_normal    = x_coord_max_normal + width_spacer_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = y_coord_min_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), ytickname = replicate(' ', yticks + 1) 
 endif else if iplot eq 7 then begin
  x_coord_min_normal    = xoffset_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = (ysize_inch - 3.*height_inch) / ysize_inch - 2.*height_spacer_normal - yoffset_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), $
  ytitle = ytitle
 endif else if iplot eq 8 then begin
  x_coord_min_normal    = x_coord_max_normal + width_spacer_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = y_coord_min_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), ytickname = replicate(' ', yticks + 1) 
 endif else if iplot eq 9 then begin
  x_coord_min_normal    = x_coord_max_normal + width_spacer_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = y_coord_min_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), ytickname = replicate(' ', yticks + 1) 
 endif else if iplot eq 10 then begin
  x_coord_min_normal    = xoffset_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = (ysize_inch - 4.*height_inch) / ysize_inch - 3.*height_spacer_normal - yoffset_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtitle = xtitle, ytitle = ytitle
 endif else if iplot eq 11 then begin
  x_coord_min_normal    = x_coord_max_normal + width_spacer_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = y_coord_min_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, ytickname = replicate(' ', yticks + 1), xtitle = xtitle 
 endif else if iplot eq 12 then begin
  x_coord_min_normal    = x_coord_max_normal + width_spacer_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = y_coord_min_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, ytickname = replicate(' ', yticks + 1), xtitle = xtitle 
 endif else if iplot eq 13 then begin
  iplot = 1
  plot, /xlog, /ylog, /nodata, fltarr(1), xstyle=5, ystyle=5
  x_coord_min_normal    = xoffset_normal
  x_coord_max_normal    = x_coord_min_normal + width_normal
  y_coord_min_normal    = (ysize_inch - height_inch) / ysize_inch - 0.*height_spacer_normal - yoffset_normal
  y_coord_max_normal    = y_coord_min_normal + height_normal
  !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
  plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
  yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
  xticklen = -0.02, yticklen=-0.02, xtickname = replicate(' ', xticks + 1), $
  ytitle = ytitle
 endif
;
 index = where (tcnt gt tcnts_size_dist[i] and tcnt le tcnte_size_dist[i], index_count)
 if index_count gt 0 then begin
  dndlog10d  = total(cnts[index,*], [1]) / (total(pflwc[index]) * adia_dlog10d)
;
;..from Cai et al. AMT (2013), Equation A4
;
  poisson_error  = (1/adia_dlog10d) * sqrt(total(cnts[index,*], [1]) / (total(pflwc[index]) * total(pflwc[index])))
  plots, adia_mid, dndlog10d, symsize = 1.0
  for j = 0, nchan - 1 do begin
   if dndlog10d[j] - poisson_error[j] gt yrange[0] then begin
    oplot, [adia_mid[j], adia_mid[j]], [dndlog10d[j] - poisson_error[j], dndlog10d[j] + poisson_error[j]]
   endif else if dndlog10d[j] + poisson_error[j] gt yrange[0] then begin
    oplot, [adia_mid[j], adia_mid[j]], [yrange[0], dndlog10d[j] + poisson_error[j]]
   endif
  endfor
;
  xyouts, x_coord_min_normal + 0.02, y_coord_max_normal - 0.02, $
   string(hhmmsss_size_dist[i], format='(i7)') + $
   string(hhmmsse_size_dist[i], format='(i7)'), /normal
;
  iplot = iplot + 1
;
 endif
;
endfor
;
ystyle = 5
plot, /xlog, /ylog, /noerase, /nodata, fltarr(1), yrange = yrange, ystyle = ystyle, $
 yticks = yticks, yminor = yminor, xrange = xrange, xstyle = xstyle, xticks = xticks, xminor = xminor, $
 xticklen = -0.02, yticklen=-0.02, xtitle = xtitle
;
resetps2x
;
end
;
;
pro get_attrib, nc_path_filename, my_var_name, my_attrib, my_return
 fptr   = ncdf_open(nc_path_filename, /nowrite)
 result = ncdf_inquire(fptr)
;
 for i = 0, result.nvars - 1 do begin
  result = ncdf_varinq(fptr, i)
  my_string = result.name
  for j = 0, result.natts - 1 do begin
   if strcmp(result.name, my_var_name, 9, /fold_case) then begin
    if strcmp(ncdf_attname(fptr, i, j), my_attrib, 9, /fold_case) then begin
     ncdf_attget, fptr, i, ncdf_attname(fptr, i, j), my_return
    endif
   endif
  endfor
 endfor
;
 ncdf_close, fptr
;
end
;
;
function convert_hhmmss_to_seconds, hhmmss
 hr1      = long (hhmmss / 10000.)
 mn1      = long (hhmmss - hr1 * 10000.) / 100
 sc1      = long (hhmmss - hr1 * 10000. - mn1 * 100.)
 return, hr1*long(3600) + mn1*long(60) + sc1
end
;
;
pro make_tcnt_kingair, nc_path_filename, tcnt
 fptr = ncdf_open (nc_path_filename, /nowrite)
 ncdf_varget, fptr, 'HOUR', hr
 ncdf_varget, fptr, 'MINUTE', mn
 ncdf_varget, fptr, 'SECOND', sc
 nsec_to_start = hr(0)*3600. + mn(0)*60. + sc(0)
 tcnt = nsec_to_start + findgen(n_elements(hr))
end
;
;
pro get_nc_variable, my_var_name, nc_path_filename, array, nsec
 fptr = ncdf_open (nc_path_filename, /nowrite)
 ncdf_varget, fptr, my_var_name, array
 ncdf_close, fptr
 nsec = n_elements(array)
end
;
;
pro get_spectra, my_var_name, nc_path_filename, spectra, nsec, nchan
 fptr = ncdf_open (nc_path_filename, /nowrite)
 ncdf_varget, fptr, my_var_name, spectra
 spectra = transpose (reform (spectra))
 nsec   = (size (spectra, /L64))[1]
 nchan  = (size (spectra, /L64))[2]
 ncdf_close, fptr
end
;
;
pro start_ps_plot, ps_path_filename
;
common plotting, xsize_inch, ysize_inch, width_inch, width_normal, height_inch, $
 height_normal, height_spacer_inch, height_spacer_normal, width_spacer_inch, width_spacer_normal, $
 yoffset_inch, yoffset_normal, xoffset_inch, xoffset_normal
;
xsize_inch            =  8.5
ysize_inch            =  11.
;
set_plot, 'ps'
!p.font    =    0
device, /times, font_size = 11, filename = ps_path_filename, /inches, yoffset = 0., $
 xoffset = 0., xsize = xsize_inch, ysize = ysize_inch, bits_per_pixel = 8, /color
;
end
;
;
pro set_p_position_yt, iplot
;
common plotting, xsize_inch, ysize_inch, width_inch, width_normal, height_inch, $
 height_normal, height_spacer_inch, height_spacer_normal, width_spacer_inch, width_spacer_normal, $
 yoffset_inch, yoffset_normal, xoffset_inch, xoffset_normal
;
common normal_coordinates, x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal
;
if iplot            eq 1 then begin
 x_coord_min_normal    = xoffset_normal
 x_coord_max_normal    = x_coord_min_normal + width_normal
 y_coord_min_normal    = (ysize_inch - 1.*height_inch) / ysize_inch - yoffset_normal
 y_coord_max_normal    = y_coord_min_normal + height_normal
 !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
endif else if iplot eq 2 then begin
 x_coord_min_normal    = xoffset_normal
 x_coord_max_normal    = x_coord_min_normal + width_normal
 y_coord_min_normal    = (ysize_inch - 2.*height_inch) / ysize_inch - yoffset_normal - 1.*height_spacer_normal
 y_coord_max_normal    = y_coord_min_normal + height_normal
 !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
endif else if iplot eq 3 then begin
 x_coord_min_normal    = xoffset_normal
 x_coord_max_normal    = x_coord_min_normal + width_normal
 y_coord_min_normal    = (ysize_inch - 3.*height_inch) / ysize_inch - yoffset_normal - 2.*height_spacer_normal
 y_coord_max_normal    = y_coord_min_normal + height_normal
 !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
endif else if iplot eq 4 then begin
 x_coord_min_normal    = xoffset_normal
 x_coord_max_normal    = x_coord_min_normal + width_normal
 y_coord_min_normal    = (ysize_inch - 4.*height_inch) / ysize_inch - yoffset_normal - 3.*height_spacer_normal
 y_coord_max_normal    = y_coord_min_normal + height_normal
 !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
endif else if iplot eq 5 then begin
 x_coord_min_normal    = xoffset_normal
 x_coord_max_normal    = x_coord_min_normal + width_normal
 y_coord_min_normal    = (ysize_inch - 5.*height_inch) / ysize_inch - yoffset_normal - 4.*height_spacer_normal
 y_coord_max_normal    = y_coord_min_normal + height_normal
 !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
endif else if iplot eq 6 then begin
 x_coord_min_normal    = xoffset_normal
 x_coord_max_normal    = x_coord_min_normal + width_normal
 y_coord_min_normal    = (ysize_inch - 6.*height_inch) / ysize_inch - yoffset_normal - 5.*height_spacer_normal
 y_coord_max_normal    = y_coord_min_normal + height_normal
 !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
endif else if iplot eq 7 then begin
 x_coord_min_normal    = xoffset_normal
 x_coord_max_normal    = x_coord_min_normal + width_normal
 y_coord_min_normal    = (ysize_inch - 7.*height_inch) / ysize_inch - yoffset_normal - 6.*height_spacer_normal
 y_coord_max_normal    = y_coord_min_normal + height_normal
 !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
endif else if iplot eq 8 then begin
 x_coord_min_normal    = xoffset_normal
 x_coord_max_normal    = x_coord_min_normal + width_normal
 y_coord_min_normal    = (ysize_inch - 8.*height_inch) / ysize_inch - yoffset_normal - 7.*height_spacer_normal
 y_coord_max_normal    = y_coord_min_normal + height_normal
 !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
endif else if iplot eq 9 then begin
 iplot = 1
 plot, /nodata, fltarr(1), xstyle=5, ystyle=5
 x_coord_min_normal    = xoffset_normal
 x_coord_max_normal    = x_coord_min_normal + width_normal
 y_coord_min_normal    = (ysize_inch - height_inch) / ysize_inch - yoffset_normal
 y_coord_max_normal    = y_coord_min_normal + height_normal
 !p.position   = [x_coord_min_normal,y_coord_min_normal,x_coord_max_normal,y_coord_max_normal]
endif
;
end