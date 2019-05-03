;------------------------------
; Plot data
data_bs = bytscl(congrid(data,dim(0)/2,dim(1)/2))

 
  range_compr_bs = bytscl(nlinearscale(congrid(range_compr,dim(0)/2,dim(1)/2),nlscalemodel))
  img2 = image(range_compr_bs,c*t_range2,y_azimut2,ASPECT_RATIO=asp*c*(prf/fs_range)/v,$
             position=[posx[1],posy[1],posx[1]+imgw,posy[1]+imgh], $
             location=winlocation1,dimensions=windim1, $
             rgb_table=colortable, axis_style=2, /interpolate, $
             xtickvalues=maketicks(c*t_range,3),ytickvalues=maketicks(y_azimut,3),$
             xtitle='range (m)', ytitle='azimut (m)', title='range compressed')

if focus_range then begin
img1 = image(data_bs,1e6*t_range2,y_azimut2,ASPECT_RATIO=asp*1e6*(prf/fs_range)/v,$
             position=[posx[0],posy[1],posx[0]+imgw,posy[1]+imgh], identifier=win1, /CURRENT, $
             rgb_table=colortable, axis_style=2, /interpolate,$
             xtickvalues=maketicks(1e6*t_range,3),ytickvalues=maketicks(y_azimut,3),$
             xtitle='range t ($\mus$)', ytitle='azimut (m)', title='radar echo')


  r_compr_h_bs = bytscl(nlinearscale(congrid(r_compr_h,dim(0)/2,dim(1)/2),nlscalemodel))
  img3 = image(r_compr_h_bs,c*t_range2,y_azimut2,ASPECT_RATIO=asp*c*(prf/fs_range)/v,$
             position=[posx[0],posy[0],posx[0]+imgw,posy[0]+imgh], $
             rgb_table=colortable, axis_style=2, /interpolate,$
             xtickvalues=maketicks(c*t_range,3),ytickvalues=maketicks(y_azimut,3),$
             xtitle='range (m)', ytitle='azimut (m)', title='range compressed!C(with hamming window)',/CURRENT, Identifier=itool)

az_compr_hh_bs = bytscl(nlinearscale(congrid(az_compr_hh,dim(0)/2,dim(1)/2),nlscalemodel))
img6 = image(az_compr_hh_bs,c*t_range2,y_azimut2,ASPECT_RATIO=asp*c*(prf/fs_range)/v,$
             position=[posx[2],posy[0],posx[2]+imgw,posy[0]+imgh], /CURRENT, $
             rgb_table=colortable, axis_style=2, /interpolate,$
             xtickvalues=maketicks(c*t_range,3),ytickvalues=maketicks(y_azimut,3),$
             xtitle='range (m)', ytitle='azimut (m)', $
             title='range and azimut compressed!C(with hamming window in both)')
end; 

; plot log-scale image
az_compr_bs = bytscl(nlinearscale(congrid(az_compr,dim(0)/2,dim(1)/2),nlscalemodel))
img4 = image(az_compr_bs,c*t_range2,y_azimut2,ASPECT_RATIO=asp*c*(prf/fs_range)/v,$
             position=[posx[2],posy[1],posx[2]+imgw,posy[1]+imgh], /Current, $
             rgb_table=colortable, axis_style=2, /interpolate,$
             xtickvalues=maketicks(c*t_range,3),ytickvalues=maketicks(y_azimut,3),$
             xtitle='range (m)', ytitle='azimut (m)', title='range+azimut compressed')

az_compr_h_bs = bytscl(nlinearscale(congrid(az_compr_h,dim(0)/2,dim(1)/2),nlscalemodel))
img5 = image(az_compr_h_bs,c*t_range2,y_azimut2,ASPECT_RATIO=asp*c*(prf/fs_range)/v,$
             position=[posx[1],posy[0],posx[1]+imgw,posy[0]+imgh], /CURRENT, $
             rgb_table=colortable, axis_style=2, /interpolate,$
             xtickvalues=maketicks(c*t_range,3),ytickvalues=maketicks(y_azimut,3),$
             xtitle='range (m)', ytitle='azimut (m)', $
             title='range and azimut compressed!C(with hamming window in azimut)')


az_compr_mlook_bs = bytscl(nlinearscale(congrid(az_compr_mlook,dim(0)/2,dim(1)/2),nlscalemodel)) 
img7 = image(az_compr_mlook_bs,c*t_range2,y_azimut2,ASPECT_RATIO=asp*c*(prf/fs_range)/v,$
             location=[scrsize(0)/2,scrsize(1)-50],dimensions=[scrsize(0)/4,scrsize(1)/2], $
             position=[0.1, 0.1, 0.9, 0.9], identifier=win2, $
             rgb_table=colortable, axis_style=2, /interpolate,$
             xtickvalues=maketicks(c*t_range,3),ytickvalues=maketicks(y_azimut,3),$
             xtitle='range (m)', ytitle='azimut (m)', $
             title='range and azimut compressed!C(with 2x azimut multilooking)')
