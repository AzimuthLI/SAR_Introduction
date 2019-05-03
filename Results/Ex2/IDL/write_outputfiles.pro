; Write TIFF files
write2tiff, data_bs,'img'+str(subex)+'_raw_' + outsuffix + '.tiff',ct=colortable
data_bs2 = bytscl(congrid(data,dim(0),dim(1)*asp))
write2tiff, data_bs2,'img'+str(subex)+'_raw2_' + outsuffix + '.tiff',ct=colortable

if focus_range then begin
  write2tiff, range_compr_bs,'img'+str(subex)+'_rg_' + outsuffix + '.tiff',ct=colortable
  range_compr_bs2 = bytscl(nlinearscale(congrid(range_compr,dim(0),dim(1)*asp),nlscalemodel))
  write2tiff, range_compr_bs2,'img'+str(subex)+'_rg2_' + outsuffix + '.tiff',ct=colortable

  write2tiff, r_compr_h_bs,'img'+str(subex)+'_rg_hamming_' + outsuffix + '.tiff',ct=colortable
  r_compr_h_bs2 = bytscl(nlinearscale(congrid(r_compr_h,dim(0),dim(1)*asp),nlscalemodel))
  write2tiff, r_compr_h_bs2,'img'+str(subex)+'_rg_hamming2_' + outsuffix + '.tiff',ct=colortable

  write2tiff, az_compr_hh_bs,'img'+str(subex)+'_rgaz_hamming_' + outsuffix + '.tiff',ct=colortable
  az_compr_hh_bs2 = bytscl(nlinearscale(congrid(az_compr_hh,dim(0),dim(1)*asp),nlscalemodel))
  write2tiff, az_compr_hh_bs2,'img'+str(subex)+'_rgaz_hamming2_' + outsuffix + '.tiff',ct=colortable
end;

write2tiff, az_compr_bs,'img'+str(subex)+'_az_' + outsuffix + '.tiff',ct=colortable
az_compr_bs2 = bytscl(nlinearscale(congrid(az_compr,dim(0),dim(1)*asp),nlscalemodel))
write2tiff, az_compr_bs2,'img'+str(subex)+'_az2_' + outsuffix + '.tiff',ct=colortable

write2tiff, az_compr_h_bs,'img'+str(subex)+'_az_hamming_' + outsuffix + '.tiff',ct=colortable
az_compr_h_bs2 = bytscl(nlinearscale(congrid(az_compr_h,dim(0),dim(1)*asp),nlscalemodel))
write2tiff, az_compr_h_bs2,'img'+str(subex)+'_az_hamming2_' + outsuffix + '.tiff',ct=colortable

write2tiff, az_compr_mlook_bs,'img'+str(subex)+'_mlook_' + outsuffix + '.tiff',ct=colortable
az_compr_mlook2_bs = bytscl(nlinearscale(congrid(az_compr_mlook,dim(0),dim(1)*asp),nlscalemodel))
write2tiff, az_compr_mlook2_bs,'img'+str(subex)+'_mlook2_' + outsuffix + '.tiff',ct=colortable


; Write pdf files
sucess = write2pdf(img2, 'chirp2D_' + strtrim(string(subex),2) + outsuffix + '.pdf', /bitmap)
sucess = write2pdf(img7, 'chirp2D_mlook_' + strtrim(string(subex),2) + outsuffix + '.pdf', /bitmap)
