pro matchdata
  hh0 = loadxdr('i_al_af_1206_hh_corr.dat')
  hh0d = size(hh0,/dimensions)  
  hh1 = loadxdr('i_1206_hh.dat')
  hh1d = size(hh0,/dimensions)
  scale = 0.1
  hh0 = congrid(hh0,hh0d[0]*scale,hh0d[1]*scale)
  hh1 = congrid(hh1,hh1d[0]*scale,hh1d[1]*scale)
    
  newton([x,y,scale],'matcharr',[0,0,1])
  
  
    ; shift images and crop range which exists only in one image
    rg_shift = arrayindices(0)
    az_shift = -arrayindices(1)
  
  if (rg_shift lt 0) then begin
    image1 = data1[0:dim(0) + rg_shift-1,*]
    image2 = data2[-rg_shift:dim(0)-1,*]
  end else begin
    image1 = data1[rg_shift:dim(0)-1,*]
    image2 = data2[0:dim(0)-1 - rg_shift,*]
  endelse
  
  if (az_shift lt 0) then begin
    image1 = image1[*,0:dim(1) + az_shift-1]
    image2 = image2[*,-az_shift:dim(1)-1]
  end else begin
    image1 = image1[*,az_shift:dim(1)-1]
    image2 = image2[*,0:dim(1)-1 - az_shift,*]
  endelse