
  data1_fft = fft(abs(data1))
  data2_fft = fft(abs(data2))
  corr_fft  = data1_fft*conj(data2_fft)
  corr      = fft(corr_fft,-1)
  
  ; find coordinates of maximum of correlation function
  corr_max     = max(abs(corr),location)
  arrayindices = array_indices(corr,location)
   
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