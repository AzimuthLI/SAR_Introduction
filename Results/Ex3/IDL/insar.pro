pro insar
  @loadparameters_and_data
  
  ; calculate cross correlation function to shift the second image axactly on the first one.
  ; The cross correlation functions is most efficiently done in frequency space via the fouriertransform of the crosscorrelationspectrum  
  @xcorrelate_and_shiftimages
  
  ; calculate interferogramm
  Interf = image1*conj(image2)
  
  ; find main frequency components in interferogramm
  ; with the "rough" assumption, that the area is quite flat.
  ; is aetna flat?!. In average, maybe yes...)  
  
  idim = size(Interf,/dimensions)
  Interf_fft = shift(fft(Interf),idim(0)/2,idim(1)/2)
  
  ; find indices which indicate the maximum of the interferogramm
  ; The max(abs()) of the interferrogram shows the maximum occurence 
  ; of spectral components.
  ; As the main fringe frequency is a result of the 
  ; different baselines / angle of viewling / tilted earth
  ; it shows the change in the frequency spectrum of the two image spectra.
  fringe_mainvalue = max(abs(Interf_fft))
  mx = max(abs(Interf_fft),location)
  indices = array_indices(abs(Interf_fft),location)
  
  ; calculate expectation value of the abs(interferrogram)*x. (similar to maximum.), but does unfortunately not work :( 
   
;  fmax_rg = findgen(idim(1))
;  fmax_az = findgen(idim(0))
;  for j=0,idim(1)-1 do fmax_rg[j] = total(abs(Interf_fft[*,j])*indgen(idim(0)))
;  fmaxrg = total(fmax_rg)/total(abs(Interf_fft))
;  for j=0,idim(0)-1 do fmax_az[j] = total(abs(Interf_fft[j,*])*indgen(idim(1)))
;  fmaxaz = total(fmax_az)/total(abs(Interf_fft))
  
  rg_shift1 = idim(0)/2-indices[0]
  az_shift1 = idim(1)/2- indices[1]
;  rg_shift2 = idim(0)/2- fmaxrg
;  az_shift2 = idim(1)/2- fmaxaz
  
  ; shift image 2 in frequency spectrum to remove main fringe frequency 
  img2_flat_fft1 = shift( fft(image2,-1) ,-rg_shift1, -az_shift1)  
  image2_flat1 = fft(img2_flat_fft1,1)
  
  ; calculate interferogram
  Interf_flat1 = image1*conj(image2_flat1)
  
  ; calculate frequency spectrum of the two images
  rg_spec1 = shift(total(abs(fft(image1))      ,2), idim[0]/2)
  rg_spec2 = shift(total(abs(fft(image2_flat1)),2), idim[0]/2)
  rg_spec2a = shift(total(abs(fft(image2)),2), idim[0]/2)
  
  ; create 'frequency' axis, actually just the pixel indices, with pixel zero in the middle
  xpixels = indgen(idim[0])-idim[0]/2
  ypixels = indgen(idim[1])-idim[1]/2
  
  ; calculate bandwidth of images in number of pixels
  signalbw1 = fwhm(xpixels,rg_spec1,POS = fwhm1x, YVAL = fwhm1y)
  signalbw2 = fwhm(xpixels,rg_spec1,POS = fwhm2x, YVAL = fwhm2y)
  
  ; the frequency "where the common bands coincide" is given by the fringe frequency.
  ; see GATELLI_94
  rg_deltaf = rg_shift1
  
  ; create range filters. When using the fringe frequency, then is the 
  ; flat-earth-removal already included.
  
  ; calculate filter bandwidth 
  filterbw = (signalbw1+signalbw2)/2-abs(rg_deltaf)
  ; calculate filter central frequencies
  filter1fc = -rg_deltaf/2
  filter2fc = +rg_deltaf/2 
  filter1 = fltarr(idim[0])
  filter2 = fltarr(idim[0])  
  filter1[where(xpixels gt filter1fc-filterbw/2 and xpixels lt filter1fc+filterbw/2),*] = 1
  filter2[where(xpixels gt filter2fc-filterbw/2 and xpixels lt filter2fc+filterbw/2),*] = 1  
  
  img1filtered = complexarr(dim(image1))
  img2filtered = complexarr(dim(image2))
  for j=0, idim[1]-1 do img1filtered[*,j] = fft(fft(image1[*,j],-1)*shift(filter1,ddim[0]/2),1)  
  for j=0, idim[1]-1 do img2filtered[*,j] = fft(fft(image2[*,j],-1)*shift(filter2,ddim[0]/2),1)
  
  ; shift image 2 in frequency spectrum to remove main fringe frequency 
  img2_flat_fil    = fft(shift( fft(img2filtered,-1) ,-rg_shift1, -az_shift1),1)    
  Interf_flat_fil = img1filtered*conj(img2_flat_fil)
  

  ; coherence calculation
  w_size = 5
  
  ; coherence of images which are just shifted
  coh12  = smooth( image1*conj(image2) ,w_size,/edge_truncate)
  coh11  = image1*conj(image1)
  coh22  = image2*conj(image2)
  denum  = sqrt(smooth(coh11,w_size,/edge_truncate)*smooth(coh22,w_size,/edge_truncate))
  coherence = abs(coh12/denum)
  coh_histo = histogram(abs(coherence),nbins=1000)
  
  ; coherence of image which are shifted and filtered
  coh12f  = smooth( img1filtered*conj(img2filtered) ,w_size,/edge_truncate)
  coh11f  = img1filtered*conj(img1filtered)
  coh22f  = img2filtered*conj(img2filtered)
  denumf  = sqrt(smooth(coh11f,w_size,/edge_truncate)*smooth(coh22f,w_size,/edge_truncate))
  coherencef = abs(coh12f/denumf)
  cohf_histo = histogram(abs(coherencef),nbins=1000)
 
  
  @plotdata
  stop
end