pro chirp_2D
; select subexercise:
;  1 for '2D-Chirp'
;  2 for 'ers_raw'
;  3 for 'rdemo'
subex = 3
; 
; 
; load parameter and (2D-) image data
@loadparameters_and_data

if (focus_range) then begin
; RANGE FOCUSING
; --------------------
; create chirped filter function with range-length of image and puls length of tau.
; the echo (=range-length of image) is obviously longer than the emitted pulse. 
; - create array for referencen function in range of range-length of image
; - select the elements in the middle
; - put there the filter-pulse (the center of the array)
; - for the fft:
; - shift the center of the pulse on the first element
; - create filter function in frequency space
tr     = createfarray(-Tau/2,1/Fs_range,Tau/2)
h_r   = conj(reverse(exp(-j*!pi*bw/tau*tr^2)))
ref_r = complexarr(N_rg) 
ref_I = round(N_rg-N_elements(h_r))/2 +indgen(N_elements(h_r))
ref_r[ref_I] = h_r
ref_r_fft = fft(shift(ref_r,round(N_rg/2)),-1)*N_rg
; create fitting time axis for ref_r (not required, but useful for plotting)
ref_rt = (findgen(N_rg)-round(N_rg/2))/Fs_range

;stop


;p = plot(real(ref_r_fft)/max(abs(ref_r_fft)),title='ref_range')
;p = plot(shift(hamming_r,dim(0)/2),/OVERPLOT,color='red')
;p = plot(h_r,/OVERPLOT,color='green')
;p = plot(ref_r,/OVERPLOT,color='blue')


; compress range line by line
range_compr = complexarr(dim)
for i=0, dim(1)-1 do range_compr[*,i] = fft(fft(data[*,i],-1)*ref_r_fft,1)

; RANGE FOCUSING WITH HAMMING FILTER
;--------------------------------------------
; define hamming window of same size as ref_r
alpha     = 0.54
; the following line works for time-space.
;h_r       = alpha + (1-alpha)*cos(2*!pi*tr/Tau)
; as we want to apply the hamming filter in frequency space, we need a frequency axis
; which fits to ref_r_fft. 
ham_r       = alpha + (1-alpha)*cos(2*!pi*tr/Tau)
hamming_r = fltarr(dim(0))

;hamming_r[ref_I] = shift(h_r,dim(0)/2)
; changed the line above to the following two:
hamming_r[ref_I] = ham_r
hamming_r = shift(hamming_r,dim(0)/2)
ref_r_fft_hamming = fft(hamming_r*shift(ref_r,round(N_rg/2)),-1)*N_rg


; compress in range with hamming windowing
r_compr_h = complexarr(dim)
for i=0, dim(1)-1 do r_compr_h[*,i] = fft(fft(data[*,i],-1)*ref_r_fft_hamming,1)
end else begin; end of range focusing
  range_compr = data 
end

;p = plot(real(ref_r_fft)/max(abs(ref_r_fft)),title='ref_range')
;p = plot(hamming_r,/OVERPLOT,color='red')
;stop


; Azimut Compression of range compressed data
;-----------------------------------------
l_sa  = lambda*R0/L_antenna   ; Diffraction limited length of synthetic antenna (-> numeric aperture)
k_az  = 2*!pi*v^2/(lambda*R0) ; define chirp rate (an 2. order approximation of the curvature of the hyperbolic shaped distance to the object) => quadratic phase behaviour  
t_max = l_sa/v                ; maximum time, of illumination of the object by the footprint
t_az = createfarray(-t_max/2,1.0/PRF,t_max/2)

;!!! IS THE NEXT LINE RIGHT?
hf_az   = conj(reverse(exp(-j*k_az*t_az^2))) ; create chirped azimut filter function
ref_az  = complexarr(dim(1))
ref_azI = (dim(1)-N_elements(hf_az))/2 +indgen(N_elements(hf_az))
ref_az[ref_azI] = hf_az       ; fit hf_az into the center of reference_az-function
; shift from the center to the first element to avoid frequency-shifts. (first-element -> time=zero)
ref_az_fft = fft(shift(ref_az,round(n_elements(ref_az)/2)),-1)*dim(1)


az_compr = complexarr(dim)
for i=0, dim(0)-1 do az_compr[i,*] = fft(fft(range_compr[i,*],-1)*ref_az_fft,1)

;p = plot(shift(real(ref_az_fft)/max(abs(ref_r_fft)),dim(0)/2),title='ref_azimut')



; Azimut Compression of range compressed data with hamming window
; ---------------------------------------------------------------
; define hamming window of same size as ref_az in frequency space (window-size=azimuth bandwidth)
alpha      = 0.54
h_az       = alpha + (1-alpha)*cos(2*!pi*t_az/t_max)
hamming_az = fltarr(dim(1))
hamming_az[ref_azI] = h_az
hamming_az = shift(hamming_az,n_elements(hamming_az)/2)

ref_az_fft_hamming = fft(hamming_az*shift(ref_az,round(n_elements(ref_az)/2)),-1)*dim(1)

; Compression of range compressed data with hamming filter in azimut direction
az_compr_h = complexarr(dim)
for i=0, dim(0)-1 do az_compr_h[i,*] = fft(fft(range_compr[i,*],-1)*ref_az_fft_hamming,1)


if (focus_range) then begin
; Compression of hamming filtered range compressed data with hamming filter in azimut direction
; ---------------------------------------------------------------
az_compr_hh = complexarr(dim)
for i=0, dim(0)-1 do az_compr_hh[i,*] = fft(fft(r_compr_h[i,*],-1)*ref_az_fft_hamming,1)
end

;---------------------------
; Multilooking
; Create ref_puls functions to split azimut range into two halfs
Naz     = N_elements(hf_az)
ref_az1 = complexarr(dim(1))
ref_az2 = complexarr(dim(1))

; split signal, reverse each part, conjungate it an put it in the right part of the output array so that it is shifted according to what fft likes.
ref_az1[0:(Naz-1)/2]        = hf_az[Naz/2:*]
ref_az2[dim(1) - Naz/2:*] = hf_az[0:(Naz-2)/2]

ref_az1_fft = fft(ref_az1,-1)*dim(1)
ref_az2_fft = fft(ref_az2,-1)*dim(1)

az_compr_1  = fltarr(dim)
az_compr_2  = fltarr(dim)

for i=0, dim(0)-1 do begin
  az_compr_1[i,*] = abs(fft(fft(range_compr[i,*],-1)*ref_az1_fft,1))
  az_compr_2[i,*] = abs(fft(fft(range_compr[i,*],-1)*ref_az2_fft,1))  
endfor
;az_compr_mlook 
az_compr_mlook = az_compr_1 + az_compr_2



; plot images
@plot_images
stop
; write data into tiff and pdf files
@write_outputfiles

end