; define hamming window in time domain
alpha    = 0.54;
hammingw = fltarr(N)
hammingw[sI] = alpha + (1-alpha)*cos(2*!pi*t[sI]/tau)

; weighten the filter signal using the hamming window
hh = h*hammingw;

; calculate hamming windowed compressed signal csh
csh = convol(s,hh,center=0,/edge_wrap)
csh = csh/max(csh)*max(cs)

; Create hamming window for frequency space
hamming_fft     = complexarr(N)
fI              = WHERE((f GE (f0-bw/2)) AND (f LE (f0+bw/2)))
hamming_fft[fI] = (alpha+(1-alpha)*cos(2*!pi*(f[fI]-f0)/bw))
 
; compress chirp in frequency space using a hamming windowed frequency spectrum
hh_fft  = h_fft*shift(hamming_fft,round(N/2)) 
csh_fft = s_fft*hh_fft

; extract phase of hamming window signal
phase_hh         = phunwrap(atan(imag(shift( hh_fft,round(N/2))),real(shift( hh_fft,round(N/2)))))
phase_csh        = phunwrap(atan(imag(shift(csh_fft,round(N/2))),real(shift(csh_fft,round(N/2)))))

; fft back to time domain
csh_ifft = fft(csh_fft,1)/N;
