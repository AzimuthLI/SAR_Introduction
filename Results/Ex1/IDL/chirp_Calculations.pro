
; create complex conjungated signal, called "filter signal"
h = complexarr(N)
h = conj(reverse(s));

; FILTERING: calculate compressed signal cs(t)
cs = shift(convol(s,h,center=0,/edge_wrap),round(N/2))

; --------------------------------------
; now the same in the frequency domain
; -------------------------------------
s_fft = fft(s)*N ; fft of signal
h_fft = fft(h)*N ; fft of filterfunction

; FILTERING: compress chirp in frequency space
cs_fft = s_fft*h_fft
 
; calculate the signal phase(f)
phase_signal =  phunwrap(atan(imag(shift(s_fft,round(N/2))),real(shift(s_fft,round(N/2)))))

; calculate phase of the signal
phase_filter = phunwrap(atan(imag(shift(h_fft,round(N/2))),real(shift(h_fft,round(N/2)))))

; calculate phase of the compressed signal
phase_compressed = phunwrap(atan(imag(shift(cs_fft,round(N/2))),real(shift(cs_fft,round(N/2)))))

; fouriertransform the frequency space compressed signal back to time-space
cs_ifft = fft(cs_fft,1)/N;

; calculate length of compressed pulse
cs_ifft_big = abs(cs_ifft) GE 0.5*max(abs(cs_ifft));
cs_tau = (t(findfirst(cs_ifft_big,/last)) - t(findfirst(cs_ifft_big))); 
