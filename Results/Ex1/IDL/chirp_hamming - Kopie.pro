pro chirp_hamming

; define positions of the plotted axes and legends
@chirp_DefinePlotPositions

; call Batchfile to create pulse
@chirp_CreatePulse

; define hamming window
alpha    = 0.54;
hammingw = fltarr(N)
hammingw[sI] = alpha + (1-alpha)*cos(2*!pi*t[sI]/tau)

; plot chirped signal
@chirp_h_plotChirpedSignal

; create complex conjungated signal -> filter signal
h  = complexarr(N)
h  = conj(reverse(s));

; weighten the filter signal using the hamming window
hh = h*hammingw;

; plot filter function
@chirp_h_plotFilterSignal

; -----------------------
; FILTERING IN TIME SPACE
; -----------------------
; calculate compressed signal cs
cs  = convol(s,h ,center=0,/edge_wrap)

; calculate hamming windowed compressed signal csh
csh = convol(s,hh,center=0,/edge_wrap)
csh = csh/max(csh)*max(cs)

; plot compressed signal
@chirp_h_plotCompressedSignal

; --------------------------------------
; now the same in the frequency domain
; -------------------------------------
s_fft = fft(s)*N ; fft of signal
h_fft = fft(conj(reverse(s)))*N ; fft of filterfunction
f = createfarray(-1/(2*tres),1/(tres*N),1/(2*tres))

; Create humming window for frequency space
hamming_fft     = complexarr(N)
fI              = WHERE((f GE (f0-bw/2)) AND (f LE (f0+bw/2)))
hamming_fft[fI] = (alpha+(1-alpha)*cos(2*!pi*(f[fI]-f0)/bw))
 
; compress chirp in frequency space
cs_fft  = s_fft*h_fft

; compress chirp in frequency space using a hamming windowed frequency spectrum
hh_fft  = h_fft*shift(hamming_fft,round(N/2)) 
csh_fft = s_fft*hh_fft
 
; plot fft of signal, (no chirp_h.. as it is exactly the same as in Exercise 1.1)
@chirp_plotFFTSignal

; plot spectrum of signal and hamming windowed signal
@chirp_h_plotFFTFilterSignal

; plot spectrum of compressed signal and compressed hamming windowed signal
@chirp_h_plotFFTCompressedSignal

; extract phases of signals.. 

f0I = WHERE(f GE (f0))
f0I = f0I[1]
phase_signal     = phunwrap(atan(imag(shift(  s_fft,round(N/2))),real(shift(  s_fft,round(N/2)))))
phase_filter     = phunwrap(atan(imag(shift(  h_fft,round(N/2))),real(shift(  h_fft,round(N/2)))))
phase_hh         = phunwrap(atan(imag(shift( hh_fft,round(N/2))),real(shift( hh_fft,round(N/2)))))
phase_compressed = phunwrap(atan(imag(shift( cs_fft,round(N/2))),real(shift( cs_fft,round(N/2)))))
phase_csh        = phunwrap(atan(imag(shift(csh_fft,round(N/2))),real(shift(csh_fft,round(N/2)))))
;plot, f/1e6, phase_signal    - phase_signal[f0I],   color='8800ff'x, Linestyle=0, xrange=[0, 40], xtitle='Frequency (MHz)', ytitle='phase of s('+Greek('omega')+')' 
;plot, f/1e6, phase_filter    - phase_filter[f0I],   color='0088ff'x, Linestyle=0, xrange=[0, 40], xtitle='Frequency (MHz)', ytitle='phase of h('+Greek('omega')+')'
;plot, f/1e6, phase_cs_fft   - phase_cs_fft[f0I],  color='0000ff'x, Linestyle=0, xrange=[0, 40], xtitle='Frequency (MHz)', ytitle='phase of g('+Greek('omega')+')'

@chirp_plotSignalphase

@chirp_h_plotFilterphase

@chirp_h_plotCompressedSignalPhase

; fft back to time domain
cs_ifft = fft(cs_fft,1)/N;
csh_ifft = fft(csh_fft,1)/N;

; print hamming windowed signal in time domain
@chirp_h_plotFFTCompressedSignalInTime

; write result to pdf-file
pdf = 'figure2.pdf'
@chirp_writePlotwindowToFile

cs_ifft_big = shift(abs(cs_ifft) GE 0.5*max(abs(cs_ifft)),round(N/2));
print, 'uncompressed pulse length: ' + string((t(findfirst(cs_ifft_big,-1)) - t(findfirst(cs_ifft_big,1)))*1e9) + ' ns'
csh_ifft_big = shift(abs(csh_ifft) GE 0.5*max(abs(csh_ifft)),round(N/2));
print, '  compressed pulse length: ' + string((t(findfirst(csh_ifft_big,-1)) - t(findfirst(csh_ifft_big,1)))*1e9) + ' ns'


stop
end

