pro chirp_multilook

; define positions of the plotted axes and legends
@chirp_DefinePlotPositions

; call Batchfile to create pulse
@chirp_CreatePulse

; do the simple Pulse calculations
@chirp_Calculations
;@chirp_plotalldata

; ----------------------------------
; Multilooking
; ----------------------------------
; Multilooking in time space
; split signal
s1 = s*(t LE 0) 
s2 = s*(t GT 0)

; split filter function (actually not required. To split s(t) OR h(t) is enough)
h1 = conj(reverse(s1))
h2 = conj(reverse(s2)) 

; calculate autocorrelation function
cs1 = shift(convol(s1,h1,center=0,/edge_wrap),round(N2))
cs2 = shift(convol(s2,h2,center=0,/edge_wrap),round(N2))

; Multilooking in frequency space
; split redundant frequency spectrum in two parts

; -> mix spectrum down to the base band, by shifting the frequency spectrum by f0.
s_ffts = shift(s_fft,-f0/max(f)*N2)
h_ffts = shift(h_fft,-f0/max(f)*N2)

; split spectrum and filter function
s1_fft = s_ffts*(f LE 0)
s2_fft = s_ffts*(f GT 0)
h1_fft = h_ffts*(f LE 0)
h2_fft = h_ffts*(f GT 0)

; compress signal in frequency space
cs1_fft = s1_fft*h1_fft
cs2_fft = s2_fft*h2_fft 

; shift baseband back to carrier and transform back to time domain
cs_mlook_ifft= shift(abs(fft(shift(cs1_fft,N2),1))+abs(fft(shift(cs2_fft,N2),1)),N2)/N

; shift spectra and calculate phases
s1_fft  = shift( s1_fft,N2)
s2_fft  = shift( s2_fft,N2)
h1_fft  = shift( h1_fft,N2)
h2_fft  = shift( h2_fft,N2)
cs1_fft = shift(cs1_fft,N2)
cs2_fft = shift(cs2_fft,N2)

; calculate phases
phase_s1_fft   = phunwrap(atan(imag(s1_fft)  , real(s1_fft)))
phase_s2_fft   = phunwrap(atan(imag(s2_fft)  , real(s2_fft)))
phase_h1_fft   = phunwrap(atan(imag(h1_fft)  , real(h1_fft)))
phase_h2_fft   = phunwrap(atan(imag(h2_fft)  , real(h2_fft)))
phase_cs1_fft   = phunwrap(atan(imag(cs1_fft)  , real(cs1_fft)))
phase_cs2_fft   = phunwrap(atan(imag(cs2_fft)  , real(cs2_fft)))


@chirp_m_plotSignal
@chirp_m_plotFilter
@chirp_m_plotCompressedSignal
@chirp_m_plotFFTSignal
@chirp_m_plotFFTFilter 
@chirp_m_plotFFTCompressedSignal
@chirp_m_plotSignalPhase
@chirp_m_plotFilterPhase
@chirp_m_plotCompressedSignalPhase

; write result to pdf-file
pdf = 'figure3.pdf'
@chirp_writePlotwindowToFile

stop
end
