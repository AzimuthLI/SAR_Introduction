pro chirp_multilook

; maybe, to keep the graphics nice..
device, retain=2, decomposed=39, SET_CHARACTER_SIZE=[16,16] 

; define time axis
t0 = -1e-6  ; (sec) from
tend = 1e-6 ; (sec) to
tres = 3e-10 ; (sec) resolution
t = createfarray(t0,tres,tend)
N = N_elements(t)
N2 = round(N/2)

; define chirped pulse
f0  = 20e6 ; (1/s) center frequency
bw  = 20e6; (1/s) bandwidth of linear chirped pulse
tau =  1.5e-6 ; (sec) puls length, with puls center at t = 0
c = 2.998e8;
w0 = f0*2*!pi;

; define hamming window
alpha = 0.54;
hammingw = alpha + (1-alpha)*cos(2*!pi*t/tau)

; define constants
i = complex(0,1)

; create pulse
s = exp(i*t*(w0 + !pi*bw*t/tau))*((t GE -tau/2) AND (t LE tau/2))

;create frequency axis
f = createfarray(-1/(2*tres),1/(tres*N),1/(2*tres))

; create complex conjungated signal
h  = conj(reverse(s));
hh = h*hammingw;

; plot chirped signal
!P.MULTI = [0, 3, 3] 
window, 0, XSIZE=1024,YSIZE=768
plot,  t*1e9 , real(s), color='ff0000'x, Background = 'ffffff'x, ytitle='signal s(t)', xtitle= 'time (ns)'
oplot, t*1e9 , imag(s), color='00ff00'x, linestyle=0
oplot, t*1e9, hammingw, color='509090'x, linestyle=0 

; plot filter function
plot,  t*1e9 , real(h), color='ff0088'x, ytitle='filter h(t)', xtitle= 'time (ns)'
oplot, t*1e9 , imag(h), color='00ff88'x, linestyle=1
oplot, t*1e9 , real(hh), color='0000ff'x, linestyle=0
oplot, t*1e9 , imag(hh), color='0000ff'x, linestyle=1


; filtering

; calculate compressed signal cs with autocorrelation
cs  = shift(convol(s,h,center=0,/edge_wrap),round(N2))
csh = shift(convol(s,hh,center=0,/edge_wrap),round(N2))
csh = csh/max(csh)*max(cs)

; plot compressed signal
 plot, t*1e9, real(cs), color='ff0000'x, yrange=[min(real(cs)),max(real(cs))], xrange=[-2e2, 2e2], ytitle='correlation g_1(t)', xtitle='time \tau (ns)'
oplot, t*1e9, imag(cs), color='00ff00'x
oplot, t*1e9, abs(cs), color='0000ff'x, thick=2
oplot, t*1e9, abs(csh), color='ff00aa'x, thick=2
 
;; iplot, t, real(shift(cs,round(N2))), imag(shift(cs,round(N2))), xrange = [-1e-7, 1e-7], color= [0, 64, 128], xtitle = "time (s)", ytitle = "real part", ztitle="imaginary part", thick=2  
;; iplot, t, fltarr(N) + max(real(cs)), real(shift(cs,round(N2))), color=[0, 0, 255], /OVERPLOT
;; iplot, t, imag(shift(cs,round(N2))), fltarr(N) + min(imag(cs)), color=[0, 0, 255], /OVERPLOT
; plot power spectrum
;; iplot, t, fltarr(N) + max(real(cs)), abs(shift(cs,round(N2))^2/max(abs(cs))), color=[255, 0, 0], /OVERPLOT, thick=2

; --------------------------------------
; now the same in the frequency domain
; -------------------------------------

s_fft = shift(fft(s)*N,N2) ; fft of signal
h_fft = shift(fft(conj(reverse(s)))*N,N2) ; fft of filterfunction

; Create humming window for frequency space
hamming_fft = (alpha+(1-alpha)*cos(2*!pi*(f-f0)/bw))*((f GE (f0-bw/2)) AND (f LE (f0+bw/2)))
 
; compress chirp in frequency space
cs_fft = s_fft*h_fft

; apply hamming window
hh_fft = h_fft*hamming_fft 
csh_fft = s_fft*hh_fft
 
; plot fft of signal
plot, f/1e6, abs(s_fft), xrange=[0, 40], color='8800ff'x, xtitle='Frequency (MHz)', ytitle='amplitude spectrum s('+Greek('omega')+')'
; plot fft of filter function 
plot, f/1e6, abs(h_fft), color='0088ff'x, Linestyle=1, xrange=[0, 40], xtitle='Frequency (MHz)', ytitle='amplitude spectrum h('+Greek('omega')+')'
oplot, f/1e6, hamming_fft/max(hamming_fft)*max(abs(h_fft)), color='509090'x, Linestyle=0, thick=2
oplot, f/1e6, abs(hh_fft), color='0088ff'x, Linestyle=0, thick=1

; plot fft of compressed signal
plot, f/1e6, abs(cs_fft), color='0000ff'x, Linestyle=0, xrange=[0, 40], xtitle='Frequency (MHz)', ytitle='amplitude spectrum g_1('+Greek('omega')+')'
oplot, f/1e6, abs(csh_fft), color='ff00aa'x, Linestyle=0,thick=2

f0I = findfirst(f GE (f0),1)
phase_s_fft   = phunwrap(atan(imag(s_fft)  , real(s_fft)))
phase_h_fft   = phunwrap(atan(imag(h_fft)  , real(h_fft)))
phase_hh_fft  = phunwrap(atan(imag(hh_fft) , real(hh_fft)))
phase_cs_fft  = phunwrap(atan(imag(cs_fft) , real(cs_fft)))
phase_csh_fft = phunwrap(atan(imag(csh_fft), real(csh_fft)))
; plot phase graphs
plot, f/1e6, phase_s_fft    - phase_s_fft[f0I],   color='8800ff'x, Linestyle=0, xrange=[0, 40], xtitle='Frequency (MHz)', ytitle='phase of s('+Greek('omega')+')' 
plot, f/1e6, phase_h_fft    - phase_h_fft[f0I],   color='0088ff'x, Linestyle=0, xrange=[0, 40], xtitle='Frequency (MHz)', ytitle='phase of h('+Greek('omega')+')'
oplot, f/1e6, phase_hh_fft  - phase_hh_fft[f0I],  color='509090'x, Linestyle=0, thick=2
plot, f/1e6, phase_cs_fft   - phase_cs_fft[f0I],  color='0000ff'x, Linestyle=0, xrange=[0, 40], xtitle='Frequency (MHz)', ytitle='phase of g('+Greek('omega')+')'
oplot, f/1e6, phase_csh_fft - phase_csh_fft[f0I], color='ff00aa'x, Linestyle=0, thick=2

; transform signal back to time space
cs_ifft  = shift(fft(shift(cs_fft,N2),1)/N,N2);
csh_ifft = shift(fft(shift(csh_fft,N2),1)/N,N2);
!P.MULTI = [7, 3, 3] & plot, t*1e9, abs(cs_ifft), color='000000'x,/NOERASE, yrange=[min(real(cs)),max(real(cs))],linestyle=2, xrange=[-2e2, 2e2]
oplot, t*1e9, abs(csh_ifft)/max(abs(csh_ifft))*max(cs), color='00aaff'x, linestyle=2,thick=2

; calculated length of compressed pulse
cs_ifft_big = abs(cs_ifft) GE 0.5*max(abs(cs_ifft));
cs_tau = (t(findfirst(cs_ifft_big,-1)) - t(findfirst(cs_ifft_big,1))); 
print, 'unfiltered uncompressed pulse length: ' + string(tau*1e9) + ' ns' 
print, 'unfiltered compressed pulse length: ' + string(cs_tau*1e9) + ' ns (resolution = ' + string(c*cs_tau) + 'm)' 
csh_ifft_big = abs(csh_ifft) GE 0.5*max(abs(csh_ifft));
csh_tau =  (t(findfirst(csh_ifft_big,-1)) - t(findfirst(csh_ifft_big,1)))
print, 'hamming filtered compressed pulse length: ' + string(csh_tau*1e9) + ' ns (resolution = ' + string(c*csh_tau) + 'm)'


;oplot, t, cs_ifft
; stop
; ----------------------------------
; Multilooking
; ----------------------------------
; Multilooking in time space
; split signal
s1 = s*(t LE 0) 
s2 = s*(t GT 0)

; plot splitted signal
window, 1, XSIZE=1024, YSIZE=768
!P.MULTI = [0, 3, 3] 
 plot, t*1e9, real(s1), color='00ff00'x, background='ffffff'x
oplot, t*1e9, real(s2), color='ff0000'x

h1 = conj(reverse(s1))
h2 = conj(reverse(s2)) 

; plot splitted filter functions
 plot, t*1e9, real(h1), color='00ff00'x, background='ffffff'x
oplot, t*1e9, real(h2), color='ff0000'x
 
; calculate autocorrelation function
cs1 = shift(convol(s1,h1,center=0,/edge_wrap),round(N2))
cs2 = shift(convol(s2,h2,center=0,/edge_wrap),round(N2))

 plot, t*1e9, abs(cs1), color='00ff00'x, xrange=[-2e2, 2e2], yrange=[0, max(abs(cs1))], ytitle='correlation g_1(t)', xtitle='time \tau (ns)'
oplot, t*1e9, abs(cs2), color='ff0000'x, linestyle=2
oplot, t*1e9, (abs(cs1)+abs(cs2))/2, color='0000ff'x, linestyle=2, thick=2
oplot, t*1e9, abs(cs)/max(abs(cs))*max(abs(cs1)), color='0000ff'x, linestyle=0, thick=2

; Multilooking in frequency space
; split redundant frequency spectrum in two parts

; -> mix spectrum down to the base band, by shifting the frequency spectrum by f0.
s_ffts = shift(s_fft,-f0/max(f)*N2)
h_ffts = shift(h_fft,-f0/max(f)*N2)

; split spectrum
s1_fft = s_ffts*(f LE 0)
s2_fft = s_ffts*(f GT 0)
h1_fft = h_ffts*(f LE 0)
h2_fft = h_ffts*(f GT 0)

 plot, f*1e-6, abs(s1_fft), color='00ff00'x, xrange=[-20, 50], ytitle='amplitude spectrum s('+greek('omega')+')', xtitle='frequency (MHz)'
oplot, f*1e-6, abs(s2_fft), color='ff0000'x
oplot, f*1e-6, abs(s_fft), color='0050ff'x
 
 plot, f*1e-6, abs(h1_fft), color='00ff00'x, xrange=[-20, 50], ytitle='amplitude spectrum s('+greek('omega')+')', xtitle='frequency (MHz)'
oplot, f*1e-6, abs(h2_fft), color='ff0000'x
oplot, f*1e-6, abs(h_fft), color='0050ff'x


cs1_fft = s1_fft*h1_fft 
cs2_fft = s2_fft*h2_fft 

; shift baseband back to carrier
cs_mlook_ifft= shift(abs(fft(shift(cs1_fft,N2),1))+abs(fft(shift(cs2_fft,N2),1)),N2)/N


 plot, f*1e-6, abs(cs1_fft),   color='00ff00'x, Linestyle=0, xrange=[-20, 50], xtitle='Frequency (MHz)', ytitle='autocorrelation function g_1('+Greek('omega')+')' 
oplot, f*1e-6, abs(cs2_fft),   color='ff0000'x, Linestyle=0
;oplot, f*1e-6, (abs(cs1_fft)+abs(cs2_fft)),   color='0000ff'x, Linestyle=2, thick=2
;oplot, f*1e-6, abs(cs_mlook_fft),   color='0000ff'x, Linestyle=0, thick=2

; plot phases
phase_s1_fft   = phunwrap(atan(imag(s1_fft)  , real(s1_fft)))
phase_s2_fft   = phunwrap(atan(imag(s2_fft)  , real(s2_fft)))
phase_h1_fft   = phunwrap(atan(imag(h1_fft)  , real(h1_fft)))
phase_h2_fft   = phunwrap(atan(imag(h2_fft)  , real(h2_fft)))
phase_cs1_fft   = phunwrap(atan(imag(cs1_fft)  , real(cs1_fft)))
phase_cs2_fft   = phunwrap(atan(imag(cs2_fft)  , real(cs2_fft)))
 plot, f*1e-6, phase_s1_fft    - phase_s1_fft[N2],   color='00ff00'x, Linestyle=0, xrange=[-20, 50], xtitle='Frequency (MHz)', ytitle='phase of s('+Greek('omega')+')' 
oplot, f*1e-6, phase_s2_fft    - phase_s2_fft[N2],   color='ff0000'x, Linestyle=0
 plot, f*1e-6, phase_h1_fft    - phase_h1_fft[N2],   color='00ff00'x, Linestyle=0, xrange=[-20, 50], xtitle='Frequency (MHz)', ytitle='phase of s('+Greek('omega')+')' 
oplot, f*1e-6, phase_h2_fft    - phase_h2_fft[N2],   color='ff0000'x, Linestyle=0
 plot, f*1e-6, phase_cs1_fft    - phase_cs1_fft[N2-1],   color='00ff00'x, Linestyle=0, xrange=[-20, 50], yrange=[-0.05,0.05], xtitle='Frequency (MHz)', ytitle='phase of s('+Greek('omega')+')' 
oplot, f*1e-6, phase_cs2_fft    - phase_cs2_fft[N2],   color='ff0000'x, Linestyle=0


!P.MULTI = [7, 3, 3]
 ;plot, t*1e9, abs(cs_mlook_ifft)/max(abs(cs_mlook_ifft)) *max(abs(cs1)), color='00aaff'x,/NOERASE, linestyle=2, xrange=[-2e2, 2e2], yrange=[0, max(abs(cs1))]
 plot, t*1e9,    (cs_mlook_ifft)/max(  (cs_mlook_ifft))*max(abs(cs1)), color='ff00aa'x, linestyle=0,thick=2,/NOERASE, xrange=[-2e2, 2e2], yrange=[0, max(abs(cs1))]

cs_mlook_ifft_I  = abs(cs_mlook_ifft) GE 0.5*max(abs(cs_ifft));
cs_mlook_tau = (t(findfirst(cs_mlook_ifft_I,-1)) - t(findfirst(cs_mlook_ifft_I,1))) 
print, '     multi look compressed pulse length: ' + string(cs_mlook_tau*1e9) + ' ns (resolution = '+string(cs_mlook_tau*c) + 'm)'

stop
end

 
FUNCTION createfarray, from,step,to
  nsamples = ((to-from)/step)
  t = fltarr(nsamples)
  t = findgen(nsamples)
  t = t*step+from 
  return, t
end 

FUNCTION real, f
  return, real_part(f)
end

FUNCTION imag, f
  return, imaginary(f)
end

function findfirst, expression, which
  index = WHERE(expression)
  if which EQ 1 then begin
    return, index[0]
  endif else begin
    return, index[N_elements(index)-1]
  endelse  
    
end

function phunwrap, ph, tolerance=tol0, maxval=maxval0

  common phunwrap_common, idlver
  if n_elements(idlver) EQ 0 then begin
      idlver = !version.release
  endif


  if n_elements(maxval0) EQ 0 then maxval = 2d*!dpi else maxval = maxval0(0)
  if n_elements(tol0) EQ 0 then tol = 0.5*maxval else tol = tol0(0)*maxval

  if n_elements(ph) LT 2 then return, ph

  sz = size(ph)
  tp = sz(sz(0)+1)

  ;; First order difference 
  case tp of 
      12: dph = [0, long(ph)-long(ph(1:*))]
      13: dph = [0, long64(ph)-long64(ph(1:*))]
      15: dph = [0, long64(ph)-long64(ph(1:*))]
      else: dph = [0, ph - ph(1:*)]
  endcase
  
  p = maxval * (fix((dph GT tol) EQ 1) - fix((dph LT (-tol)) EQ 1))
  if idlver GT 5.25 then begin
      ;; Use built-in version if available
      r = total(p, /cumulative) 
  endif else begin
      ;; .. if not, then use the lame FOR loop
      r = p
      for i = 1L, n_elements(r)-1 do $
        r(i) = r(i) + r(i-1)
  endelse

  return, ph+r
end