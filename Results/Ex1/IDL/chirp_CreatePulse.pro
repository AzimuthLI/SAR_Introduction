; this is a batch-file, called by "@CreateChirpedPulse"

; define time axis
t0   = -1e-6     ; (sec) from
tend =  1e-6     ; (sec) to
tres =  3e-10    ; (sec) resolution

; define parameters chirped pulse
f0  = 20e6       ; (Hz)  Center frequency of pulse
bw  = 20e6       ; (Hz)  Bandwidth of linear chirped pulse
tau =  1.5e-6    ; (sec) Puls length, with puls center at t = 0

bc  = bw/tau; frequency change rate
w0  = f0*2*!pi ; (1/s) center frequency

i = complex(0,1) ; complex unity
c = 2.998e8;

fmin = 0         ; lower limit of frequency axis of plots
fmax = 45e6      ; upper limit of frequency axis of plots
fmin_phi = -20e6    ; lower limit of frequency axis of plots

; create time axis
t = createfarray(t0,tres,tend)
N = N_elements(t)
N2 = round(N/2)

; create frequency axis
f = createfarray(-1/(2*tres),1/(tres*N),1/(2*tres))
; find Index of center frequency
f0I = findfirst(f GE (f0))
f00I = findfirst(f GE (0))
frange = where(f ge fmin and f le fmax)
frange_phi = where(f ge fmin_phi and f le fmax)

; create chirped pulse
s = complexarr(N)
sI    = WHERE((t GE -tau/2) AND (t LE tau/2))
s[sI] =   exp(i*t(sI)*(w0 + !pi*bc*t(sI)))
