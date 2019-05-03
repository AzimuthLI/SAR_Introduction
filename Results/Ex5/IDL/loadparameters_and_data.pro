
hh1 = loadxdr('data/rvog_t1_hh.dat')
vv1 = loadxdr('data/rvog_t1_vv.dat')
xx1 = loadxdr('data/rvog_t1_xx.dat')
hh2 = loadxdr('data/rvog_t2_hh.dat')
vv2 = loadxdr('data/rvog_t2_vv.dat')
xx2 = loadxdr('data/rvog_t2_xx.dat')
hh3 = loadxdr('data/rvog_t3_hh.dat')
vv3 = loadxdr('data/rvog_t3_vv.dat')
xx3 = loadxdr('data/rvog_t3_xx.dat')

ddim = size(hh1,/dim)
;constants for sim_data_rvog

H        = 2.880e3        ;(sensor height, m)
lambda   = 0.24          ;(wavelength, L-band, m)
B_12     = -10           ;(horz baseline btw Im1&2, m)
B_13     = -20           ;(horz baseline btw Im1&3, m)
W        = 100e6         ;(bandwidth in range, Hz)
grng_res = 0.5           ;(ground rng pixel spacing, m)
theta    = 45*!DTOR      ;(angle of incidence to img centre, assume constant for small area) (rads)
c        = 3e8           ;(speed of light, m/s)
R0       = H/cos(theta)  ;broadside range (m)
alpha    = 0             ;local slope 
k        = 4*!PI/lambda  ;wave number

pscl = sin(theta)

; calculate size and position of plotting window
wwidth  = 0.36
scrsize = get_screen_size()
wscr = scrsize[0]
hscr = scrsize[1]
wdim = [wwidth*wscr,wwidth*wscr*sqrt(2)]
wloc = [wscr - wwidth,0]


