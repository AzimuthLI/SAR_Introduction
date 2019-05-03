; define (natural) constants
c = 2.998d8          ;speed of light (m/s)
j = complex(0.0,1.0) ; complex unity

; get screen dimensions
scrsize = get_screen_size()

; load parameters for data
case subex of
1: begin
    R0        = 1000d    ;(m)   closest range
    v         = 70d       ;(m/s) sensor velocity
    lambda    = 0.0566d  ;(m)   wavelength
    l_antenna = 2d        ;(m)   antenna length
    bw        = 50d6     ;(Hz)  band width
    tau       = 5d-6     ;(s)   pulse width
    fs_range  = 100d6    ;(Hz)  sampling frequency in range
    prf       = 400d      ;(Hz)  pulse repetition frequency
    datfile   = 'data/chirp_2d_test.dat'
    outsuffix   = 'new_h'
    posx = [0,1,2]/3.0+0.04
    imgw = 0.29
    posy = [0,1]/2.0+0.06
    imgh = 0.38
    winlocation1 = [scrsize(0)/2,scrsize(1)-50]
    windim1      = [scrsize(0)/2,scrsize(1)-100]
    asp = 1
    colortable = 13
    focus_range = 1
    nlscalemodel = 1        
   end
2: begin
    rg_delay  = 0.00551973   ;(s)   range to first pulse
    ts_rg     = 5.274261d-8  ;(s)   sampling period (in range)
    Res_az    = 6d            ;(m)   azimuth resolution

    R0        = c*rg_delay/2 ;(m)   closest range 
    v         = 7055.3079    ;(m/s) sensor velocity
    lambda    = 0.056666     ;(m)   wavelength
    l_antenna = 2*res_az     ;(m)   antenna length
    bw        = -1.55404d7    ;(Hz)  bandwidth
    tau       = 3.71d-5      ;(s)   pulse width
    fs_range  = 1/ts_rg      ;(Hz)  sampling frequency in range  
    prf       = 1678.712     ;(Hz)  pulse repetition frequency
    datfile   = 'data/ers_raw_demo.dat'    
    outsuffix   = 'new_h'
    posx = [0,1,2]/3.0+0.09
    imgw = 0.22
    posy = [0,0.75]/2.0+0.13
    imgh = 0.36
    winlocation1 = [scrsize(0)/2,scrsize(1)-50]
    windim1      = [scrsize(0)/2,scrsize(1)*2/3]        
    asp = 0.25
    colortable = 0
    focus_range = 1
    nlscalemodel = 3        
   end

3: begin
    rg_delay  = 9.4d-6       ;(s)   range delay R0=c*t/2. (i.e. range to first pulse, s)
    res_az    = 1.5          ;(m)   resolution in azimuth

    ; the following data is not required, as the data is already range compressed
    ;bw        =             ;(Hz)  bandwidth
    ;tau       =              ;(s)   pulse width
    ;fs_range  =              ;(Hz)  sampling frequency in range  
    R0        = c*rg_delay/2 ;(m)   closest range 
    v         = 72.6         ;(m/s) sensor velocity
    lambda    = 0.056666     ;(m)   wavelength
    l_antenna = 2*res_az     ;(m)   antenna length
    prf       = 952.38       ;(Hz)  pulse repetition frequency

    fs_range  = 200d6        ; actually: not known!! But set to an arbitrary value...    
    datfile   = 'data/rdemo040689_cmp.dat'
    outsuffix   = 'new_h'
    focus_range = 0
    
    posx = [0,0,1]/3.0+0.09
    imgw = 0.25
    posy = [0,0.87]/2.0+0.11
    imgh = 0.3
    winlocation1 = [scrsize(0)/2,scrsize(1)-50]
    windim1      = [scrsize(0)/2,scrsize(1)*2/3]        
    asp = 1.0/16
    colortable = 0    
    nlscalemodel = 5        
   end
endcase


; open and load image
; -----------------------
openr, unit, datfile,/xdr,/get_lun

; get dimensions of binary data (1st 8 bytes of file)
dim = lonarr(2)
readu, unit, dim

; create complex array which receives the data and load data
data = complexarr(dim)
readu, unit,data

; close file and free file unit
close, unit
free_lun, unit

; set parameters which contain the image size (px)
N_rg = dim(0)
N_az = dim(1)

; generate range and azimut axis
;if (focus_range) then begin 
  t_range   = (findgen(dim(0))-dim(0)/2.0)/fs_range
;end
t_range2  = (findgen(dim(0)/2)-dim(0)/4.0)/fs_range*2
y_azimut  = (findgen(dim(1))-dim(1)/2.0)*v/prf
y_azimut2 = (findgen(dim(1)/2)-dim(1)/4.0)*v/prf*2
  