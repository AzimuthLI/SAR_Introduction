c        = 2.998d8     ;speed of light (m/s)
Fs_range = 22.5e6      ;range sampling frequency (Hz)
bw_range = 20e6        ;range bandwidth (Hz)

; --------- Image 1 opening and reading ---------
openr, lun, 'data/i_xsar125141_bild1.dat',/xdr,/get_lun
dim1=lonarr(2)
readu, lun, dim1
data1 = complexarr(dim1)
readu, lun, data1
close, lun
free_lun, lun

; --------- Image 2 opening and reading ---------
openr, lun, 'data/i_xsar125141_bild2.dat',/xdr,/get_lun
dim2=lonarr(2)
readu, lun, dim2
data2 = complexarr(dim2)
readu, lun, data2
close, lun
free_lun, lun


ddim = dim1
