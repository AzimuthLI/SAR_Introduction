pro pol2

;COMPUTE AND VISUALIZE PAULI IMAGE (CAN SMOOTH TO IMPROVE THE RESULTS)
;STORE BYTE-SCALED DATA IN 3-D ARRAY AND
;VISUALIZE WITH TV,ARRAY,TRUE=3

device,retain=2,decomposed=0

OPENR, uh, '/home/pisc_io/Documents/i_al_af_1206_hh_corr.dat', /xdr, /get_lun
dim = lonarr(2)
READU, uh, dim
hh = complexarr(dim)
READU, uh, hh
OPENR, uv, '/home/pisc_io/Documents/i_al_af_1206_vv_corr.dat', /xdr, /get_lun
vv = complexarr(dim)
READU, uv, vv
OPENR, ux, '/home/pisc_io/Documents/i_al_af_1206_xx_corr.dat', /xdr, /get_lun
xx = complexarr(dim)
READU, ux, xx

;PAULI IMAGE
r = complexarr(dim)
r=(hh-vv)
g = complexarr(dim)
g=(2.*xx)
b = complexarr(dim)
b=(hh+vv)
window, 1,title = "R channel"
tv, bytscl(abs(r),0,2.*mean(abs(r)))
window, 2, title = "G channel"
tv, bytscl(abs(g),0,2.*mean(abs(g)))
window, 3, title = "B channel"
tv, bytscl(abs(b),0,2.*mean(abs(b)))

;STORE BYTE-SCALED DATA IN 3-D ARRAY
im=complexarr(dim(0),dim(1),3)
im[*,*,0]=abs(r)*255/(2.*mean(abs(r))) ;normalization to improve the visualization
im[*,*,1]=abs(g)*255/(2.*mean(abs(g)))
im[*,*,2]=abs(b)*255/(2.*mean(abs(b)))
window,4, title = "Pauli image"
tv,abs(im),true=3

;SMOOTH TO IMPROVE THE RESULTS
im_s7=smooth(im,[7,7,1])
window,5, title = "Pauli image smoothed (winsz=7x7)"
tv,abs(im_s7),true=3
im_s3=smooth(im,[3,3,1])
window,6, title = "Pauli image smoothed(winsz=3x3)"
tv,abs(im_s3),true=3

stop

CLOSE, uh
 free_lun, uh
CLOSE, uv
 free_lun, uv
CLOSE, ux
 free_lun, ux

end
