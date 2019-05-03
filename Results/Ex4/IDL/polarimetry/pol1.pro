pro pol1
;POLSAR INTRODUCTION COURSE
;VISUALIZE HH CHANNEL,VV CHANNEL,XX CHANNEL

device,retain=2,decomposed=0

OPENR, uh, '/home/pisc_io/Documents/i_al_af_1206_hh_corr.dat', /xdr, /get_lun
dim = lonarr(2)
READU, uh, dim
hh = complexarr(dim)
READU, uh, hh
Rg=dim(1)
Az=dim(0)

OPENR, uv, '/home/pisc_io/Documents/i_al_af_1206_vv_corr.dat', /xdr, /get_lun
dim = lonarr(2)
READU, uv, dim
vv = complexarr(dim)
READU, uv, vv

OPENR, ux, '/home/pisc_io/Documents/i_al_af_1206_xx_corr.dat', /xdr, /get_lun
dim = lonarr(2)
READU, ux, dim
xx = complexarr(dim)
READU, ux, xx

window, 1, xsize=Rg/3,ysize=Az/3, title = "HH channel"
tv, bytscl(congrid(abs(hh),Rg/3,Az/3), min(abs(hh)),2.5*mean(abs(hh)))
window, 2, xsize=Rg/3,ysize=Az/3, title = "VV channel"
tv, bytscl(congrid(abs(vv),Rg/3,Az/3), min(abs(vv)),2.5*mean(abs(vv)))
window, 3, xsize=Rg/3,ysize=Az/3, title = "XX channel"
tv, bytscl(congrid(abs(xx),Rg/3,Az/3), min(abs(xx)),2.5*mean(abs(xx)))

stop

CLOSE, uh
 free_lun, uh
CLOSE, uv
 free_lun, uv
CLOSE, ux
 free_lun, ux


end