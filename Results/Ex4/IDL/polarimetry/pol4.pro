pro pol4
;Compute coherency matrix.
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

;PAULI SCATTERING VECTOR
im=complexarr(dim(0),dim(1),3)
im[*,*,0]=hh+vv
im[*,*,1]=hh-vv
im[*,*,2]=2.*xx

;coherency matrix.
coher=complexarr(dim(0),dim(1),3,3)
v=complexarr(3)
for i=0,dim(0)-1 do begin
    for j=0,dim(1)-1 do begin
    	v[0]=im[i,j,0]
	v[1]=im[i,j,1]
	v[2]=im[i,j,2]
	coher(i,j,*,*)=(v)#transpose(conj(v))
    end
end

coher_s=smooth(coher,[7,7,1,1])
coher_s=smooth(coher,[7,7,1,1])

stop

CLOSE, uh
 free_lun, uh
CLOSE, uv
 free_lun, uv
CLOSE, ux
 free_lun, ux

end