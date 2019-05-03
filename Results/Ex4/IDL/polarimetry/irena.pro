pro irena

; 5. Compute eigenvalues and eigenvectors.  Both analytically and using
;    built-in IDL function (LA_EIGENQL). 
;    Compare eigenvalues from built-in and analytic solutions (compute
;    relative % difference for ~1000 pixels and other statistics)
;    Compare (statisics) T-matrix re-created using eigenvectors and
;    eigenvalues (both solutions).
;    Re-create using eqn 5.74 in Irena's thesis.
;        
;         Tip: READ IDL help of LA_EIGENQL and REVERSE functions carefully!
;        
;   


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
im=complexarr(3,dim(0),dim(1))
im[0,*,*]=(hh+vv)
im[1,*,*]=(hh-vv)
im[2,*,*]=(2.*xx)
im=im/(sqrt(2))
;coherency matrix.
t=complexarr(3,3,dim(0),dim(1))
v=complexarr(3)
for i=0,dim(0)-1 do begin
    for j=0,dim(1)-1 do begin
    	v[0]=im[0,i,j]
	v[1]=im[1,i,j]
	v[2]=im[2,i,j]
	t(*,*,i,j)=transpose(v)##(conj(v))
    end
end
t=smooth(t,[1,1,7,7])
;provare per singoli valori di i e j;
;poi generalizzare con un ciclo for tutti e due i metodi
;quindi fare la differenza relativa su 1000 pixel (dopo aver controllato ke l'ordine degli autovalori sia lo stesso con if,<,>)

i=200
j=300

a=-(t[0,0,i,j]+t[1,1,i,j]+t[2,2,i,j])
b=-((abs(t[0,1,i,j]))^2.+(abs(t[0,2,i,j]))^2.+(abs(t[1,2,i,j]))^2.-t[0,0,i,j]*t[1,1,i,j]-t[0,0,i,j]*t[2,2,i,j]-t[1,1,i,j]*t[2,2,i,j])
c=-(t[0,0,i,j]*t[1,1,i,j]*t[2,2,i,j]-t[0,0,i,j]*(abs(t[1,2,i,j]))^2.-t[2,2,i,j]*(abs(t[0,1,i,j]))^2.-t[1,1,i,j]*(abs(t[0,2,i,j]))^2.+t[0,1,i,j]*conj(t[0,2,i,j])*t[1,2,i,j]+conj(t[0,1,i,j])*t[0,2,i,j]*conj(t[1,2,i,j]))
m=2.*(a)^3.-9.*a*b+27.*c
k=(a)^2.-3.*b
n=(m)^2.-4.*(k^3.)
w1=complex(-0.5,sqrt(3.)/2.)
w2=conj(w1)

l=fltarr(3)
l[0]=-(a+((m+sqrt(n))/2.)^(1/3.)+((m-sqrt(n))/2.)^(1/3.))
l[1]=-(a+w2*(((m+sqrt(n))/2.)^(1/3.))+w1*(((m-sqrt(n))/2.)^(1/3.)))
l[2]=-(a+w1*(((m+sqrt(n))/2.)^(1/3.))+w2*(((m-sqrt(n))/2.)^(1/3.)))

l=l[sort(l)]/3.
print,"eigenvalues analitically computed",l

eigen_value=fltarr(3)
eigen_value=la_eigenql(t[*,*,i,j],eigenvectors=eigen_vectors)
print,'la_eigenql eigenvalue'
print,eigen_value
print,'la_eigenql eigenvectors'
print,eigen_vectors
stop
CLOSE, uh
 free_lun, uh
CLOSE, uv
 free_lun, uv
CLOSE, ux
 free_lun, ux


end
