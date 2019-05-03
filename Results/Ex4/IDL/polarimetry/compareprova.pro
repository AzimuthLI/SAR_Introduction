pro compareprova
 
; 5. Compute eigenvalues and eigenvectors.  Both analytically and using
;    built-in IDL function (LA_EIGENQL). 
;    Compare eigenvalues from built-in and analytic solutions (compute
;    relative % difference for ~1000 pixels and other statistics)


device,retain=2,decomposed=0

OPENR, uh, '/home/pisc_io/Documents/IDL/polarimetry/i_al_af_1206_hh_corr.dat', /xdr, /get_lun
dim = lonarr(2)
READU, uh, dim
hh = complexarr(dim)
READU, uh, hh
OPENR, uv, '/home/pisc_io/Documents/IDL/polarimetry/i_al_af_1206_vv_corr.dat', /xdr, /get_lun
vv = complexarr(dim)
READU, uv, vv
OPENR, ux,'/home/pisc_io/Documents/IDL/polarimetry/i_al_af_1206_xx_corr.dat', /xdr, /get_lun
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


; 5. Compute eigenvalues and eigenvectors analytically 
dim0=20
n1=(dim(0)/2)-dim0/2
n2=n1+dim0-1
dim1=50
m1=(dim(1)/2)-dim1/2
m2=m1+dim1-1
l=complexarr(3,dim(0),dim(1))
w1=complex(-0.5,sqrt(3.)/2.)
w2=CONJ(w1)
for i=n1,n2 do begin;10 do begin;
    for j=m1,m2 do begin;10 do begin;
print,i,j,' l'
a=-(t[0,0,i,j]+t[1,1,i,j]+t[2,2,i,j])
b=-((abs(t[0,1,i,j]))^2+(abs(t[0,2,i,j]))^2+(abs(t[1,2,i,j]))^2-t[0,0,i,j]*t[1,1,i,j]-t[0,0,i,j]*t[2,2,i,j]-t[1,1,i,j]*t[2,2,i,j])
c=-(t[0,0,i,j]*t[1,1,i,j]*t[2,2,i,j]-t[0,0,i,j]*(abs(t[1,2,i,j]))^2-t[2,2,i,j]*(abs(t[0,1,i,j]))^2-t[1,1,i,j]*(abs(t[0,2,i,j]))^2+t[0,1,i,j]*conj(t[0,2,i,j])*t[1,2,i,j]+conj(t[0,1,i,j])*t[0,2,i,j]*conj(t[1,2,i,j]))
m=2*(a)^3-9*a*b+27*c
k=(a)^2-3*b
n=(m)^2-4*(k^3)
l[0,i,j]=-(a+((m+sqrt(n))/2.)^(1/3.)+((m-sqrt(n))/2.)^(1/3.))
l[1,i,j]=-(a+w2*(((m+sqrt(n))/2.)^(1/3.))+w1*(((m-sqrt(n))/2.)^(1/3.)))
l[2,i,j]=-(a+w1*(((m+sqrt(n))/2.)^(1/3.))+w2*(((m-sqrt(n))/2.)^(1/3.)))
l1=complexarr(3,dim(0),dim(1))
l1=l[*,i,j]
l1=l1[sort(l1)]
l[*,i,j]=l1
endfor
endfor
l=float(l/3.)


u=complexarr(3,3,dim(0),dim(1))
for i=0,dim(0)-1 do begin;10 do begin;
    for j=0,dim(1)-1 do begin;10 do begin;
print,i,j,' u'
for z=0,2 do begin
    qn=complex((conj(t[0,1,i,j])*(t[2,2,i,j]-l[z,i,j])-conj(t[0,2,i,j])*t[1,2,i,j])); N.B. complex(), qn , qd ,  q=qn/qd!!!
    qd=complex((conj(t[0,2,i,j])*(t[1,1,i,j]-l[z,i,j])-conj(t[0,1,i,j])*conj(t[1,2,i,j])))
    q=qn/qd 
    u[*,z,i,j]= [-(t[1,2,i,j]/conj(t[0,1,i,j]))-((t[1,1,i,j]-l[z])*q/conj(t[0,1,i,j])), q, 1]
    u_norm=sqrt(total(abs(u[*,z,i,j])^2))
    u[*,z,i,j]= u[*,z,i,j]/u_norm ;NORMALIZATION!!!!!
endfor
endfor
endfor
;


;Compute eigenvalues and eigenvectors using
;built-in IDL function (LA_EIGENQL). 
eigen_value=complexarr(3,dim(0),dim(1))
for i=0,dim(0)-1 do begin;10 do begin;
    for j=0,dim(1)-1 do begin;10 do begin;
       eigen_value[*,i,j]=la_eigenql(t[*,*,i,j],eigenvectors=eigen_vectors)
endfor
endfor


;Compare eigenvalues from built-in and analytic solutions (compute
;    relative % difference for ~1000 pixels and other statistics)

;diff=complexarr(3,dim(0),dim(1))
diff=complexarr(3,dim0,dim1)

for i=0,dim0-1 do begin;10 do begin;
    for j=0,dim1-1 do begin;10 do begin;
       diff[*,i,j]= (abs(l[*,i+n1,j+m1]-eigen_value[*,i+n1,j+m1]))*100./l[*,i+n1,j+m1]
endfor
endfor
print,'relative difference'
print,mean(abs(diff))


stop

end