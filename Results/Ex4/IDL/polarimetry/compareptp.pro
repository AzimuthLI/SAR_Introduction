pro compareptp

; 5. Compute eigenvalues and eigenvectors.  Both analytically and using
;    built-in IDL function (LA_EIGENQL). 
;    Re-create using eqn 5.74 in Irena's thesis.
;        
;         Tip: READ IDL help of LA_EIGENQL and REVERSE functions carefully!
;        
;   

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

i=100
j=100
w1=complex(-0.5,sqrt(3)/2.)
w2=CONJ(w1)
l=complexarr(3)
a=-(t[0,0,i,j]+t[1,1,i,j]+t[2,2,i,j])
b=-((abs(t[0,1,i,j]))^2+(abs(t[0,2,i,j]))^2.+(abs(t[1,2,i,j]))^2-t[0,0,i,j]*t[1,1,i,j]-t[0,0,i,j]*t[2,2,i,j]-t[1,1,i,j]*t[2,2,i,j])
c=-(t[0,0,i,j]*t[1,1,i,j]*t[2,2,i,j]-t[0,0,i,j]*(abs(t[1,2,i,j]))^2-t[2,2,i,j]*(abs(t[0,1,i,j]))^2-t[1,1,i,j]*(abs(t[0,2,i,j]))^2+t[0,1,i,j]*conj(t[0,2,i,j])*t[1,2,i,j]+conj(t[0,1,i,j])*t[0,2,i,j]*conj(t[1,2,i,j]))
m=2*(a)^3-9*a*b+27*c
k=(a)^2-3*b
n=(m)^2-4*(k^3)
l[0]=-(a+((m+sqrt(n))/2.)^(1/3.)+((m-sqrt(n))/2.)^(1/3.))
l[1]=-(a+w2*(((m+sqrt(n))/2.)^(1/3.))+w1*(((m-sqrt(n))/2.)^(1/3.)))
l[2]=-(a+w1*(((m+sqrt(n))/2.)^(1/3.))+w2*(((m-sqrt(n))/2.)^(1/3.)))
l=float(l[sort(l)]/3.)
print,'analytical derivation of eigenvalue'
 print,l

u=complexarr(3,3)
for z=0,2 do begin
    qn=complex((conj(t[0,1,i,j])*(t[2,2,i,j]-l[z])-conj(t[0,2,i,j])*t[1,2,i,j])); N.B. complex(), qn , qd ,  q=qn/qd!!!
    qd=complex((conj(t[0,2,i,j])*(t[1,1,i,j]-l[z])-conj(t[0,1,i,j])*conj(t[1,2,i,j])))
    q=qn/qd 
    u[*,z]= [-(t[1,2,i,j]/conj(t[0,1,i,j]))-((t[1,1,i,j]-l[z])*q/conj(t[0,1,i,j])), q, 1]
    u_norm=sqrt(total(abs(u[*,z])^2))
    u[*,z]= u[*,z]/u_norm ;NORMALIZATION!!!!!
endfor
print,'analytical derivation of eigenvectors'
print,u

        
;Compute eigenvalues and eigenvectors using
;built-in IDL function (LA_EIGENQL). 
eigen_value=complexarr(3)
eigen_value=la_eigenql(t[*,*,i,j],eigenvectors=eigen_vectors)
print,'la_eigenql eigenvalue'
print,eigen_value
print,'la_eigenql eigenvectors'
print,eigen_vectors

;    Compare (statisics) T-matrix re-created using eigenvectors and
;    eigenvalues (both solutions).

lambda_a=fltarr(3,3)
for k=0,2 do begin
         	lambda_a[k,k]=l[k]
       end
recreated_a=(u)#lambda_a#transpose(conj(u))
print,'original coherence matrix'
print,t[*,*,i,j]
print,"coherence matrix re-created using eigenvectors and eigenvalue (eigenvalues and eigenvectors computed analytically)"
print,recreated_a
print,"Max relative error (%)= ", max((abs(recreated_a)-abs(t(*,*,i,j)))/abs(t(*,*,i,j)))*100.


lambda_b=fltarr(3,3)
for k=0,2 do begin
         	lambda_b[k,k]=eigen_value[k]
       end
recreated_b=conj(eigen_vectors)#lambda_b#transpose((eigen_vectors))
print,"coherence matrix re-created using eigenvectors and eigenvalue(eigenvalues and eigenvectors using built-in IDL function (LA_EIGENQL))"
print,recreated_b
print,"Max relative error (%)= ", max((abs(recreated_b)-abs(t(*,*,i,j)))/abs(t(*,*,i,j)))*100.


;    Re-create using eqn 5.74 in Irena's thesis.

ta1=complexarr(3,3)
ta2=complexarr(3,3)
ta3=complexarr(3,3)
t_a=complexarr(3,3)
 ta1=l[0]*u[*,0]#transpose(conj(u[*,0])) ;
 ta2=l[1]*u[*,1]#transpose(conj(u[*,1]))
 ta3=l[2]*u[*,2]#transpose(conj(u[*,2]))
t_a=ta1+ta2+ta3
print,"coherence matrix re-created using  5.74 in Irena's thesis (eigenvalues and eigenvectors computed analytically)"
print,t_a
print,"Max relative error (%)= ", max((abs(t_a)-abs(t(*,*,i,j)))/abs(t(*,*,i,j)))*100.

t1=complexarr(3,3)
t2=complexarr(3,3)
t3=complexarr(3,3)
t_b=complexarr(3,3)
 t1=eigen_value[0]*conj(eigen_vectors[*,0])#transpose((eigen_vectors[*,0]))
 t2=eigen_value[1]*conj(eigen_vectors[*,1])#transpose((eigen_vectors[*,1]))
 t3=eigen_value[2]*conj(eigen_vectors[*,2])#transpose((eigen_vectors[*,2]))
t_b=t1+t2+t3
print,"coherence matrix re-created using 5.74 in Irena's thesis (eigenvalues and eigenvectors using built-in IDL function (LA_EIGENQL))"
print,t_b
print,"Max relative error (%)= ", max((abs(t_b)-abs(t(*,*,i,j)))/abs(t(*,*,i,j)))*100.

stop
;OUTPUT
; 
; analytical derivation of eigenvalue
;       5559.28      57438.7      174501.
; analytical derivation of eigenvectors
; (   -0.0791045,   -0.0288829)(    0.0291266,   -0.0537290)(     0.994572,      0.00000)
; (     0.226690,     0.314260)(    -0.150332,     0.905965)(    0.0804990,      0.00000)
; (     0.916547,    0.0519919)(    -0.255856,    -0.295686)(    0.0659261,      0.00000)
; % Loaded DLM: LAPACK.
; la_eigenql eigenvalue
;       5559.31      57438.7      174501.
; la_eigenql eigenvectors
; (    0.0842099,      0.00000)(  -0.00893325,   -0.0604586)(    -0.934245,    -0.341115)
; (     0.387489,      0.00000)(     0.646806,    -0.651931)(    0.0470937,    0.0652861)
; (    -0.918020,     -0.00000)(     0.272192,    -0.280721)(   -0.0658203,  -0.00373372)


; original coherence matrix
; (      155726.,      0.00000)(     -29212.2,     -30431.9)(      11154.9,     -1891.49)
; (     -29212.2,      30431.9)(      75142.8,      0.00000)(     -3477.46,     -490.278)
; (      11154.9,      1891.49)(     -3477.46,      490.280)(      6629.76, -7.34913e-05)


; coherence matrix re-created using eigenvectors and eigenvalue (eigenvalues and eigenvectors computed analytically)
; (      155726.,  0.000244141)(     -29212.2,     -30431.9)(      11154.9,     -1891.49)
; (     -29212.2,      30431.9)(      75142.8, -0.000488281)(     -3477.47,     -490.274)
; (      11154.9,      1891.49)(     -3477.47,      490.274)(      6629.72,      0.00000)
; Max relative error (%)=   9.73266e-05


; coherence matrix re-created using eigenvectors and eigenvalue(eigenvalues and eigenvectors using built-in IDL function (LA_EIGENQL))
; (      155726.,      0.00000)(     -29212.2,     -30431.9)(      11154.9,     -1891.49)
; (     -29212.2,      30431.9)(      75142.8, -0.000976562)(     -3477.46,     -490.279)
; (      11154.9,      1891.49)(     -3477.46,      490.279)(      6629.75,  0.000122070)
; Max relative error (%)=   1.72627e-05
; coherence matrix re-created using  5.74 in Irena's thesis (eigenvalues and eigenvectors computed analytically)
; (      155726.,  0.000244141)(     -29212.2,     -30431.9)(      11154.9,     -1891.49)
; (     -29212.2,      30431.9)(      75142.8, -0.000488281)(     -3477.47,     -490.274)
; (      11154.9,      1891.49)(     -3477.47,      490.274)(      6629.72,      0.00000)
; Max relative error (%)=   9.73266e-05
; coherence matrix re-created using 5.74 in Irena's thesis (eigenvalues and eigenvectors using built-in IDL function (LA_EIGENQL))
; (      155726.,      0.00000)(     -29212.2,     -30431.9)(      11154.9,     -1891.49)
; (     -29212.2,      30431.9)(      75142.8, -0.000976562)(     -3477.46,     -490.279)
; (      11154.9,      1891.49)(     -3477.46,      490.279)(      6629.75,  0.000122070)
; Max relative error (%)=   1.72627e-05






end