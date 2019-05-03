pro ex_eigenvalue

n=3

t=[[-3,1,5],[1,0,-2],[5,-2,4]]
print,'matrix'
print,t
eigen_value=la_eigenql(t,eigenvectors=eigen_vectors)
print,'la_eigenql eigenvalue'
print,eigen_value
print,'la_eigenql eigenvectors'
print,eigen_vectors

a=-(t[0,0]+t[1,1]+t[2,2])
b=-((abs(t[0,1]))^2.+(abs(t[0,2]))^2+(abs(t[1,2]))^2.-t[0,0]*t[1,1]-t[0,0]*t[2,2]-t[1,1]*t[2,2])
c=-(t[0,0]*t[1,1]*t[2,2]-t[0,0]*(abs(t[1,2]))^2.-t[2,2]*(abs(t[0,1]))^2.-t[1,1]*(abs(t[0,2]))^2.+t[0,1]*conj(t[0,2])*t[1,2]+conj(t[0,1])*t[0,2]*conj(t[1,2]))

m=2*(a)^3.-9*a*b+27*c
k=(a)^2.-3.*b
n=(m)^2.-4.*(k^3.)

w1=complex(-0.5,sqrt(3)/2.)
w2=CONJ(w1)
l=complexarr(3)
l[0]=-(a+((m+sqrt(n))/2.)^(1/3.)+((m-sqrt(n))/2)^(1/3.))
l[1]=-(a+w2*(((m+sqrt(n))/2.)^(1/3.))+w1*(((m-sqrt(n))/2.)^(1/3.)))
l[2]=-(a+w1*(((m+sqrt(n))/2.)^(1/3.))+w2*(((m-sqrt(n))/2.)^(1/3.)))
l=float(l/3.)
print,'analytical derivation of eigenvalue'
print,l

u=fltarr(3,3)
for k=0,2 do begin
    q1=conj(t[0,1])*(t[2,2]-l[k])-conj(t[0,2])*t[1,2]
    q2=conj(t[0,2])*(t[1,1]-l[k])-conj(t[0,1])*conj(t[1,2])
    q=q1/q2
    u[k,*]= transpose([-(t[1,2]/conj(t[0,1]))-((t[1,1]-l[k])*q/conj(t[0,1])), q, 1])
endfor
print,'analytical derivation of eigenvectors'
print,u

ta1=fltarr(3,3)
ta2=fltarr(3,3)
ta3=fltarr(3,3)
tta=fltarr(3,3)
ta1=l[0]*(u[*,0]#transpose(conj(u[*,0])))
ta2=l[1]*(u[*,1]#transpose(conj(u[*,1])))
ta3=l[2]*(u[*,2]#transpose(conj(u[*,2])))
; ta1=l[0]*(u[0,*]##transpose(conj(u[0,*])))
; ta2=l[1]*(u[1,*]##transpose(conj(u[1,*])))
; ta3=l[2]*(u[2,*]##transpose(conj(u[2,*])))
ta=ta1+ta2+ta3
print,"coherence matrix re-created using  5.74 in Irena's thesis (eigenvalues and eigenvectors computed analytically)"
print,ta

tb1=complexarr(3,3)
tb2=complexarr(3,3)
tb3=complexarr(3,3)
tb=complexarr(3,3)
tb1=eigen_value[0]*(eigen_vectors[*,0]#transpose(conj(eigen_vectors[*,0])))
tb2=eigen_value[1]*(eigen_vectors[*,1]#transpose(conj(eigen_vectors[*,1])))
tb3=eigen_value[2]*(eigen_vectors[*,2]#transpose(conj(eigen_vectors[*,2])))
; tb1=eigen_value[0]*(eigen_vectors[0,*]##transpose(conj(eigen_vectors[0,*])))
; tb2=eigen_value[1]*(eigen_vectors[1,*]##transpose(conj(eigen_vectors[1,*])))
; tb3=eigen_value[2]*(eigen_vectors[2,*]##transpose(conj(eigen_vectors[2,*])))
tb=tb1+tb2+tb3
print,"coherence matrix re-created using 5.74 in Irena's thesis (eigenvalues and eigenvectors using built-in IDL function (LA_EIGENQL))"
print,tb


stop
; maxErr=0d
; for i=0,n-1 do begin
; ;m*eigenvectors=eigenvalue*eigenvectors
;     alhs=mat ## eigen_vectors[*,i]
;     arhs=eigen_value[i]*eigen_vectors[*,i]
;     maxErr=maxErr>(Abs(alhs-arhs))
; end
; print,'la_eigenql max error'
; print,maxErr
; 
; ; Reduce to tridiagonal form  
; q = mat    ; make a copy  
; LA_TRIRED, q, d, e  
;   
; ; Compute eigenvalues and eigenvectors  
; eigenvalues = d  
; eigenvectors = q  
; LA_TRIQL, eigenvalues, e, eigenvectors  
; PRINT, 'LA_TRIQL eigenvalues:'  
; PRINT, eigenvalues  
; PRINT, 'LA_TRIQL eigenvectors:'  
; PRINT, eigenvectors  
; 
; 
; coef=[-41,55,-15,1]
; lambda=complexarr(3)
; lambda=imsl_zeropoly(coef)



end