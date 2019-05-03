pro eigenvalue

n=3
id=identity(n)
t=3*id
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
l=l/3.
print,'analytical derivation of eigenvalue'
print,l
u=complexarr(3,3)
for i=0,2 do begin
    q=(conj(t[0,1])*(t[2,2]-l[i])-conj(t[0,2])*t[1,2])/(conj(t[0,2])*(t[1,1]-l[i])-conj(t[0,1])*conj(t[1,2]))help,q

    u[*,i]= [-(t[1,2]/conj(t[0,1]))-((t[1,1]-l[i])/conj(t[0,1]))*q, q, 1]
endfor
print,'analytical derivation of eigenvectors'
print,u






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
; ;s=(n-lambda*M)
; coef=[-9,27,-9,1]
; lambda=complexarr(3)
; ;lambda=imsl_zeropoly(coef)



stop
end