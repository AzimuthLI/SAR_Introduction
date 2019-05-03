function eigenvaluesT3, T3, eigenvectors=eigvect
;  analytical solution, taken from appendix A6 from the book:
; "Polarimetric Radar imaging" from Jong-Sen Lee and Eric Pottier

a0 = reform(T3[*,*,0,0]) 
b0 = reform(T3[*,*,1,1]) 
c0 = reform(T3[*,*,2,2])
z1 = reform(T3[*,*,1,0])
z2 = reform(T3[*,*,2,0])
z3 = reform(T3[*,*,2,1])

detT3 = a0*b0*c0 - c0*z1*conj(z1) - b0*z2*conj(z2) + z1*conj(z2)*z3 + conj(z1)*z2*conj(z3) - a0*z3*conj(z3)
trT3  = a0 + b0 + c0

A = a0*b0 + a0*c0 + b0*c0 - z1*conj(z1) - z2*conj(z2)  - z3*conj(z3)
B = a0^2 - a0*b0 + b0^2 - a0*c0 - b0*c0 + c0^2 + 3*(z1*conj(z1) + z2*conj(z2) + z3*conj(z3))
C = 27*detT3 - 9*A*trT3 + 2*trT3^3 + sqrt((27*detT3 - 9*A*trT3 + 2*trT3^3)^2-4*B^3)


; a factor of 0.5 was removed here in contrast to the book.
lambda1 = (trT3/3 + 2.0^(1.0/3)*B/(3*C^(1.0/3)) + C^(1.0/3)/(3*2.0^(1.0/3)))
lambda2 = (trT3/3 - complex(1, sqrt(3))*B/(3*2.0^(2.0/3)*C^(1.0/3)) - complex(1,-sqrt(3))*C^(1.0/3)/(6*2.0^(1.0/3)))
lambda3 = (trT3/3 - complex(1,-sqrt(3))*B/(3*2.0^(2.0/3)*C^(1.0/3)) - complex(1, sqrt(3))*C^(1.0/3)/(6*2.0^(1.0/3)))


eigenvals = real([[[lambda3]],[[lambda2]],[[lambda1]]])
ddim = size(T3)
ddim = ddim [1:2]

ev_order = intarr([ddim,3])
; sort eigenvalues in descending order
for i=0,ddim(0)-1 do for j=0,ddim(1)-1 do begin
  ev_order[i,j,*] = reverse(sort(eigenvals[i,j,*]))
  eigenvals[i,j,*] = eigenvals[i,j,ev_order[i,j,*]]
end

;lambda1 = eigenvals[*,*,0]
;lambda2 = eigenvals[*,*,1]
;lambda3 = eigenvals[*,*,2]


;eigenvectors
; check if row-columns are ordered correctly!!
if ARG_PRESENT(eigvect) then begin
  (Scope_varfetch(scope_varname(eigvect,level=-1),/ENTER,level=-1)) = complexarr([ddim,3,3])
  print, 'returning eigenvectors as rows (here is an error in the code or in the book page 383 of Jong-Sen-Lee)'
  ; test with the matrix: T[0,0,*,*] =[[-1, 1, 3],[1,2, 0],[3, 0, 2]]
  eigvect[*,*,*,0] = [(lambda1-c0)/conj(z2)+((lambda1-c0)*conj(z1)+conj(z2)*z3)*conj(z3)/(((b0-lambda1)*conj(z2)-conj(z1)*conj(z3))*conj(z2)),((lambda1-c0)*conj(z1)+conj(z2)*z3)/((b0-lambda1)*conj(z2)-conj(z1)*conj(z3)),1+complexarr(ddim)]
  eigvect[*,*,*,1] = [(lambda2-c0)/conj(z2)+((lambda2-c0)*conj(z1)+conj(z2)*z3)*conj(z3)/(((b0-lambda2)*conj(z2)-conj(z1)*conj(z3))*conj(z2)),((lambda2-c0)*conj(z1)+conj(z2)*z3)/((b0-lambda2)*conj(z2)-conj(z1)*conj(z3)),1+complexarr(ddim)]
  eigvect[*,*,*,2] = [(lambda3-c0)/conj(z2)+((lambda3-c0)*conj(z1)+conj(z2)*z3)*conj(z3)/(((b0-lambda3)*conj(z2)-conj(z1)*conj(z3))*conj(z2)),((lambda3-c0)*conj(z1)+conj(z2)*z3)/((b0-lambda3)*conj(z2)-conj(z1)*conj(z3)),1+complexarr(ddim)]
  
  ; sort eigenvectors according to descending order of eigenvalues
  for i=0,ddim(0)-1 do for j=0,ddim(1)-1 do eigvec[i,j,*,*] = eigvec[i,j,ev_order[i,j,*],*]  
end  

return, eigenvals
end