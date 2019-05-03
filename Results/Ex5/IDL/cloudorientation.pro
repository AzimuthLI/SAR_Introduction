function cloudorientation, c, significance=significance
; c: complex data set cloud

  N = 8; not be lower than 4 
  
  ; rotate dataset in N steps by pi
  cx = fltarr(N)
  phi = findgen(N)/N*!pi
  for j=0,N-1 do begin
    crx   = real(c*exp(complex(0,1)*phi[j]))    
    cx[j] = sqrt(total((crx-mean(crx))^2))
  end
  r0 = mean(cx)
  ; project cx(phi) on cos-sin-orthonormalsystem   
  c1 = total((cx-r0)*cos(2*phi))/(0.5*N)
  c2 = total((cx-r0)*sin(2*phi))/(0.5*N)   
  ; calculate angle of long half axis against x-axis
  ;phix = +0.5*atan(c1,c2)-!pi/4
  phix = -0.5*atan(c2,c1)
  
  ; calculate minima and maxia width of rotation 
  crx = real(c*exp(complex(0,1)*(!pi/2-phix)))  
  b   = sqrt(total((crx-mean(crx))^2))
  crx = real(c*exp(complex(0,1)*(-phix)))  
  a   = sqrt(total((crx-mean(crx))^2))
  if max(cx) gt a*1.01 then print, 'max(cx) > cxmax: ', max(cx), a 
  if min(cx) lt b*0.99 then print, 'max(cx) > cxmax: ', min(cx), b

  significance = (a-b)/(a+b)
  
;  p = plot(phi,cx,'b',/o)
;  p = plot(phi, r0+c1*cos(2*phi)+c2*sin(2*phi),':b',/o)
;  q = plot(real(c),imag(c),'*k',aspect_ratio=1)
;  q = plot(mean(real(c))+[0,r0*cos(phix)],mean(imag(c))+[0,r0*sin(phix)],'b',/o,axis_style=0)
  return, phix
end
