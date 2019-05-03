function xbragg_boundary, rla         
; create boundary where inversion is possible
Nangles = 100
Nbeta   = 30
Neps    = 30
H_high     = fltarr(Nangles,Nbeta)
H_low      = fltarr(Nangles,Neps)
alpha_high = fltarr(Nangles,Nbeta)
alpha_low  = fltarr(Nangles,Neps)

beta1_bd = (findgen(Nbeta))/(Nbeta-1)*!pi/2
eps_bd   =  ((findgen(Neps)+1)/Neps)^2*100 + 1

for i=0,Nangles-1 do begin
  theta = float(i)/(Nangles-1)*(max(rla)-min(rla)) + min(rla)   
  res = Xbragg([100], beta1_bd, theta)
  H_high[i,*] = res[0,*]
  alpha_high[i,*] = res[1,*] 
  res = Xbragg(eps_bd, [!pi/2], theta)
  
  H_low[i,*] = res[0:Neps-1]
  alpha_low[i,*] = res[Neps:2*Neps-1]
endfor
      
return, {alpha_high:alpha_high, H_high:H_high, alpha_low:alpha_low, H_low:H_low}
end