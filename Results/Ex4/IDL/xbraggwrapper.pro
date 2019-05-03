function Xbraggwrapper, x
  epsilon = transpose([x(0,*)])
  beta1   = transpose([x(1,*)])
  H0      = scope_varfetch('H0_fetch',Level=-1)
  alpha0  = scope_varfetch('alpha0_fetch',Level=-1)
  theta_aoi = scope_varfetch('theta_aoi_fetch',Level=-1)

  res = Xbragg(transpose([x(0,*)]),transpose([x(1,*)]),theta_aoi)
  return, [res[0:n_elements(epsilon)-1,*]-H0,res[n_elements(epsilon):(n_elements(epsilon)*2-1),*]-alpha0]  
end