function matchwrapper, scale, return_img=return_img
  dataA = scope_varfetch('dataA',Level=-1)
  dataB = scope_varfetch('dataB',Level=-1)

  dimA  = size(dataA,/dimensions)
  dataA = congrid(dataA,round(dimA[0]*scale),round(dimA[1]*scale))
     
  data1 = matcharr(dataA, dataB)
  data2 = matcharr(dataB, dataA)
  c = corr(data1,data2)
  c_max     = max(abs(c),location)
  ; find coordinates of maximum of correlation function
  arrIndices = array_indices(c,location)
  
  if keyword_set(return_img) then begin     
    (scope_varfetch('deltaxy',Level=-1,/ENTER)) = arrindices
    (scope_varfetch('data1',Level=-1,/ENTER)) = data1
    (scope_varfetch('data2',Level=-1,/ENTER)) = data2
  end
  
  return, c_max
  
end