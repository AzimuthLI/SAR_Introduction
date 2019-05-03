  function  matcharr, dataA, dataB
    dimA = size(dataA,/dimensions)
    dimB = size(dataB,/dimensions)
    atype = size(A,/type)
    btype = size(B,/type)
    ; put both data into a array of same size:
    ddim = [max([dimA[0],dimB[0]],xmI), max([dimA[1],dimB[1]],ymI)]
    data = fltarr(ddim)
    ;data2 = fltarr(ddim)
    data[0:dimA[0]-1,0:dimA[1]-1] = abs(dataA)
    ;data2[0:dimB[0]-1,0:dimB[1]-1] = abs(dataB)
    return, data
  end
