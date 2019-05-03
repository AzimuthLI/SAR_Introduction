function plotucirc, data, color, _EXTRA = extras
if n_params() eq 2 then begin  
  LOADCT, 1, NCOLORS=256, RGB_TABLE=ctab
  case color of 
  'r': c = shift(ctab,0,1)
  'g': c = shift(ctab,0,2)
  'b': c = shift(ctab,0,0)
  end
end else begin
  c = 0
  color = 'w'
end
;stop
N = n_elements(data)
phi = findgen(48)/47*2*!pi
h = hist_2D(real(reform(data,N)),imag(reform(data,N)),min1=-1,min2=-1,max1=1,max2=1,bin1=3./sqrt(N),bin2=3./sqrt(N))
dimh = size(h,/dim)
x = findgen(dimh[0])/(dimh[0]-1)*2-1
y = findgen(dimh[1])/(dimh[1]-1)*2-1
p = image(alog(h+0.01),x,y, _EXTRA = extras, rgb_table=c)
p1 = plot(1.01*cos(phi),1.01*sin(phi),color,AXIS_STYLE=0, _EXTRA = extras,title='',overplot=p)
;p = plot([-1,1],[0,0],'-w1',AXIS_STYLE=0,overplot=p1)
;p = plot([0,0],[-1,1],'-w1',AXIS_STYLE=0,overplot=p1)
return, p
end