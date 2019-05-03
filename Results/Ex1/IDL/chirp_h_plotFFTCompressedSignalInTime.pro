p3a.select
; set the transparency of the two old abs-values
p3d.transparency = 90
p3e.transparency = 90
p3a.thick=0.1
p3b.thick=0.1
p3c.thick=0.1
p3d.thick=0.1

; add new ifft pulses
p3f = plot(t*1e6, abs(shift(csh_ifft,round(N/2)))/max(abs(csh_ifft))*max(cs), color='#00aa00',thick=0.33, name = '|csh(t), fft|',/OVERPLOT)
legend3.add, p3f

; add arrow and text
fw = fwhm(t*1e6,abs(shift(csh_ifft,round(N/2)))/max(abs(csh_ifft))*max(cs),Pos=l3x,yval=l3y)
t3 = text(l3x[0]*1.3,l3y[1],'$\tau_p = $'+string(format='(%"%4.2f")',fw)+'$ \mus$',/DATA,alignment=1, font_size=10,target=p3a)
l3 = marrow(l3x,l3y,width=0.01,height=50,/DATA,target=p3a)
p3a.xrange = [-0.25, 0.25]