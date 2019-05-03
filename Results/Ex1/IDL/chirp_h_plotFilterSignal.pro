; plot filter signal
p2a.select

; make old signal 90% transparent
p2a.TRANSPARENCY = 90
p2b.TRANSPARENCY = 90
p2a.thick=0.1
p2b.thick=0.1

; add hamming windowed filter function
p2c = plot(t*1e6, real(hh), color='#0000ff', linestyle=0,/OVERPLOT, name='real Hamming filter')
p2d = plot(t*1e6, imag(hh), color='#cccc00', linestyle=0,/OVERPLOT, name='imag Hamming filter')
legend2.add, p2c
legend2.add, p2d
legend2.position = legend2.position - [0.056, 0.15]

; add arrow and text
fw = fwhm(t*1e6,abs(hh),Pos=l2x,yval=l2y)
t2 = text(mean(l2x),1.2,'$\tau_p$'+string(format='(%"%4.2f")',fw)+'$ \mus$',/DATA,alignment=0.5, font_size=10,target=p2a)
l2 = marrow(l2x,[1.1,1.1],width=0.05,height=0.05,/DATA,target=p2a)
p2d.yrange = [-1.2, 1.5]
