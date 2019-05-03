; plot frequency spectrum of compressed signal
p6a.select
p6a.transparency = 90
p6a.thick=0.1

; add hamming windowed frequency spectrum of compressed signal
p6b = plot(f/1e6, abs(shift(csh_fft,round(N/2))), color='#bb00bb',/OVERPLOT, name = '$\itcs_h\rm(\it\nu\rm)$')
p6b.title = 'Compressed signal $\its\rm(\it\nu\rm)$'
; extract id of title
settitlefontsize, p6b, 10
legend6 = legend(target=[p6a, p6b],position=[lxp[2],lyp[1]-0.16], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)

; add arrow and text
fw = fwhm(f/1e6,abs(shift(csh_fft,round(N/2))),Pos=l6x,yval=l6y)
t6 = text(mean(l6x),l6y[0]*1.1,'$bw=$'+string(format='(%"%4.1f")',fw)+'$ MHz$',/DATA,alignment=0.5, font_size=10,target=p6a)
l6 = marrow(l6x,[l6y[0],l6y[0]],width=1,height=5000,/DATA,target=p6a)
p6a.xrange = [0, 40]