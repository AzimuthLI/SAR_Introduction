; plot frequency spectrum of signal 
p5a.select

; make old spectrum transparent and add name
p5a.transparency = 90
p5a.thick=0.1
p5a.name = '!5h!3(!9n!3)'

; add hamming windowed signal
p5b = plot(f/1e6, abs(shift(hh_fft,round(N/2))), color='#0000ff', name='filter $\ith\rm(\nu)\cdot\itw\rm(\nu)$',/OVERPLOT)
p5c = plot(f/1e6, (hamming_fft)/max(hamming_fft)*max(abs(h_fft)), color='#ff8800', name='hamming $\itw\rm(\nu)$',/OVERPLOT)
legend5 = legend(target=[p5a, p5b, p5c],position=[lxp[1],lyp[1]-0.12], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)

; add arrow and text
fw = fwhm(f/1e6,abs(shift(hh_fft,round(N/2))),Pos=l5x,yval=l5y)
t5 = text(mean(l5x),l5y[1],'$bw=$'+string(format='(%"%4.1f")',fw)+'$ MHz$',/DATA,alignment=0.5, font_size=10,target=p5a)
l5 = marrow(l5x,[l5y[0],l5y[0]],width=1,height=40,/DATA,target=p5a)
p5a.xrange = [0, 40]
