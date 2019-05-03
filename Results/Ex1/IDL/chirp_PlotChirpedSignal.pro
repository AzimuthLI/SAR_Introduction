p1a = plot(t*1e6, real(s), color='#ff0000', location=[0,0], $
      dimensions=[1000, 700], position=[x0g[0], y0g[0],x1g[0],y1g[0]])
p1b = plot(t*1e6, imag(s), color='#00ff00', linestyle=0,/OVERPLOT)
p1a.xtitle = 'time !5t!3 (!9m!3s)'
p1a.ytitle = 'Amplitude of Signal !5s!3(!5t!3) !C (a.u.)'
p1a.name   = 'real s(t) '
p1a.XTICKVALUES = [-1,0, 1]
p1a.YTICKVALUES = [-1,0,   1]
 p1a.xticklen=0.025
p1a.yticklen=0.025
p1a.xminor = 4
p1a.yminor = 0
p1b.name   = 'imag s(t) '
legend1 = legend(target=[p1a, p1b],position=[lxp[0],lyp[0]-0.19], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)


fw = fwhm(t*1e6,abs(s),Pos=l1x,yval=l1y)
;l1 = polyline(l1x,[1.1,1.1],arrow_style=3,/DATA,thick=1,arrow_size=0.015)
t1 = text(mean(l1x),1.2,'$\tau_p = $'+string(format='(%"%4.2f")',fw)+'$ \mus$',/DATA,alignment=0.5, font_size=10)
p1a.yrange = [-1.2, 1.5]