p3a = plot(t*1e6, real(cs), color='#bb00bb', yrange=[min(real(cs)),max(real(cs))], $
     /CURRENT, position=[x0g[2], y0g[0],x1g[2],y1g[0]],name = '\phi(\nu)') 
p3b = plot(t*1e6, imag(cs), color='#88aa00',/OVERPLOT)
p3c = plot(t*1e6, abs(cs),  color='#000000',/OVERPLOT, linestyle=0)
p3a.xtitle = 'time !5t!3(!9m!3s)'
p3a.ytitle = 'Absolut Value !Cof Compressed Signal !5cs!3(!5t!3) !C (a.u.)'
p3a.name   = 'real cs(t) ' 
p3a.XTICKVALUES = [-0.2,0, 0.2]
p3a.YTICKVALUES = [-2000, 0, 2000, 4000] 
p3a.xminor = 1
p3a.yminor = 0
p3a.xticklen=0.025
p3a.yticklen=0.025
p3b.name   = 'imag !5cs!3(!5t!3)'
p3c.name   = 'abs !5cs!3(!5t!3)'
legend3 = legend(target=[p3a, p3b, p3c],position=[lxp[2],lyp[0]], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)

