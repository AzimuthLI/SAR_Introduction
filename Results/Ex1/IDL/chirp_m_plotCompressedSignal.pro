p3a = plot(t*1e6, abs(cs1), color='#00ff00', /CURRENT, position=[x0g[2], y0g[0],x1g[2],y1g[0]],name = '$|cs_1(t)|$',thick=1.5) 
p3b = plot(t*1e6, abs(cs2), color='#0000ff',/OVERPLOT, name='$|cs_2(t)|$',thick=1)
p3c = plot(t*1e6, abs(cs)/2,  color='#ff0000',/OVERPLOT,name='$|cs(t)|$')
p3d = plot(t*1e6, (abs(cs1)+abs(cs2))/2,  color='#ffff00',/OVERPLOT,name='$|cs_{m=2}(t)|$',thick=0.33)
p3e = plot(t*1e6, abs(cs_mlook_ifft)/max((cs_mlook_ifft))*max(abs(cs1)),  color='#ffbb00',/OVERPLOT,name='$|cs_{m=2}(t),fft|$')

p3a.xtitle = 'time !5t!3(!9m!3s)'
p3a.ytitle = 'Absolut Value !Cof Compressed Signal !5cs!3(!5t!3) !C (a.u.)'
p3a.XTICKVALUES = [-0.2,0, 0.2]
p3a.YTICKVALUES = [0, 1000, 2000] 
p3a.xminor = 1
p3a.yminor = 0
p3a.xticklen=0.025
p3a.yticklen=0.025
legend3 = legend(target=[p3a, p3b, p3c, p3d,p3e],position=[lxp[2],lyp[0]], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)


cs_mlook_ifft_I  = abs(cs_mlook_ifft) GE 0.5*max(abs(cs_ifft));
cs_mlook_tau = (t(findfirst(cs_mlook_ifft_I,/last)) - t(findfirst(cs_mlook_ifft_I)))
print, '     multi look compressed pulse length: ' + strtrim(string(cs_mlook_tau*1e9),2) + ' ns (resolution = '+strtrim(string(cs_mlook_tau*c),2) + 'm)'

fwm = fwhm(t*1e6, abs(cs_mlook_ifft)/max((cs_mlook_ifft))*max(abs(cs1)), pos=l3xa,yval=l3ya) 
fws1 = fwhm(t*1e6, abs(cs)/2, pos=l3xb,yval=l3yb)
 
t3a = text(-0.01,l3ya[0]*1.2,'$\tau_p=$'+string(format='(%"%4.2f")',fwm)+'$ \mus$',/DATA,alignment=1, font_size=10,target=p3a)
l3a = marrow(l3xa,l3ya*1.1,width=0.01,height=50,/DATA,target=p3a)

t3b = text(mean(l3xb),l3yb[0]*0.7,'$\tau_p=$'+string(format='(%"%4.2f")',fws1)+'$ \mus$',/DATA,alignment=1, font_size=10,target=p3a)
l3b = marrow(l3xb,l3yb*0.9,width=0.01,height=50,/DATA,target=p3a)

p3a.xrange = [-0.2, 0.3]
p3a.yrange = [0, 2500]