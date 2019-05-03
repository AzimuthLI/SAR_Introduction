all_phases = [phase_h1_fft[where(f gt 0 and f le 50e6)]-phase_h1_fft[f00I],phase_h2_fft[where(f le 0 and f gt -20e6)]-phase_h2_fft[f00I],phase_filter[frange_phi]-phase_filter[f0I]]
yran = [min(all_phases), max(all_phases)]

p8a = plot(f[where(f gt 0)]/1e6,phase_h1_fft[where(f gt 0)]-phase_h1_fft[f00I], color='#00ff00', Linestyle=0, yrange = yran, $
           /CURRENT, position=[x0g[1], y0g[2],x1g[1],y1g[2]],name='$\phi_{h1}(\nu)$')
p8b = plot(f[where(f le 0)]/1e6, phase_h2_fft[where(f le 0)]-phase_h2_fft[f00I], color='#0000ff', name='$\phi_{h2}(\nu)$',/OVERPLOT)
p8c = plot(f/1e6, phase_filter-phase_filter[f0I], color='#ff0000', name='$\phi_h(\nu)$',/OVERPLOT)           
p8a.xtitle='Frequency !9n !3 (MHz)'
p8a.ytitle='Phase !9f!C !3 of Signal (rad)'
p8a.XTICKVALUES = [0, 20, 40]
p8a.xminor = 1
p8a.yminor = 0
p8a.xticklen=0.025
p8a.yticklen=0.025
legend8 = legend(target=[p8a, p8b,p8c],position=[lxp[1],lyp[2]-0.12], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)
p8a.xrange=[-20, 50]