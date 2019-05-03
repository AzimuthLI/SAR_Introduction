all_phases = [phase_s1_fft[where(f gt 0 and f le 50e6)]-phase_s1_fft[f00I],phase_s2_fft[where(f le 0 and f gt -20e6)]-phase_s2_fft[f00I],phase_signal[frange_phi]-phase_signal[f0I]]
yran = [min(all_phases), max(all_phases)]

p7a = plot(f[where(f gt 0)]/1e6,phase_s1_fft[where(f gt 0)]-phase_s1_fft[f00I], color='#00ff00', Linestyle=0, yrange = yran, $
           /CURRENT, position=[x0g[0], y0g[2],x1g[0],y1g[2]],name='$\phi_1(\nu)$')
p7b = plot(f[where(f le 0)]/1e6, phase_s2_fft[where(f le 0)]-phase_s2_fft[f00I], color='#0000ff', name='$\phi_2(\nu)$',/OVERPLOT)
p7c = plot(f/1e6, phase_signal-phase_signal[f0I], color='#ff0000', name='$\phi(\nu)$',/OVERPLOT)           
p7a.xtitle='Frequency !9n !3 (MHz)'
p7a.ytitle='Phase !9f!C !3 of Signal (rad)'
p7a.XTICKVALUES = [0, 20, 40]
p7a.xminor = 1
p7a.yminor = 0
p7a.xticklen=0.025
p7a.yticklen=0.025
legend7 = legend(target=[p7a, p7b,p7c],position=[lxp[0],lyp[2]-0.12], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)
p7a.xrange=[-20, 50]