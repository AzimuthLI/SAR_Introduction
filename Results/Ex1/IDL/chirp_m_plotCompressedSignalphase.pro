all_phases = [phase_compressed[frange_phi]-phase_compressed[f0I], phase_cs1_fft[frange_phi]-phase_cs1_fft[f0I], phase_cs2_fft[frange_phi]-phase_cs2_fft[f0I]]
yran = [min(all_phases), max(all_phases)]*1.1

p9a = plot(f/1e6,phase_cs1_fft-phase_cs1_fft[f0I], color='#00ff00', yrange = yran, thick=1,$
           /CURRENT, position=[x0g[2], y0g[2],x1g[2],y1g[2]],name='$\phi_{c1}(\nu)$')
p9b = plot(f/1e6, phase_cs2_fft-phase_cs1_fft[f0I], color='#0000ff', name='$\phi_{c2}(\nu)$',/OVERPLOT,thick=1)
p9c = plot(f/1e6, phase_compressed-phase_compressed[f0I], color='#ff0000', name='$\phi_{cs}(\nu)$',/OVERPLOT,thick=0.33,linestyle=2)           
p9a.xtitle='Frequency !9n !3 (MHz)'
p9a.ytitle='Phase !9f!C !3 of Signal (rad)'
p9a.XTICKVALUES = [0, 20, 40]
p9a.xminor = 1
p9a.yminor = 0
p9a.xticklen=0.025
p9a.yticklen=0.025
legend9 = legend(target=[p9a, p9b,p9c],position=[lxp[2],lyp[2]-0.12], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)
p9a.xrange =[-20, 50]