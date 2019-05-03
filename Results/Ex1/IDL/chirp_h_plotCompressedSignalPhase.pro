p9a.select
p9b = plot(f/1e6, phase_csh - phase_csh[f0I], color='#ffbb00',name = '$\phi_{ch}(\nu)$',/OVERPLOT,tick=0.33)
legend9 = legend(target=[p9a, p9b],position=[lxp[2],lyp[2]-0.12], font_size=8, sample_width=0.03, HORIZONTAL_ALIGNMENT='right', vertical_spacing=0.005, shadow=0)
p9b.yrange=[-0.045,0.06]
p9b.ytickvalues=[-0.04,0,0.04]
