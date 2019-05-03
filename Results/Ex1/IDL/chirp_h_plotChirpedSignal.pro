; plot chirped signal
p1a.select

; add hamming window / envelope
p1c = plot(t*1e6, hammingw, color='#0000ff', linestyle=0,/OVERPLOT,name='Hamming Window')
p1c.YRANGE = [-1.2, 1.5]
legend1.add, p1c
legend1.position = legend1.position + [-0.04, 0.02]