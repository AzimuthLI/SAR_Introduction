; Exercise 1.1
; define procedure chirpP (same name as file-name.)
pro chirpSimplePulse

; define positions of the plotted axes and legends
@chirp_DefinePlotPositions

; call Batchfile to create pulse
@chirp_CreatePulse

; do all calculations
@chirp_Calculations


;----------------------------------------
; PLOT all Simple Pulse Data
;----------------------------------------
@chirp_plotalldata

; add arrow and text
fw  = fwhm(t*1e6,abs(cs),Pos=l3x,yval=l3y)
t3 = text(-0.05,l3y[1],'$\tau_p = $'+string(format='(%"%4.2f")',fw)+'$ \mus$',/DATA,alignment=1, font_size=10,target=p3a)
l3 = marrow(l3x,l3y,width=0.01,height=100,/DATA,target=p3a)
p3a.xrange = [-0.25,0.25]

; write plot window into a pdf file
pdf = 'figure1.pdf'
@chirp_writePlotwindowToFile
stop
end