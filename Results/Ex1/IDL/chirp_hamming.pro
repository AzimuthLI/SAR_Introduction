pro chirp_hamming

; define positions of the plotted axes and legends
@chirp_DefinePlotPositions

; call Batchfile to create pulse
@chirp_CreatePulse

; calculate and plot simple pulse as reference
@chirp_Calculations
@chirp_plotalldata

; calculate hamming windowed signal..
@chirp_h_Calculations
@chirp_h_plotalldata

; write result to pdf-file
pdf = 'figure2.pdf'
@chirp_writePlotwindowToFile

cs_ifft_big = shift(abs(cs_ifft) GE 0.5*max(abs(cs_ifft)),round(N/2));
print, 'uncompressed pulse length: ' + string((t(findfirst(cs_ifft_big,/last)) - t(findfirst(cs_ifft_big)))*1e9) + ' ns'
csh_ifft_big = shift(abs(csh_ifft) GE 0.5*max(abs(csh_ifft)),round(N/2));
print, '  compressed pulse length: ' + string((t(findfirst(csh_ifft_big,/last)) - t(findfirst(csh_ifft_big)))*1e9) + ' ns'


stop
end

