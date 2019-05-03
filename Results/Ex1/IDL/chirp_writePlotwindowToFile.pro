error_status = 0

Catch, error_status
cd, current=pwd 

if (error_status eq 0) then p1a.save, pdf, /LANDSCAPE, /CENTIMETERS, PAGE_SIZE=[21,29.7], width=29.7
if (error_status eq 0) then print, 'graphics written to file '''+pwd+'\'''+pdf+'''.' $   
else begin 
  catch, /cancel
  print, 'error in writting file.'
end