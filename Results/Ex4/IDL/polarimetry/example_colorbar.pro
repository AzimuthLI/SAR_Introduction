@colorbar_define

pro example_colorbar

array=fltarr(100,200)

array=findgen(100)#findgen(200) 

loadct,39

device, decompose=0
WINDOW,1,xsize=200+320,ysize=500,title='Title'
tv,congrid(bytscl(array,0,20000),200,500)
colorbar = Obj_New("COLORBAR", Format='(F8.1)',/vertical)
colorbar->Draw
colorbar->SetProperty, Range=[0,20000], /Erase, /Draw
Obj_Destroy, colorbar
!p.background=255
device, decompose=1

end
