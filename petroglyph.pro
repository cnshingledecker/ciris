PRO petroglyph
file = '/home/cns/losalamos/abundance.csv'


limit = ','

big_data = READ_CSV(file, N_TABLE_HEADER=2, DELIMITER=limit)

x = 1.0e-4 ; cm
y = 1.0e-5 ; cm
z = 1.0e-5 ; cm

volume = x*y*z 

time = big_data.FIELD1
fluence = big_data.FIELD2
sp1 = big_data.FIELD3
sp2 = big_data.FIELD4

n1 = sp1/(1E20*volume) ;number density of O atoms
n2 = sp2/(1e20*volume) ;number density of O3 atoms
 ; fit = curvefit(fluence,n2)
;print, fit


PRINT, time[0],fluence[0],sp1[0],sp2[0]

aspect_ratio=1.5
xsize=9
ysize=xsize/aspect_ratio
set_plot, 'ps'
device,filename='O3_new.eps',encapsulated=1
device, xsize=xsize,ysize=ysize

plot1 = PLOT(fluence,n2,/XLOG, XTITLE='Fluence (H$^{+}$/cm$^{2}$)',YTITLE='[O$_{3}$](10$^{20}$ cm$^{-3}$)',FONT_NAME="Hershey 3",FONT_SIZE=30)
;plot2 = PLOT(SMOOTH(n2,7),/OVERPLOT)
;plot2 = PLOT(time,sp1, /XLOG, /OVERPLOT)

leg = LEGEND(TARGET=[plot1] , /AUTO_TEXT_COLOR,label='O$_{3}$', /RELATIVE, POSITION=[0.3,0.9],FONT_NAME="Hershey 3",FONT_SIZE=30)
device, /close
set_plot, 'x'

;PRINT, sp1[0]




END
