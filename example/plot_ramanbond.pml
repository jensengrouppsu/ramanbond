
bg white
set ray_trace_mode, 3
set ray_trace_color, black
viewport 800, 600
color black, name c
color yellow, name s
set stick_color, white
run ramanbond/plotatombond.py
load Ag8_co2.xyz
plotatombond('Ag8_co2_mode.p')
set_view (     1.000000000,    0.000000000,    0.000000000,     0.000000000,    1.000000000,    0.000000000,     0.000000000,    0.000000000,    1.000000000,     0.000000000,    0.000000000,  -50.000000000,     0.000000000,    0.000000000,    0.000000000,    40.000000000,  100.000000000,  -20.000000000 )
ray 1600, 1200
