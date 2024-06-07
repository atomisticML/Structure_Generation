#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {perspective
  right -6.83*x up 7.62*y
  direction 50.00*z
  location <0,0,50.00> look_at <0,0,0>}


light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}
// no fog
#declare simple = finish {phong 0.7}
#declare pale = finish {ambient 0.5 diffuse 0.85 roughness 0.001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 specular 0.5 }
#declare jmol = finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.05 diffuse 0.3 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.01 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.050;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     translate LOC}
#end

cylinder {< -2.49,  -3.48,  -1.74>, < -1.52,  -2.66,  -6.82>, Rcell pigment {Black}}
cylinder {< -2.48,   2.07,  -0.85>, < -1.51,   2.88,  -5.93>, Rcell pigment {Black}}
cylinder {<  2.16,   1.92,   0.01>, <  3.13,   2.73,  -5.07>, Rcell pigment {Black}}
cylinder {<  2.15,  -3.63,  -0.88>, <  3.12,  -2.81,  -5.96>, Rcell pigment {Black}}
cylinder {< -2.49,  -3.48,  -1.74>, < -2.48,   2.07,  -0.85>, Rcell pigment {Black}}
cylinder {< -1.52,  -2.66,  -6.82>, < -1.51,   2.88,  -5.93>, Rcell pigment {Black}}
cylinder {<  3.12,  -2.81,  -5.96>, <  3.13,   2.73,  -5.07>, Rcell pigment {Black}}
cylinder {<  2.15,  -3.63,  -0.88>, <  2.16,   1.92,   0.01>, Rcell pigment {Black}}
cylinder {< -2.49,  -3.48,  -1.74>, <  2.15,  -3.63,  -0.88>, Rcell pigment {Black}}
cylinder {< -1.52,  -2.66,  -6.82>, <  3.12,  -2.81,  -5.96>, Rcell pigment {Black}}
cylinder {< -1.51,   2.88,  -5.93>, <  3.13,   2.73,  -5.07>, Rcell pigment {Black}}
cylinder {< -2.48,   2.07,  -0.85>, <  2.16,   1.92,   0.01>, Rcell pigment {Black}}
atom(<  1.45,  -1.10,  -1.70>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #0
atom(<  1.33,  -1.25,  -3.89>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #1
atom(< -1.86,  -0.11,  -3.79>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #2
atom(<  0.50,   2.02,  -2.12>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #3
atom(<  1.86,   2.24,  -5.12>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #4

// no constraints
