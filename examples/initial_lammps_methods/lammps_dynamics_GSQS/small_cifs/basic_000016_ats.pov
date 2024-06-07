#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {perspective
  right -7.30*x up 7.57*y
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

cylinder {< -2.75,  -2.36,  -1.84>, < -1.95,  -1.68,  -6.03>, Rcell pigment {Black}}
cylinder {< -2.73,   2.93,  -0.99>, < -1.93,   3.60,  -5.18>, Rcell pigment {Black}}
cylinder {<  2.67,   2.76,   0.01>, <  3.47,   3.43,  -4.18>, Rcell pigment {Black}}
cylinder {<  2.66,  -2.53,  -0.84>, <  3.46,  -1.86,  -5.03>, Rcell pigment {Black}}
cylinder {< -2.75,  -2.36,  -1.84>, < -2.73,   2.93,  -0.99>, Rcell pigment {Black}}
cylinder {< -1.95,  -1.68,  -6.03>, < -1.93,   3.60,  -5.18>, Rcell pigment {Black}}
cylinder {<  3.46,  -1.86,  -5.03>, <  3.47,   3.43,  -4.18>, Rcell pigment {Black}}
cylinder {<  2.66,  -2.53,  -0.84>, <  2.67,   2.76,   0.01>, Rcell pigment {Black}}
cylinder {< -2.75,  -2.36,  -1.84>, <  2.66,  -2.53,  -0.84>, Rcell pigment {Black}}
cylinder {< -1.95,  -1.68,  -6.03>, <  3.46,  -1.86,  -5.03>, Rcell pigment {Black}}
cylinder {< -1.93,   3.60,  -5.18>, <  3.47,   3.43,  -4.18>, Rcell pigment {Black}}
cylinder {< -2.73,   2.93,  -0.99>, <  2.67,   2.76,   0.01>, Rcell pigment {Black}}
atom(< -2.08,   1.14,  -1.84>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #0
atom(<  1.86,  -0.72,  -0.88>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #1
atom(<  1.45,   0.45,  -2.79>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #2
atom(< -0.20,  -2.21,  -1.37>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #3
atom(< -1.83,  -1.31,  -4.03>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #4

// no constraints
