#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {perspective
  right -5.60*x up 7.70*y
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

cylinder {< -2.47,  -2.73,  -1.53>, < -1.67,  -2.07,  -5.68>, Rcell pigment {Black}}
cylinder {< -2.45,   2.66,  -0.66>, < -1.66,   3.33,  -4.81>, Rcell pigment {Black}}
cylinder {<  1.15,   2.55,   0.01>, <  1.95,   3.21,  -4.14>, Rcell pigment {Black}}
cylinder {<  1.14,  -2.85,  -0.86>, <  1.93,  -2.19,  -5.01>, Rcell pigment {Black}}
cylinder {< -2.47,  -2.73,  -1.53>, < -2.45,   2.66,  -0.66>, Rcell pigment {Black}}
cylinder {< -1.67,  -2.07,  -5.68>, < -1.66,   3.33,  -4.81>, Rcell pigment {Black}}
cylinder {<  1.93,  -2.19,  -5.01>, <  1.95,   3.21,  -4.14>, Rcell pigment {Black}}
cylinder {<  1.14,  -2.85,  -0.86>, <  1.15,   2.55,   0.01>, Rcell pigment {Black}}
cylinder {< -2.47,  -2.73,  -1.53>, <  1.14,  -2.85,  -0.86>, Rcell pigment {Black}}
cylinder {< -1.67,  -2.07,  -5.68>, <  1.93,  -2.19,  -5.01>, Rcell pigment {Black}}
cylinder {< -1.66,   3.33,  -4.81>, <  1.95,   3.21,  -4.14>, Rcell pigment {Black}}
cylinder {< -2.45,   2.66,  -0.66>, <  1.15,   2.55,   0.01>, Rcell pigment {Black}}
atom(<  0.76,   2.16,  -4.39>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #0
atom(< -0.25,   2.28,  -2.70>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #1
atom(<  1.16,  -2.28,  -2.17>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #2
atom(< -1.28,   0.51,  -4.66>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #3
atom(<  1.28,  -0.36,  -1.84>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #4

// no constraints
