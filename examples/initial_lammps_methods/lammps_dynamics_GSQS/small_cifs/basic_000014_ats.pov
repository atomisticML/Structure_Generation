#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {perspective
  right -9.15*x up 7.16*y
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

cylinder {< -3.79,  -2.68,  -2.14>, < -2.75,  -1.81,  -7.58>, Rcell pigment {Black}}
cylinder {< -3.78,   2.54,  -1.31>, < -2.74,   3.41,  -6.74>, Rcell pigment {Black}}
cylinder {<  3.32,   2.31,   0.01>, <  4.36,   3.18,  -5.42>, Rcell pigment {Black}}
cylinder {<  3.30,  -2.91,  -0.83>, <  4.34,  -2.04,  -6.26>, Rcell pigment {Black}}
cylinder {< -3.79,  -2.68,  -2.14>, < -3.78,   2.54,  -1.31>, Rcell pigment {Black}}
cylinder {< -2.75,  -1.81,  -7.58>, < -2.74,   3.41,  -6.74>, Rcell pigment {Black}}
cylinder {<  4.34,  -2.04,  -6.26>, <  4.36,   3.18,  -5.42>, Rcell pigment {Black}}
cylinder {<  3.30,  -2.91,  -0.83>, <  3.32,   2.31,   0.01>, Rcell pigment {Black}}
cylinder {< -3.79,  -2.68,  -2.14>, <  3.30,  -2.91,  -0.83>, Rcell pigment {Black}}
cylinder {< -2.75,  -1.81,  -7.58>, <  4.34,  -2.04,  -6.26>, Rcell pigment {Black}}
cylinder {< -2.74,   3.41,  -6.74>, <  4.36,   3.18,  -5.42>, Rcell pigment {Black}}
cylinder {< -3.78,   2.54,  -1.31>, <  3.32,   2.31,   0.01>, Rcell pigment {Black}}
atom(<  2.39,   0.58,  -2.27>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #0
atom(<  2.28,  -2.02,  -1.36>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #1
atom(< -2.97,   1.85,  -4.81>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #2
atom(<  0.30,   1.57,  -1.71>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #3
atom(< -1.90,  -0.23,  -6.95>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #4

// no constraints
