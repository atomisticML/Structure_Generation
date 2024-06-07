#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {perspective
  right -6.14*x up 6.97*y
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

cylinder {< -2.48,  -3.20,  -1.43>, < -1.65,  -2.50,  -5.77>, Rcell pigment {Black}}
cylinder {< -2.47,   1.45,  -0.69>, < -1.64,   2.14,  -5.03>, Rcell pigment {Black}}
cylinder {<  1.29,   1.33,   0.01>, <  2.12,   2.02,  -4.33>, Rcell pigment {Black}}
cylinder {<  1.28,  -3.32,  -0.73>, <  2.11,  -2.63,  -5.07>, Rcell pigment {Black}}
cylinder {< -2.48,  -3.20,  -1.43>, < -2.47,   1.45,  -0.69>, Rcell pigment {Black}}
cylinder {< -1.65,  -2.50,  -5.77>, < -1.64,   2.14,  -5.03>, Rcell pigment {Black}}
cylinder {<  2.11,  -2.63,  -5.07>, <  2.12,   2.02,  -4.33>, Rcell pigment {Black}}
cylinder {<  1.28,  -3.32,  -0.73>, <  1.29,   1.33,   0.01>, Rcell pigment {Black}}
cylinder {< -2.48,  -3.20,  -1.43>, <  1.28,  -3.32,  -0.73>, Rcell pigment {Black}}
cylinder {< -1.65,  -2.50,  -5.77>, <  2.11,  -2.63,  -5.07>, Rcell pigment {Black}}
cylinder {< -1.64,   2.14,  -5.03>, <  2.12,   2.02,  -4.33>, Rcell pigment {Black}}
cylinder {< -2.47,   1.45,  -0.69>, <  1.29,   1.33,   0.01>, Rcell pigment {Black}}
atom(<  0.05,   0.10,  -3.74>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #0
atom(< -1.53,  -0.55,  -4.75>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #1
atom(<  1.53,   1.93,  -4.32>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #2
atom(< -0.51,  -1.36,  -1.88>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #3
atom(< -1.41,   0.75,  -2.13>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #4

// no constraints
