#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {perspective
  right -5.72*x up 6.76*y
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

cylinder {< -2.20,  -3.08,  -1.50>, < -1.46,  -2.47,  -5.32>, Rcell pigment {Black}}
cylinder {< -2.18,   1.48,  -0.76>, < -1.45,   2.09,  -4.59>, Rcell pigment {Black}}
cylinder {<  1.99,   1.34,   0.01>, <  2.72,   1.95,  -3.81>, Rcell pigment {Black}}
cylinder {<  1.98,  -3.22,  -0.72>, <  2.71,  -2.61,  -4.55>, Rcell pigment {Black}}
cylinder {< -2.20,  -3.08,  -1.50>, < -2.18,   1.48,  -0.76>, Rcell pigment {Black}}
cylinder {< -1.46,  -2.47,  -5.32>, < -1.45,   2.09,  -4.59>, Rcell pigment {Black}}
cylinder {<  2.71,  -2.61,  -4.55>, <  2.72,   1.95,  -3.81>, Rcell pigment {Black}}
cylinder {<  1.98,  -3.22,  -0.72>, <  1.99,   1.34,   0.01>, Rcell pigment {Black}}
cylinder {< -2.20,  -3.08,  -1.50>, <  1.98,  -3.22,  -0.72>, Rcell pigment {Black}}
cylinder {< -1.46,  -2.47,  -5.32>, <  2.71,  -2.61,  -4.55>, Rcell pigment {Black}}
cylinder {< -1.45,   2.09,  -4.59>, <  2.72,   1.95,  -3.81>, Rcell pigment {Black}}
cylinder {< -2.18,   1.48,  -0.76>, <  1.99,   1.34,   0.01>, Rcell pigment {Black}}
atom(< -1.33,   0.68,  -3.87>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #0
atom(<  1.30,  -0.78,  -3.57>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #1
atom(<  0.82,   1.83,  -3.11>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #2
atom(<  1.15,   0.51,  -1.28>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #3
atom(< -0.54,  -1.41,  -3.19>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #4

// no constraints
