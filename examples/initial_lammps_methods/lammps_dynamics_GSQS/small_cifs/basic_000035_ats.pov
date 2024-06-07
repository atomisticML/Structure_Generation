#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {perspective
  right -6.98*x up 6.27*y
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

cylinder {< -2.00,  -2.28,  -1.49>, < -0.98,  -1.42,  -6.86>, Rcell pigment {Black}}
cylinder {< -1.99,   2.12,  -0.78>, < -0.96,   2.98,  -6.16>, Rcell pigment {Black}}
cylinder {<  2.30,   1.98,   0.01>, <  3.32,   2.85,  -5.36>, Rcell pigment {Black}}
cylinder {<  2.29,  -2.42,  -0.69>, <  3.31,  -1.56,  -6.07>, Rcell pigment {Black}}
cylinder {< -2.00,  -2.28,  -1.49>, < -1.99,   2.12,  -0.78>, Rcell pigment {Black}}
cylinder {< -0.98,  -1.42,  -6.86>, < -0.96,   2.98,  -6.16>, Rcell pigment {Black}}
cylinder {<  3.31,  -1.56,  -6.07>, <  3.32,   2.85,  -5.36>, Rcell pigment {Black}}
cylinder {<  2.29,  -2.42,  -0.69>, <  2.30,   1.98,   0.01>, Rcell pigment {Black}}
cylinder {< -2.00,  -2.28,  -1.49>, <  2.29,  -2.42,  -0.69>, Rcell pigment {Black}}
cylinder {< -0.98,  -1.42,  -6.86>, <  3.31,  -1.56,  -6.07>, Rcell pigment {Black}}
cylinder {< -0.96,   2.98,  -6.16>, <  3.32,   2.85,  -5.36>, Rcell pigment {Black}}
cylinder {< -1.99,   2.12,  -0.78>, <  2.30,   1.98,   0.01>, Rcell pigment {Black}}
atom(< -1.93,   0.37,  -1.14>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #0
atom(< -1.05,  -0.16,  -3.93>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #1
atom(<  0.71,   0.23,  -5.72>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #2
atom(<  0.74,  -1.25,  -2.07>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #3
atom(<  0.35,  -1.59,  -4.88>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #4

// no constraints
