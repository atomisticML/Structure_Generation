#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {perspective
  right -5.70*x up 6.94*y
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

cylinder {< -2.71,  -2.91,  -1.65>, < -1.66,  -2.02,  -7.18>, Rcell pigment {Black}}
cylinder {< -2.70,   2.42,  -0.80>, < -1.64,   3.31,  -6.33>, Rcell pigment {Black}}
cylinder {<  1.66,   2.28,   0.01>, <  2.71,   3.17,  -5.52>, Rcell pigment {Black}}
cylinder {<  1.64,  -3.05,  -0.84>, <  2.70,  -2.16,  -6.37>, Rcell pigment {Black}}
cylinder {< -2.71,  -2.91,  -1.65>, < -2.70,   2.42,  -0.80>, Rcell pigment {Black}}
cylinder {< -1.66,  -2.02,  -7.18>, < -1.64,   3.31,  -6.33>, Rcell pigment {Black}}
cylinder {<  2.70,  -2.16,  -6.37>, <  2.71,   3.17,  -5.52>, Rcell pigment {Black}}
cylinder {<  1.64,  -3.05,  -0.84>, <  1.66,   2.28,   0.01>, Rcell pigment {Black}}
cylinder {< -2.71,  -2.91,  -1.65>, <  1.64,  -3.05,  -0.84>, Rcell pigment {Black}}
cylinder {< -1.66,  -2.02,  -7.18>, <  2.70,  -2.16,  -6.37>, Rcell pigment {Black}}
cylinder {< -1.64,   3.31,  -6.33>, <  2.71,   3.17,  -5.52>, Rcell pigment {Black}}
cylinder {< -2.70,   2.42,  -0.80>, <  1.66,   2.28,   0.01>, Rcell pigment {Black}}
atom(<  0.07,   1.92,  -5.56>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #0
atom(< -0.40,  -1.92,  -1.85>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #1
atom(< -0.23,   0.49,  -4.06>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #2
atom(<  0.69,   1.84,  -2.27>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #3
atom(<  0.60,  -0.24,  -2.28>, 1.39, rgb <0.75, 0.75, 0.90>, 0.0, ase2) // #4

// no constraints
