clear;
close all;
clc;

b      = 40;
a0_t   = 2*pi;
a0_r   = 2*pi;
c_t    = 6;
c_r    = 6;
aero_t = 0;
aero_r = 0;
geo_t  = 7;
geo_r  = 7;
N      = 10;

[e, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N)
