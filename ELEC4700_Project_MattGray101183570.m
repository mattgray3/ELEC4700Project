% Milestone 1 Matt Gray 101183570 ELEC 4700 Project 
set(0,'defaulttaxesfontsize', 20)
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'DefaultLineLineWidth', 2)
set(0, 'Defaulttaxeslinewidth', 2)

set(0, 'DefaultFigureWindowStyle', 'docked')

c_c = 299792458;
c_eps_0 = 8.8542149e-12;
c_eps_0_cm = c_eps_0/100; 
c_mu_0 = 1/c_eps_0/c_c^2;
c_q = 1.60217653e-19;
c_hb = 1.05457266913e-34;
c_h = c_hb*2*pi;

InputParasL.E0 = 1e5;
InputParasL.we = 0;
InputParasL.t0 = 2e-12;
InputParasL.wg = 5e-13;
InputParasL.phi = 0;
InputParasR = 0;

n_g = 3.5;
vg = c_c/n_g * 1e2;
lambda = 1550e-9;



