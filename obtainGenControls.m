function [m,f]=obtainGenControls( delta,omega,e,vg,thetag,pg,qg,OMEGA_S)

%OBTAINGENCONTROLS obtains the steady-states value  of  generator controls.
%  [m0,f0]=obtainGenControls( d_vec,xd_vec,...
%   xprime_vec,delta0,omega0,e0,thetag0,vg0,pg0,qg0,OMEGA_S)
%  obtains the steady-state values of generator mechanical power input $m$
% and generator internal field voltage $f$, accroding to the equilibrium point
% of equations  (1b), (1c) in CDC 2016 paper. 
% The power flow solutions of theta, v, pg, and generator steady-states
% delta, omega, and e must be available.
%
% Description of Outputs:
% 1. m:  generator mechanical power input,  size(G,1)
% 2. f: generator internal field voltage in pu voltage, size(G,1)
%
% Description of inputs:
% 1. d_vec: vector of damping coefficients, size(G,1)
% 2. xd_vec: vector of direct axis synchronous reactance pu, size(G,1)
% 3. xprime_vec: direct axis transient reactance pu, size(G,1)
% 4. delta: vector of generator internal angles in radians, size(G,1)
% 5. omega: vector of generator internal frequencies in radians per second,
% size(G,1)
% 6. e: vector of generator electromotive force in pu volts, size(G,1)
% 7. thetag: vector of generator terminal voltage angles in radians, size(G,1)
% 4. vg: vector of generator terminal voltage magnitudes in pu volts, size(G,1)
% 5. pg: vector of generator real power generated output in pu watts, size(G,1)
% 6. qg: vector of generator reactive power generated output in pu vars, size(G,1)
% 7. OMEGA_S: the steady-state frequency
% 
% See also: obtainGenStates
%
% Required: 
% 1. Fix equation numbers. 

global d_vec xd_vec    xprime_vec

m=d_vec.*(omega-OMEGA_S)+pg; 
f=(xd_vec./xprime_vec).*e- ( (xd_vec-xprime_vec)./(xprime_vec)).*vg.*cos(delta-thetag);