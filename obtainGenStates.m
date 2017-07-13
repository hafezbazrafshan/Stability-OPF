function [ delta, e] = obtainGenStates(vg, thetag, pg,qg )
%OBTAINGENSTATES obtains the steady-states value  of  generator states
%  [ delta, e] = obtainGenStates(xq_vec,xprime_vec, theta, v, pg,qg ) 
%  obtains the steady-state values of generator internal angle delta 
% and generator electromotive force e, accroding to the equilibrium point
% of equations  (2a), (2b) in CDC 2016 paper. 
% The power flow solutions of theta, v, pg, and qg must be available.
%
% Description of Outputs:
% 1. delta:  generators internal phase angle in radians,  size(G,1)
% 2. e: generators electromotive force in pu voltage, size(G,1)
%
% Description of Inputs:
% 1. xq_vec: vector of quadrature axis synchronous reactance pu, size(G,1)
% 2. xprime_vec: vector of  direct axis transient reactance pu, size(G,1)
% 3. theta: vector of terminal voltage angles in radians, size(G,1)
% 4. v: vector of terminal voltage magnitudes in pu voltage, size(G,1)
% 5. pg: vector of real power generated output in pu watts, size(G,1)
% 6. qg: vector of reactive power generated output in pu vars, size(G,1)
%
% See also obtainGenControls
%
% Required: 
% 1. Fix equation numbers. 

global G
global xq_vec xprime_vec
delta=zeros(G,1);
e=zeros(G,1);
for i=1:G  
    x(1) = delta(1) ; x(2) = e(1);
   problem.options = optimoptions('fsolve','Display','off');
   problem.objective= @(x) genStateFunctions( x,xq_vec(i), xprime_vec(i), vg(i),thetag(i), pg(i), qg(i));
   problem.x0 = [0; 0]; 
   problem.solver='fsolve'; 
   x=fsolve(problem); 
   delta(i)= x(1); 
   e(i) = x(2); 
  
end

end

