function [ gx,ga,gu ] = gFunctionJacob(z,...
     m_vec, d_vec, tau_vec, xd_vec,...
    xprime_vec)
%GFUNCTIONJACOB calculates the jacobian of g(x,a,u); 
% [gx,ga,gu]= gFunctionJacob(z, m_vec, d_vec, tau_vec, xd_vec,xprime_vec)
% calculates the jacobian of g(x,a,u). 
% 
% Description of Outputs: 
% 1. gx: the jacobian of g with respect to x, size(3*G,3*G); 
% 2. ga: the jacobian of g with respect to a, size(3*G, 2*N+2*G)
% 3. gu: the jacobian of g with respect to u, size(3*G, 2*G)
% 
% Description of Inputs: 
% 1. z: vector of z=(x,a,u) combining states, algebraic and control
% variables
% 2. m_vec: vector of generator inertia constants in pu, size(G,1). 
% 3. d_vec: vector of damping coefficients, size(G,1)
% 4. tau_vec: vector of direct axis transient open-circuit time constant,
% size(G,1). 
% 5. xd_vec:  vector of direct axis synchronous reactance pu, size(G,1).
% 6. xprime_vec: direct axis transient reactance pu, size(G,1).
%
% See also gFunctionJacobVectorized


global N G  gen_set ...
    deltaIdx omegaIdx eIdx ...
    thetaIdx vIdx pgIdx qgIdx mIdx fIdx


V=z(vIdx);
theta=z(thetaIdx);
Vg=V(gen_set);
thetag=theta(gen_set);
VgIdx=vIdx(gen_set);
thetagIdx=thetaIdx(gen_set);

delta=z(deltaIdx);



gz=zeros(3*G, 2*N+7*G); 
gznCols=2*N+7*G;
gznRows=3*G;

g1Idx=1:G; 
g2Idx=g1Idx(end)+1:g1Idx(end)+G;
g3Idx=g2Idx(end)+1:g2Idx(end)+G;

gz(sub2ind( [gznRows gznCols], g1Idx, omegaIdx)) =1;
gz( sub2ind([gznRows gznCols], g2Idx, omegaIdx))=- d_vec./m_vec; 
gz(sub2ind([gznRows gznCols], g2Idx, mIdx))=1./m_vec; 
gz(sub2ind([gznRows gznCols], g2Idx, pgIdx))=-1./m_vec;
gz(sub2ind([gznRows gznCols], g3Idx, deltaIdx))= - ((xd_vec-xprime_vec)./xprime_vec).*Vg.*sin(delta-thetag).*(1./tau_vec);
gz(sub2ind([gznRows gznCols], g3Idx, eIdx))=-(1./tau_vec).*(xd_vec./xprime_vec);
gz(sub2ind([gznRows gznCols], g3Idx, fIdx))=1./tau_vec;
gz(sub2ind([gznRows gznCols], g3Idx, VgIdx))=(1./tau_vec).*( (xd_vec-xprime_vec)./xprime_vec).*cos(delta-thetag); 
gz(sub2ind([gznRows gznCols], g3Idx, thetagIdx)) =(1./tau_vec).*( (xd_vec-xprime_vec)./xprime_vec).*Vg.*sin(delta-thetag); 



gx=gz(:, [deltaIdx; omegaIdx; eIdx]); 
ga=gz(:,[vIdx;thetaIdx;pgIdx;qgIdx]); 
gu=gz(:,[mIdx;fIdx]); 

end