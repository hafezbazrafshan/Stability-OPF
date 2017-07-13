function [ gx,ga,gu ] = gFunctionJacobVectorized(z)
%GFUNCTIONJACOBVECTORIZED calculates the jacobian of g(x,a,u); 
% [gx,ga,gu]= gFunctionJacobVectorized(z, m_vec, d_vec, tau_vec, xd_vec,xprime_vec)
% calculates the jacobian of g(x,a,u). 
% 
% Description of Outputs: 
% 1. gx: the jacobian of g with respect to x, size(4*G,4*G); 
% 2. ga: the jacobian of g with respect to a, size(4*G, 2*N+2*G)
% 3. gu: the jacobian of g with respect to u, size(4*G, 2*G)
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
% 7. Tch_vec
% 8. freqR_vec
% See also gFunctionJacob

global   m_vec d_vec tau_vec xd_vec xprime_vec Tch_vec freqR_vec

global N G  gen_set ...
    deltaIdx omegaIdx eIdx mIdx...
    thetaIdx vIdx pgIdx qgIdx prefIdx fIdx 

M=diag(m_vec);
D=diag(d_vec);
T=diag(tau_vec);
Xprime=diag(xprime_vec); 
Xd=diag(xd_vec);
Tch=diag(Tch_vec); 
freqR=diag(freqR_vec); 


V=z(vIdx);
theta=z(thetaIdx);
Vg=V(gen_set);
V_g=diag(Vg);
thetag=theta(gen_set);
VgIdx=vIdx(gen_set);
thetagIdx=thetaIdx(gen_set);

delta=z(deltaIdx);


gz=zeros(4*G, 2*N+8*G); 



g1Idx=(1:G).'; 
g2Idx=(g1Idx(end)+1:g1Idx(end)+G).';
g3Idx=(g2Idx(end)+1:g2Idx(end)+G).';
g4Idx=(g3Idx(end)+1:g3Idx(end)+G).';

gz(g1Idx, omegaIdx) =eye(G);
gz( g2Idx, omegaIdx)=- inv(M)*D;
gz(g2Idx,mIdx)=inv(M); 
gz(g2Idx,pgIdx)=-inv(M);
gz(g3Idx,deltaIdx)= -inv(T)*inv(Xprime)*(Xd-Xprime)*V_g*diag(sin(delta-thetag)); 
gz(g3Idx,eIdx)=-inv(T)*inv(Xprime)*Xd;
gz(g3Idx, fIdx)=inv(T);
gz(g3Idx,VgIdx)=inv(T)*inv(Xprime)*(Xd-Xprime)*diag(cos(delta-thetag));
gz(g3Idx,thetagIdx) =inv(T)*inv(Xprime)*(Xd-Xprime)*V_g*diag(sin(delta-thetag)); 

gz(g4Idx,mIdx)=- inv(Tch);
gz(g4Idx,omegaIdx)=-inv(Tch)*inv(freqR);

gz(g4Idx,prefIdx)=inv(Tch);


gx=gz(:, [deltaIdx; omegaIdx; eIdx; mIdx]); 
ga=gz(:,[vIdx;thetaIdx;pgIdx;qgIdx]); 
gu=gz(:,[prefIdx;fIdx]); 




end

