function [K, trCostEstimate,gammaEstimate] = LQRstep( zS, z0, alpha,networkS)
%LQRstep calculates the required feedback gain for the LQR step based on
%the desired steady-state zS and the previous steady-state z0. 
% The formulation solves an SDP per equation (18) but with zS known. 
% 
% Description of Outputs: 
% 1. K: the LQR gain, size(2*G, 3*G) 
% 2. trCostEstimate: is the estimate of the transient cost
% 
% Description of Inputs:
% 1. zS: the desired steady-state z=(x,a,u), size(2*N+7*G,1).
% 2. z0: the previous steady-state z0=(x0,a0,u0), size(2*N+7*G,1).
% 3. alpha

global Sbase N G L...
    deltaIdx omegaIdx eIdx mIdx...
    thetaIdx vIdx pgIdx qgIdx prefIdx fIdx...
 tau_vec xd_vec xq_vec xprime_vec d_vec m_vec Tch_vec freqR_vec
 



%% Obtaining the jacobians:
 [ gx,ga,gu ] = gFunctionJacobVectorized(z0);
[ hx, ha,hu ] = hFunctionJacobVectorized( z0 );
Asys=gx-ga*inv(ha)*hx;
Bsys=gu-ga*inv(ha)*hu;

%% costs:


delta0=z0(deltaIdx); 
omega0=z0(omegaIdx); 
e0=z0(eIdx); 
m0=z0(mIdx); 

x0=[delta0;omega0;e0;m0];



deltas=zS(deltaIdx); 
omegas=zS(omegaIdx); 
es=zS(eIdx); 
ms=zS(mIdx);
vs=zS(vIdx);
thetas=zS(thetaIdx); 
pgs=zS(pgIdx); 
qgs=zS(qgIdx); 
prefs=zS(prefIdx); 
fs=zS(fIdx); 
xs=[deltas;omegas;es;ms];






cvx_begin quiet
cvx_solver sdpt3
variables   gama2
variable P(4*G,4*G) symmetric
variable Y(2*G,4*G) 
expression Qinv(4*G,4*G) 
expression Rinv(2*G,2*G) 
minimize( gama2) ;
subject to:







[Qinv,Rinv]=QinvRinv(pgs,qgs,alpha, Qinv, Rinv,networkS);


-[-gama2, xs.'-x0.'; xs-x0, -P] ==semidefinite(4*G+1);

 -[Asys*P + P*Asys.' + Bsys*Y + Y.' * Bsys.', P , Y.';
    P, -Qinv, zeros(4*G, 2*G); 
    Y, zeros(2*G, 4*G), -Rinv] == semidefinite(10*G);

P==semidefinite(4*G);
   
cvx_end

K=-Rinv*Bsys.'*conj(inv(P));
trCostEstimate=gama2;

gammaEstimate=(xs-x0).'*inv(P)*(xs-x0);



end

