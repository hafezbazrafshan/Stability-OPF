function [vgs,pgsNonSlack, thetaSlack, K,ssCost, trCostEstimate ] = augmentedOPF( z0,networkS,...
    deltaploadg,deltaqloadg,deltaploadl,deltaqloadl, alpha, Tlqr)
%AUGMENTEDOPF  implements augmented opf per equation (18) CDC 2016. 
%    [vs,thetas, pgs,qgs, K ] = augmentedOPF( z0,...
 %   deltaploadg,deltaqloadg,deltaploadl,deltaqloadl) implements the
 %   augmetned OPF based on linear approximation of a known equilibrium z0
 % for the power systems described by nonlinear equations g(x,a,u) and
 % h(x,a,u). 
 %
 % Description of Outputs: 
 % 1. vs: the calculated optimal steady-state voltage magnitude, size(N,1).
 % 2. thetas: the calculated optimal steady-state voltage angle in radians,
 % size(N,1). 
 % 3. pgs: the calculated optimal steady-state real power injection
 % (setpoints) in pu Watts, size(G,1). 
 % 4. qgs: the calculated optimal steady-state reactive power injection 
% in pu Vars, size(G,1). 
% 5. K: the calculated optimal linear feedback gain, size(2*G,3*G)
% 6. ssCost: is the calculated steady-state cost of real power generation
% 7. trCostEstimate: is the gama---an estimate of the transient cost
% 
% Description of Inputs: 
% 1. z0: the  equilibrium point used for linearization
% 2. deltaploadg: the difference of the new desired  real power to the initial load level
% for generator nodes, size(G,1).
% 3. deltaqloadg: the difference of the new desired  reactive power to the initial load level
% for generator nodes, size(G,1).
% 4. deltaploadl: the difference of the new desired  real power to the initial load level
% for load nodes, size(L,1).
% 5. deltaqloadl: the difference of the new desired  reactive power to the initial load level
% for generator nodes, size(L,1)
% 6. alpha: alpha
% 7. Tlqr: is the Tlqr factor for transient control
% See also approxOPF, LQRstep
%
% Required:
% 

global Sbase N G L...
    deltaIdx omegaIdx eIdx mIdx...
    thetaIdx vIdx pgIdx qgIdx prefIdx fIdx...
    
global gen_set


 
    
%% Obtaining the jacobians:
 [ gx,ga,gu ] = gFunctionJacobVectorized(z0);
[ hx, ha,hu ] = hFunctionJacobVectorized( z0);
Asys=gx-ga*inv(ha)*hx;
Bsys=gu-ga*inv(ha)*hu;






%% costs:
c2k=networkS.gencost(:,5).*Sbase.^2;
c1k=networkS.gencost(:,6).*Sbase; 
c0k=networkS.gencost(:,7);

delta0=z0(deltaIdx); 
omega0=z0(omegaIdx); 
e0=z0(eIdx); 
m0=z0(mIdx);

x0=[delta0;omega0;e0;m0];


h1Idx=1:G;
h2Idx=G+1:2*G;
h3Idx=2*G+1:3*G;
h4Idx=3*G+1:4*G;
h5Idx=4*G+1:4*G+L;
h6Idx=4*G+L+1:4*G+2*L;
deltaD=zeros(h6Idx(end),1);


deltaD(h3Idx)=-deltaploadg;
deltaD(h4Idx)=-deltaqloadg;
deltaD(h5Idx)=-deltaploadl;
deltaD(h6Idx)=-deltaqloadl;



cvx_begin quiet
cvx_solver sdpt3
variables zs(2*N+8*G)   gama2
variable P(4*G,4*G) symmetric
variable Y(2*G,4*G) 
expression Qinv(4*G,4*G) 
expression Rinv(2*G,2*G) 

deltas=zs(deltaIdx); 
omegas=zs(omegaIdx); 
es=zs(eIdx); 
ms=zs(mIdx);
vs=zs(vIdx);
thetas=zs(thetaIdx); 
pgs=zs(pgIdx); 
qgs=zs(qgIdx); 
prefs=zs(prefIdx); 
fs=zs(fIdx); 

xs=[deltas;omegas;es;ms];

minimize( c2k.'*square(pgs) + c1k.'*pgs+c0k.'*ones(G,1)+ (Tlqr)*gama2) ;
subject to:




omegas==omega0;

-pi<=thetas<=pi





% Qinv(sub2ind([4*G,4*G], omegaIdx, omegaIdx)) =0.0001*(-alpha*pgs + network.gen(:,9)./Sbase);
% Qinv(sub2ind([4*G,4*G], deltaIdx, deltaIdx)) =0.0001*( -alpha*pgs + network.gen(:,9)./Sbase);
% Qinv(sub2ind([4*G, 4*G], eIdx,eIdx)) =0.01* (-alpha*qgs+ network.gen(:,4)./Sbase);
% Qinv(sub2ind([4*G,4*G],mIdx,mIdx))=0.005* (-alpha*pgs+ network.gen(:,9)./Sbase);
% 
% 
% Rinv(sub2ind([2*G, 2*G], 1:G,1:G))  =0.005*(-alpha*pgs + network.gen(:,9)./Sbase);
% Rinv(sub2ind([2*G, 2*G], G+1:2*G, G+1:2*G)) =0.01*(-alpha*qgs + network.gen(:,4)./Sbase);
[Qinv,Rinv]=QinvRinv(pgs,qgs,alpha, Qinv, Rinv,networkS);


zeros(4*G,1) == [gx, ga, gu]*( zs-z0); 
deltaD==[hx,ha,hu]*(zs-z0);

networkS.bus(:,13)<= vs<=networkS.bus(:,12);
networkS.gen(:,5)./Sbase<=qgs<=networkS.gen(:,4)./Sbase;
networkS.gen(:,10)./Sbase <= pgs <= networkS.gen(:,9)./Sbase; 

slackIdx=find(networkS.bus(:,2)==3); % finds the slack bus
thetas(slackIdx)==0;

-[-gama2, xs.'-x0.'; xs-x0, -P] ==semidefinite(4*G+1);

 -[Asys*P + P*Asys.' + Bsys*Y + Y.' * Bsys.', P , Y.';
    P, -Qinv, zeros(4*G, 2*G); 
    Y, zeros(2*G, 4*G), -Rinv] == semidefinite(10*G);
% 
P==semidefinite(4*G);
   
cvx_end

K=-Rinv*Bsys.'*conj(inv(P));

ssCost=c2k.'*(pgs.^2) + c1k.'*pgs+c0k.'*ones(G,1);
trCostEstimate=gama2;

vgs=vs(gen_set);
pgsNonSlack=pgs(networkS.bus(gen_set,2)==2);
thetaSlack=thetas(slackIdx);
end

