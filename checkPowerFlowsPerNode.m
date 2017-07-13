
function [checkpf, checkEqs,realGen_check, reactiveGen_check, ...
    realLoad_check,reactiveLoad_check]...
    = checkPowerFlowsPerNode(VS,thetaS,pgS,qgS, pdS,qdS)
% CHECKPOWERFLOWS Validates given power flow solution.
% [checkpf,checkEqs,realGen_check,...
%     reactiveGen_check, realLoad_check,...
%     reactiveLoad_check] = checkPowerFlows(VS, thetaS, pgS, qgS, pdS, qdS )  
%  validates given power flow solution.  This is a node by node
%  implementation. See checkPowerFlows for a vectorized implementation.
%  
% Description of outputs:
%  1. checkpf:  is a scalar binary which equals  1 when all power flow
%  equations are satisfied with absolute accuracy 1e-3. 
%  2. checkEqs: is a vector of size(2*N,1), expressing the difference with zero
% for any of the equations in CDC 2016 paper equations (2c)-(3b)
%  3. realGen_check: is a vector of size(G,1), expressing the difference with
%  zero for equations (2c)
%  4. reactiveGen_check: is a vector of size(G,1), expressing the difference
%  with zero for equations (2d)
% 5. realLoad_check: is a vector of size(L,1), expressing the difference with
% zero for equations (3a)
% 6. reactiveLoad_check: is a vector of size(L,1), expressing the difference
% with zero for equations (3d)
%
% Description of inputs:
%  1. VS: the steady-state voltage magnitude power flow solution
% 2. thetaS: the steady-state voltage phase power flow solution (in Radians)
% 3. pgS: the generator real power set points set points and the calculated pg for slack bus 
% 4. qgS: the generator reactive power inputs
% 5. pdS: the steady-state real power loads used to run the power flow
% 6. qdS: the steady-state reactive power loads used to run the power flow
%
% See also checkPowerFlows
%
% Required modifications:
% 1. Fix equation number references.

global N G L node_set gen_set...
    load_set Ymat Gmat Bmat



realGen_check=ones(G,1);
reactiveGen_check=ones(G,1);

% checking the first 2G equations 
% [manual equations (16a), (16b)]
for ii=1:length(gen_set)
idx=gen_set(ii); % network index
realGen_check(ii)=+pdS(idx)- pgS(ii) +...
    VS(idx).*(Gmat(idx,:)*(VS.*cos(thetaS(idx)-thetaS))+Bmat(idx,:)*(VS.*sin(thetaS(idx)-thetaS)));
reactiveGen_check(ii) =+qdS(idx)-qgS(ii) +...
  VS(idx).*(-Bmat(idx,:)*(VS.*cos(thetaS(idx)-thetaS))+Gmat(idx,:)*(VS.*sin(thetaS(idx)-thetaS)));
end


% checking the 2L equations for power flows
% [manual (16c), (16d)]
realLoad_check=ones(L,1);
reactiveLoad_check=ones(L,1);
for ii=1:length(load_set)
idx=load_set(ii); % network index
realLoad_check(ii)=pdS(idx)+...
    VS(idx).*(Gmat(idx,:)*(VS.*cos(thetaS(idx)-thetaS))+Bmat(idx,:)*(VS.*sin(thetaS(idx)-thetaS)));
reactiveLoad_check(ii) =+qdS(idx) +...
  VS(idx).*(-Bmat(idx,:)*(VS.*cos(thetaS(idx)-thetaS))+Gmat(idx,:)*(VS.*sin(thetaS(idx)-thetaS)));
end

checkEqs=[realGen_check; reactiveGen_check; realLoad_check;reactiveLoad_check];

if sum(abs(checkEqs))<1e-3
disp('Power flow equations satisfied'); 
checkpf=1;
else
disp('Power flow equations NOT satisfied'); 
checkpf=0;
end
end

