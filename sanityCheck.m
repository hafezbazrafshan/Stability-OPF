function [ sanitycheck1,sanitycheck2,sanitycheck3 , success] = sanityCheck(...
    deltaVec, omegaVec, eVec, mVec, ...
    thetaVec, vVec, pgVec, qgVec,...
    prefVec, fVec, ...
    ploadVec,qloadVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

global  N G L  gen_set load_set ...
    deltaIdx omegaIdx eIdx mIdx  ...
 xq_vec xprime_vec ...
 fsample n_samples...


disp('Performing sanity checks on dynamical simulation'); 
% sanity check 1: checking power flow equatios at each time-step:
% check power flow equations
checkpfVec=zeros(1,n_samples); 
checkEqsVec=zeros(2*N,n_samples); 
realGen_checkVec=zeros(G,n_samples); 
reactiveGen_checkVec=zeros(G,n_samples); 
realLoad_checkVec=zeros(L,n_samples); 
reactiveLoad_checkVec=zeros(L,n_samples); 
for tt=1:n_samples
    
pload=ploadVec(:,tt); 
qload=qloadVec(:,tt);
[checkpfVec(:,tt), checkEqsVec(:,tt),realGen_checkVec(:,tt), reactiveGen_checkVec(:,tt), ...
    realLoad_checkVec(:,tt),reactiveLoad_checkVec(:,tt)]=...
   checkPowerFlows(vVec(:,tt),thetaVec(:,tt),pgVec(:,tt),qgVec(:,tt), pload,qload);


end

if all(checkpfVec==1)
    disp(' Sanity check 1: in transient simulations, power flow equations were all satisfied at every time step.'); 
    sanitycheck1=1;
else 
    disp(' Sanity check 1: FAILED, power flow equaitons were not satisfied at every time step.'); 
        sanitycheck1=0;
end

% Sanity check 2:  checking the entire algebraic equations (kind of
% redundant for power flow equations that are already checked but OK)
h1Idx=1:G;
h2Idx=G+1:2*G;
h3Idx=2*G+1:3*G;
h4Idx=3*G+1:4*G;
h5Idx=4*G+1:4*G+L;
h6Idx=4*G+L+1:4*G+2*L;
d=zeros(h6Idx(end),1);

h1Vec=zeros(G,n_samples); 
h2Vec=zeros(G,n_samples); 
h3Vec=zeros(G,n_samples); 
h4Vec=zeros(G,n_samples); 
h5Vec=zeros(L,n_samples); 
h6Vec=zeros(L,n_samples); 
sanity2Vec=zeros(2*N+2*G,n_samples); 

for tt=1:n_samples
[h1,h2,h3,h4,h5,h6] = hFunctionVectorized(deltaVec(:,tt),eVec(:,tt),vVec(:,tt),thetaVec(:,tt),...
    pgVec(:,tt),qgVec(:,tt));
pload=ploadVec(:,tt); 
qload=qloadVec(:,tt);
ploadg=pload(gen_set);
qloadg=qload(gen_set);
ploadl=pload(load_set);
qloadl=qload(load_set);
d(h3Idx)=-ploadg;
d(h4Idx)=-qloadg;
d(h5Idx)=-ploadl;
d(h6Idx)=-qloadl;   
hz=[h1;h2;h3;h4;h5;h6];
sanity2Vec(:,tt)=hz-d;

end

if max(max(abs(sanity2Vec)))<1e-4
    disp('Sanity check 2: in transient simulations, all algebraic equations were satisfied.'); 
    sanitycheck2=1;
else 
    disp('Sanity check 2: FAILED, some algebraic equations were not satisfied'); 
    sanitycheck2=0;
end

% sanity check 3: comparison of calculated state derivative and numerical
% derivative:

err1=abs(deltaDotVec(:,1:end-1) - numericalDerivative(deltaVec,1/fsample)); 
err2=abs(omegaDotVec(:,1:end-1)- numericalDerivative(omegaVec,1/fsample));
err3=abs(eDotVec(:,1:end-1) - numericalDerivative(eVec,1/fsample)); 
 err4=abs(mDotVec(:,1:end-1)- numericalDerivative(mVec, 1/fsample)); 

err1=err1(:,2:end); 
err2=err2(:,2:end); 
err3=err3(:,2:end);
err4=err4(:,2:end);

if sum(max([err1;err2;err3;err4])>1e-1)<0.2*n_samples
    disp('Sanity check 3: numerical and computed time derivatives for states match'); 
    sanitycheck3=1;
else 
      disp('Sanity check 3: FAILED, numerical and computed time derivatives for states DO NOT match'); 
      sanitycheck3=0;
end


if sum([sanitycheck1,sanitycheck2,sanitycheck3])==3
    success=1; 
else 
    success=0;
end

end

