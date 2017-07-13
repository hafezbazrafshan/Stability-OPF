function [deltaVec, omegaVec, eVec, mVec,...
    thetaVec, vVec, pgVec, qgVec, ...
    prefVec,fVec, ...
    pdVec, qdVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec] = retrieve_output( t, ZNEW , NoiseVector)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global ControlMode
% system constants
global OMEGA_S Sbase N G L node_set gen_set load_set Ymat Gmat Bmat Cg...
    yff_vec yft_vec  ytf_vec ytt_vec

% indices
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx   zIdx 

% machine 
global  tau_vec xd_vec xq_vec xprime_vec d_vec m_vec Tch_vec freqR_vec ...
    

% dynamical simulations
global Tfinal Tpert fsample n_samples n_pertSamples Mass...
pertSet pPertValues qPertValues NoiseVarianceSet 


% initial conditions
global x0 omega0 delta0 e0 m0...
    a0 v0 theta0 pg0 qg0...
    u0 pref0 f0...
    pd0 qd0...
    vg0 thetag0...
    z0 network


% 0plus conditions?
global x0plus omega0plus delta0plus e0plus m0plus...
    a0plus v0plus theta0plus pg0plus qg0plus...
    u0plus pref0plus f0plus...
    pd0plus qd0plus...
    vg0plus thetag0plus...
    z0plus 

% next time-slot conditions
global xS omegaS deltaS eS mS...
    aS vS thetaS pgS qgS...
    uS prefS fS...
    pdS qdS...
    vgS thetagS...
    zS networkS
 
 % global 
global  KLQRstep






disp('Retrieving the output');
% retrieving time-varying quantities: 
deltaVec=ZNEW(deltaIdx,:); 
omegaVec=ZNEW(omegaIdx,:); 
eVec=ZNEW(eIdx,:); 
mVec=ZNEW(mIdx,:);
thetaVec=ZNEW(thetaIdx,:); 
vVec=ZNEW(vIdx,:); 
pgVec=ZNEW(pgIdx,:); 
qgVec=ZNEW(qgIdx,:); 



% retrieving load 
pdVec=zeros(N,n_samples); 
qdVec=zeros(N,n_samples);





prefVec=zeros(G,n_samples); 
fVec=zeros(G,n_samples);



   for tt=1:n_samples
 [prefVec(:,tt), fVec(:,tt)]=   control_law('LQR',deltaVec(:,tt),omegaVec(:,tt),eVec(:,tt),mVec(:,tt),...
                    vVec(:,tt), thetaVec(:,tt), pgVec(:,tt),qgVec(:,tt));
    
end



% computing xdot:
deltaDotVec=zeros(G,n_samples); 
omegaDotVec=zeros(G,n_samples); 
eDotVec=zeros(G,n_samples); 
mDotVec=zeros(G,n_samples); 

for tt=1:n_samples
[ deltaDotVec(:,tt), omegaDotVec(:,tt), eDotVec(:,tt) , mDotVec(:,tt)] = gFunctionVectorized( ...
    deltaVec(:,tt), omegaVec(:,tt), eVec(:,tt),mVec(:,tt),...
     vVec(gen_set,tt), thetaVec(gen_set,tt),pgVec(:,tt),prefVec(:,tt),fVec(:,tt));
end

for tt=1:n_samples
    [pdVec(:,tt),qdVec(:,tt)]=loadPert('Transient', t(tt),pd0,qd0,pertSet,pPertValues, qPertValues,Tpert,Tfinal,NoiseVector) ;   
end







end

