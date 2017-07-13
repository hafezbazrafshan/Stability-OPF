function [ fz] = runDynamics(t,z,zDot,NoiseVector)
%odeFunction15i puts together the equations of the descriptor system based
%on the requirements of the MATLAB ode15i solver. 
%  [ fz] = odeFunction15i(t,znew,zDot, K,z0, zS, pload0, qload0, Tpert,
%  nodepert_set,pertvalue_set) puts the DAE equations of $\dot{x}=g(x,a,u)$ and
%  $d=h(x,a,u)$ of CDC 2016 equations (4) and (5) in the following form:
% $0 = fz= f(t,z, \dot{z}) = [E\dot{z} -gz; hz-d]$
% 
% Description of Outputs: 
% 1. fz: a vector of size(2N+2*G,1) at time t fz= f(t,z, zDot) = [MASS*zDot -gz; hz-d]$
% 
% Description of Inputs: 
% 1. t: the time instant for the functions to be evaluated
% 2. z: the variable z, size(2*N+5*G,1).
% 3. zDot: the derivative of variable z, size(2*N+5*G,1).
% 4. K: the linear feedback gain that relates the controls to states,
% size(2*G, 2*G).
% 5. z0: the initial conditions to start the simulations, where variable z
% starts from
% 6. zS: the desired condition to end the simulations, where we hope that
% the controller u=Kx will take variable z to. 
% 7. pload0: the initial real load condition, size(G,1). 
% 8. qload0; the initial reactive load condition, size(G,1).
% 9. Tpert: the perturbation time instant
% 10. nodepert_set: the set of nodes where load pertubration occurs
% 12. pertvalue_set: the set of load perturbation values corresponding to
% nodepert_set
% 
% Required:
%

% system constants [these do not change]
global ControlMode

global OMEGA_S Sbase N G L node_set gen_set load_set Ymat Gmat Bmat Cg...
    yff_vec yft_vec  ytf_vec ytt_vec

%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx   zIdx 

% machine [these do not change]
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


% 0plus conditions
global x0plus omega0plus delta0plus e0plus m0plus...
    a0plus v0plus theta0plus pg0plus qg0plus...
    u0plus pref0plus f0plus...
    pd0plus qd0plus...
    vg0plus thetag0plus...
    z0plus...
    deltaDot0plus omegaDot0plus eDot0plus mDot0plus

% next time-slot conditions
global xS omegaS deltaS eS mS...
    aS vS thetaS pgS qgS...
    uS prefS fS...
    pdS qdS...
    vgS thetagS...
    zS networkS
 
 % global 
global  KLQRstep





delta=z(deltaIdx);
omega=z(omegaIdx);
e=z(eIdx);
m=z(mIdx); 
theta=z(thetaIdx);
v=z(vIdx);
pg=z(pgIdx);
qg=z(qgIdx);
thetag=theta(gen_set);
vg=v(gen_set);




%

if strcmp(ControlMode,'LQR')
    
[pref,f]=control_law('LQR',delta,omega,e,m,...
                    v, theta, pg,qg);
end


        















h1Idx=(1:G).';
h2Idx=(G+1:2*G).';
h3Idx=(2*G+1:3*G).';
h4Idx=(3*G+1:4*G).';
h5Idx=(4*G+1:4*G+L).';
h6Idx=(4*G+L+1:4*G+2*L).';
d=zeros(h6Idx(end),1);

[pload,qload]=loadPert('Transient',t,pd0,qd0,pertSet, pPertValues,qPertValues, Tpert,Tfinal,NoiseVector);
ploadg=pload(gen_set);
qloadg=qload(gen_set);
ploadl=pload(load_set);
qloadl=qload(load_set);
d(h3Idx)=-ploadg;
d(h4Idx)=-qloadg;
d(h5Idx)=-ploadl;
d(h6Idx)=-qloadl;


[ deltaDot, omegaDot, eDot, mDot] =gFunctionVectorized( ...
    delta, omega, e,m,...
    vg,thetag, pg,pref,f);


gz=[deltaDot; omegaDot; eDot;mDot];


[h1,h2,h3,h4,h5,h6] = hFunctionVectorized(delta,e,v,theta,...
    pg,qg);

hz=[h1;h2;h3;h4;h5;h6];
    fz=[[zDot(deltaIdx); zDot(omegaIdx); zDot(eIdx); zDot(mIdx)] - gz;hz-d];




end

