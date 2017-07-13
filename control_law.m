function [pref,f]=control_law(ControlMode,delta,omega,e,m,...
                    v, theta, pg,qg)
      
                
           
% system constant variables
global network Sbase N G L node_set gen_set load_set Ymat Gmat Bmat Cg...
    
% index variables 
global  deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx xIdx  aIdx uIdx OMEGA_S...
    
    

% machine constants
global tau_vec xd_vec xq_vec xprime_vec d_vec m_vec Tch_vec freqR_vec ...
    

% dynamical simulation parameters
global  Tfinal Tpert fsample n_samples n_pertSamples Mass...
    

% disturbance parameters
global pertSet pPertValues qPertValues


%initial conditions
global x0 a0 u0 delta0 omega0 e0 m0 v0 theta0 pg0 qg0 pref0 f0...
   vg0 thetag0  pd0 qd0 % some other useful initial conditions

%initial conditions at zero plus
global x0plus a0plus u0plus delta0plus omega0plus e0plus m0plus v0plus theta0plus pg0plus qg0plus pref0plus f0plus...
    vg0plus thetag0plus pd0plus qd0plus...
    deltaDot0plus omegaDot0plus eDot0plus mDot0plus...% % some other useful initial conditions at zero plus
xDot0plus aDot0plus

% new equilibrium 
global xS aS uS deltaS omegaS eS mS vS thetaS pgS qgS prefS fS...
    vgS thetagS pdS qdS...
    deltaDotS omegaDotS eDotS mDotS...
    networkS

% LQR control
global KLQRstep





                
     switch ControlMode
         
         
         case 'OpenLoop'
             
             
             pref=prefS;
             f=fS;
             
         case 'LQR'
             
             u=KLQRstep*[delta- deltaS; omega- omegaS; e-eS;m-mS]+[prefS;fS];
             
             
             
             pref=u(1:G); 
            f=u(G+1:end);


             end

             
             

     end
             
             
     
     
