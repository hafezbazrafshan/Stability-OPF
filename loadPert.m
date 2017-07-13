function [pload,qload, noise_vector]=loadPert(timeMod, t,pload0,qload0,pertSet, pPertValues, qPertValues,Tpert,Tfinal,noise_vector)
%LOADPERT produces a new load starting at Tpert
%   [pload,qload]=loadPert(timeMod, t,pload0,qload0,nodepert_set, pertvalue_set, Tpert)
%   produces a new  real and reactive power load at Tpert 
% 
% Description of Outputs: 
% 1. pload: the new real power load of size(N,1);
% 2. qload: the new reactive power load of size(N,1);
% 
% Description of Inputs:
% 1. timeMod: a switch option for 'Steady-State' or 'Transient'. 
% 2. t: the time instant (a scalar). If timeMod is set to 'Steady-State' is not used.
% 3. pload0: the initial real power load of size(N,1); 
% 4. qload0: the intial reactive power load of size(N,1); 
% 5. nodepert_set: set of nodes in which we would like the load to be
% perturbed. 
% 6. pertvalue_set: set of pertubration values corresponding to
% nodepert_set
% 7. Tpert: time of load perturbation. 
%
% Required:
% 1. Extend to Tpert vector

global NoiseVarianceSet fsample

switch timeMod
    
    case 'Steady-State'
 N=length(pload0);
pstep=zeros(N,1);
pstep(pertSet)=pPertValues;
qstep=zeros(N,1);
qstep(pertSet)=qPertValues;

pload= pload0 +pstep;
qload=qload0+qstep;
        
    case 'Transient'
% N=length(pload0);
% step=zeros(N,1);
% step(pertSet)=pPertValues;
% Tramp=Tfinal/2;
% 
noise=zeros(size(noise_vector(:,1)));
if t>0
n1=floor(t*fsample);
n2=ceil(t*fsample);
alpha=n2-t*fsample;
noise=alpha*noise_vector(:,n1+1)+(1-alpha)*noise_vector(:,n2+1);
end
% 
% if t<=Tramp
% pload=pload0+step.*((t-Tpert)./(Tramp-Tpert)).*(t>Tpert)'+noise.*(t>Tpert)';
% else 
%     pload=pload0+step.*(t>Tpert).'+noise.*(t>Tpert)';
% %     qload=qload0+step.*(t>Tpert).'+noise.*(t>Tpert)';
% end
% 
% qload=qload0;


 N=length(pload0);
pstep=zeros(N,1);
pstep(pertSet)=pPertValues;
qstep=zeros(N,1);
qstep(pertSet)=qPertValues;

pload= pload0 +pstep+noise;
qload=qload0+qstep;
        
end



end

