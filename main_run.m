function main_run(casefile,stcontrol,lfcontrol,alpha)
%% Guaranteeing that matpower is on the filepath:
% The code uses MATPOWER, and it is important to make sure matpower is
% enabled.
clc;
close all;
currentdirectory=pwd;
cd('..'); 
try
cd('matpower6.0/'); 
matpower_directory=pwd;
cd(currentdirectory); 
addpath(matpower_directory); 
disp('MATPOWER was sucessfully added to the path'); 
catch 
disp('ERROR: unable to find MATPOWER')
end


% The following is being issue since CVX uses nargin.
%  Matlab has depcrecated nargin and is using nargchk.
% The warnings are annoying so I turned them off.
warning('off','MATLAB:nargchk:deprecated')

global ControlMode
ControlMode=lfcontrol;

global SteadyStateMode
SteadyStateMode=stcontrol;


%% Defining some global variables
% Some of these variables refer to the network. For example, Sbase, Ymat
% (the bus admittance matrix), etc.  Some others are specific indices
% within a vector.  Some are LQR parameters.  
% These variables are declared global to avoid extra function
% arguments. 

% system constants [these do not change]
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


if strcmp(ControlMode,'AGC')
global y0 y0plus yS yDot0plus yDotS yIdx...
    ParticipationFactors NumberOfAreas AreaSet TieLineFromSet TieLineToSet...
    ACE0plus PscheduledS GensPerArea BusesPerArea...
    K_I K_ACE K_PG K_Pflow K_sumPG

end


matpower_options=mpoption('out.all',0); % this suppresses MATPOWER print output

%% 1.  Importing case39 with machine constants embedded.
% The purpose of this step is to have a so-called base network.
% network=runpf('casefiles/case39wmac_con', matpower_options);
casestr=['casefiles/',casefile'];
network=runpf(casestr, matpower_options);

% Derive bus admittance matrix and relevant network information


[ N,G,L,Ymat, Gmat, Bmat,...
    node_set, gen_set, load_set,...
   Cg, yff_vec, yft_vec, ytf_vec, ytt_vec] = networkParams( network );

Sbase=network.baseMVA;
OMEGA_S=2*pi*60; % synchronous frequency


%% Areas:
if strcmp(ControlMode,'AGC')
AreaSet=unique(network.bus(:,7)); 
NumberOfAreas=length(AreaSet); 
BusesPerArea=cell(NumberOfAreas,1); 
GensPerArea=cell(NumberOfAreas,1);
TieLineFromSet=cell(NumberOfAreas,1);
TieLineToSet=cell(NumberOfAreas,1);
for ii=1:NumberOfAreas
    BusesPerArea{ii,1}=find(network.bus(:,7)==AreaSet(ii)); 
    GensPerArea{ii,1}=find(network.bus(gen_set,7)==AreaSet(ii));
    TieLineFromSet{ii,1}= find(and(network.bus(network.branch(:,1),7)==AreaSet(ii),network.bus(network.branch(:,2),7)~=AreaSet(ii)));
    TieLineToSet{ii,1}=find(and(network.bus(network.branch(:,1),7)~=AreaSet(ii),network.bus(network.branch(:,2),7)==AreaSet(ii)));
end
end


%% 2.  populating the initialized steady-state variables from the MATPOWER case file:
v0=network.bus(:,8); % eighth column of bus matrix is voltage magnitude solution
theta0=degrees2radians(network.bus(:,9)); % nineth column of bus matrix is voltage phase solution
pg0=network.gen(:,2)./Sbase; % second column of gen matrix is real power set points (or solution for slack bus)
qg0=network.gen(:,3)./Sbase; % third column of gen matrix is reactive power solutions
pd0=network.bus(:,3)./Sbase; % third column of bus matrix is real power demand
qd0=network.bus(:,4)./Sbase; % fourth column of bus matrix is reactive power demand

% Verifying the initial power flow solution:
 [checkpf, checkEqs,realGen_check, reactiveGen_check, ...
    realLoad_check,reactiveLoad_check]=...
   checkPowerFlows(v0,theta0,pg0,qg0, pd0,qd0);
if checkpf==1
    disp('Initial power flow solution was correct'); 
else 
    disp('Initial power flow solution was incorrect'); 
    disp('Check case file'); 
    pause;
end

%% 3.  Adding transient parameters to the network:
% The original MATPOWER case file does not include transient parameters.
% The imported case file has it embedded. 
% retrieve machine constants:
Sbase2=network.mac_con(:,3); 
tau_vec=network.mac_con(:,9);
xd_vec=network.mac_con(:,6).*Sbase./Sbase2;
xq_vec=network.mac_con(:,11).*Sbase./Sbase2;
xprime_vec=network.mac_con(:,7).*Sbase./Sbase2;
d_vec=network.mac_con(:,17).*Sbase./Sbase2;
m_vec=network.mac_con(:,16)/(pi*60).*Sbase2./Sbase;


Tch_vec=0.5*ones(G,1); 
freqR_vec=0.02*ones(G,1)./(2*pi); 


%% 4. Obtain generator internal angles and electromotive force from the power flow solution
vg0=v0(gen_set);
thetag0=theta0(gen_set);

[ delta0, e0]=obtainGenStates(vg0, thetag0, pg0, qg0 );
omega0=repmat(OMEGA_S,G,1); % creating a vector of OMEGA_S of size(G,1), for all generator nodes.

%% 5.  Obtaining generator steady-state controls from the power flow solutions and steady-state of states
[m0,f0]=obtainGenControls(delta0,omega0,e0,vg0,thetag0,pg0,qg0, OMEGA_S);

pref0=m0;
if strcmp(ControlMode,'AGC')
y0=zeros(G,1); 
end

%% 6.  Defining the indices of vector z for dynamical simulation:

% *****states*****x
% delta size(G,1)
% omega size(G,1)
% e size(G,1)
deltaIdx=(1:G).';
omegaIdx=(deltaIdx(end)+1:deltaIdx(end)+G).';
eIdx=(omegaIdx(end)+1:omegaIdx(end)+G).';
mIdx=(eIdx(end)+1:eIdx(end)+G).';



%*****algebraic variables*****a
% theta size(N,1)
% v size(N,1)
% pg size(G,1)
% qg size(G,1)
vIdx=(mIdx(end)+1:mIdx(end)+N).';
thetaIdx=(vIdx(end)+1:vIdx(end)+N).';
pgIdx=(thetaIdx(end)+1:thetaIdx(end)+G).';
qgIdx=(pgIdx(end)+1:pgIdx(end)+G).';

if strcmp(ControlMode,'AGC')
yIdx=(qgIdx(end)+1:qgIdx(end)+G).';
end

%*****control variables*******u
%pref and f
prefIdx=(qgIdx(end)+1:qgIdx(end)+G).';
fIdx=(prefIdx(end)+1:prefIdx(end)+G).';



zIdx=(1:fIdx(end)).'; % should equal 2*N+8*G
z0=zeros(length(1:fIdx(end)),1); 

z0(deltaIdx)=delta0;
z0(omegaIdx)=omega0;
z0(eIdx)=e0;
z0(mIdx)=m0;

z0(thetaIdx)=theta0;
z0(vIdx)=v0;
z0(pgIdx)=pg0;
z0(qgIdx)=qg0;


z0(prefIdx)=pref0;
z0(fIdx)=f0;



%% 7.  Introducing new load for the next OPF time-slot:
disp('Perturbing the load'); 

pRatio=0.06; 
qRatio=0.06;
pertSet=find(or(network.bus(:,3)>0, network.bus(:,4)>0));
pPertValues=pRatio*pd0(pertSet);
qPertValues=qRatio*qd0(pertSet); 
NoiseVarianceSet=0.*network.bus(:,3)/Sbase;
for jj=1:length(pertSet)
   messageSTR=['Modified real and reactive power consumption of bus No. ', num2str(pertSet(jj)), ' respectively by ',...
       num2str(pPertValues(jj)), ' puWatts.', ' and ',  num2str(qPertValues(jj)), ' puVars.']; 
   disp(messageSTR);
    pause(0.1); 
end


% new steady-state conditions
[pdS,qdS]=loadPert('Steady-State',[],pd0,qd0,pertSet,pPertValues,qPertValues,[],[],[]);

networkS=network;
networkS.gen(:,9)=1.1*network.gen(:,9);
networkS.bus(:,3)= pdS.*Sbase; 
networkS.bus(:,4)=qdS.*Sbase; 

 %% 8. Setting LQR parameters
% alpha=0.8;
Tlqr=1;


%% 9.  Solving the augmentedOPF for the next time-slot


% setting matpower options need in subsequent load-flow
mpopt = mpoption('model', 'AC', 'pf.alg', 'NR', 'verbose', 0, 'out.all',0); 

if strcmp(SteadyStateMode,'LQR-OPF')
    % to obtain Vg, pgS, qgS
dpgS=pdS(gen_set)-pd0(gen_set);
dqdgS=qdS(gen_set)-qd0(gen_set);
dpdlS=pdS(load_set)-pd0(load_set);
dqdlS=qdS(load_set)-qd0(load_set);
[vgS,pgNonSlack,thetaSlack, ~, ~, ~ ] = augmentedOPF( z0,networkS,...
    dpgS,dqdgS,dpdlS,dqdlS,alpha, Tlqr); % we understand the slack bus angle is fixed at zero
networkS.gen(:,6)=vgS;
networkS.gen(networkS.bus(gen_set,2)==2,2)=pgNonSlack.*Sbase;
networkS.bus(networkS.bus(:,2)==3,9)=thetaSlack;
% run matpower power flow:
[networkS,successflag]=runpf(networkS,mpopt);
elseif strcmp(SteadyStateMode,'OPF')
  networkS=  runopf(networkS,mpopt);
  [networkS,successflag]=runpf(networkS,mpopt); % just reassuring
end    

%% 10. Obtain a true equilibrium for the next time-slot
vS= networkS.bus(:,8);
vgS=vS(gen_set);
thetaS= degrees2radians(networkS.bus(:,9));
thetagS=thetaS(gen_set); 
pgS=networkS.gen(:,2)./Sbase; 
qgS=networkS.gen(:,3)./Sbase;

c2k=networkS.gencost(:,5).*Sbase.^2;
c1k=networkS.gencost(:,6).*Sbase; 
c0k=networkS.gencost(:,7);

ssCost=c2k.'*square(pgS) + c1k.'*pgS+c0k.'*ones(G,1);
    


% check new power flow
 [checkpf2, checkEqs2,realGen_check2, reactiveGen_check2, ...
    realLoad_check2,reactiveLoad_check2]=...
   checkPowerFlows(vS,thetaS,pgS,qgS, pdS,qdS);

if checkpf2==1
    disp('The power flow solution for the second time slot is correct'); 
else 
    disp('The power flow solution for the second time slot in incorrect'); 
    disp('Check the new network conditions and MATPOWER runpf successflag'); 
    pause;
end

[ deltaS,eS]=obtainGenStates(vgS, thetagS, pgS, qgS );
omegaS=repmat(OMEGA_S,G,1);

% Obtaining generator steady-state controls
[mS,fS]=obtainGenControls( deltaS,omegaS,eS,vgS,thetagS,pgS,qgS, OMEGA_S);
prefS=mS;
zS=[deltaS; omegaS; eS;mS;  vS; thetaS; pgS; qgS; prefS; fS]; 




%% Run new LQR with setpoints determined to reach the true equilibrium
[KLQRstep,trCost2,gammaEstimate] = LQRstep( zS, z0,alpha,networkS);

%% 12. Define the MASS matrix (The E matrix in $E\dot{x}$ descriptor systems)
if strcmp(ControlMode,'LQR')
    Mass=zeros(4*G+2*N+2*G,4*G+2*N+2*G);
Mass(sub2ind(size(Mass), [deltaIdx,omegaIdx,eIdx,mIdx],[deltaIdx,omegaIdx,eIdx,mIdx]))=1;
elseif strcmp(ControlMode,'AGC')
Mass=zeros(4*G+2*N+2*G+G,4*G+2*N+2*G+G);
Mass(sub2ind(size(Mass), [deltaIdx,omegaIdx,eIdx,mIdx],[deltaIdx,omegaIdx,eIdx,mIdx]))=1;
Mass(sub2ind(size(Mass), [yIdx],[yIdx]))=1;
end
%% 13. Set dynamical simulation parameters:
options = odeset('Mass',Mass,'MassSingular', 'yes','MStateDependence','none', ...
    'RelTol',1e-7,'AbsTol',1e-6,'Stats','off');

Tfinal=5*60;
Tpert=0; % NEEDS TO BE CLOSE TO t=0 for increased accuracy 
fsample = 10;
n_samples = Tfinal * fsample+1;
n_pertSamples = max(Tpert,0) * fsample+1;
t = 0:1/fsample:Tfinal;
NoiseVector=repmat(NoiseVarianceSet,1,n_samples).*randn(length(NoiseVarianceSet),n_samples);





%% Simulation intial conditions at zero plus
display('Configuring 0plus intial conditions due to disturbance');
[pd0plus,qd0plus]=loadPert('Transient', 1,pd0,qd0,pertSet, pPertValues, qPertValues,Tpert,Tfinal,NoiseVector)


delta0plus=delta0; 
omega0plus=omega0; 
e0plus=e0; 
m0plus=m0;

% solve for v0plus, theta0plus, pg0plus, qg0plus
options2=optimset('MaxFunEvals',100000, 'MaxIter',100000,'disp','iter');
fun=@(x)algebraic( x,delta0plus,e0plus, pd0plus, qd0plus);
[x,fval] = fsolve(fun, [v0;theta0],options2);
algIdxV=1:N;
algIdxTheta=N+1:2*N;
v0plus=x(algIdxV); 
theta0plus=x(algIdxTheta); 
Xprime=diag(xprime_vec); 
Xq=diag(xq_vec);
vg0plus=v0plus(gen_set);
thetag0plus=theta0plus(gen_set); 
E0plus=diag(e0plus);
Vg0plus=diag(vg0plus);
pg0plus=inv(Xprime)*E0plus*Vg0plus*sin(delta0plus-thetag0plus)+...
   (1/2)*inv(Xq)*inv(Xprime)*(Xprime-Xq)*Vg0plus*Vg0plus*sin(2*(delta0plus-thetag0plus));
qg0plus=inv(Xprime)*E0plus*Vg0plus*cos(delta0plus-thetag0plus)-...
    (1/2)*inv(Xq)*inv(Xprime)*(Xprime+Xq)*Vg0plus*vg0plus+ ...
    (1/2)*inv(Xq)*inv(Xprime)*(Xprime-Xq)*Vg0plus*Vg0plus*cos(2*(delta0plus-thetag0plus));

%% AGC:

if strcmp(ControlMode,'AGC')
     ParticipationFactors=cell(NumberOfAreas,1);
    y0plus=zeros(G,1);
     PscheduledS=zeros(NumberOfAreas,1);
     K_I=1;
     K_ACE=	1;
     K_PG=1;
     K_sumPG=1;
     K_Pflow=1;
     for ii=1:NumberOfAreas
          [pFromS,~]=determineLineFlows(TieLineFromSet{ii,1},vS,thetaS);
               [~,pToS]=determineLineFlows(TieLineToSet{ii,1},vS,thetaS); 
               PscheduledS(ii,1)=sum(pFromS)-sum(pToS);
ParticipationFactors{ii,1}=pgS(GensPerArea{ii,1})./sum(pgS(GensPerArea{ii,1}));

     end

[ yDot0plus, ACE0plus, Pmeasured0plus, OmegaMeasured0plus] = agcParams(omega0plus,y0plus, v0plus,theta0plus, pg0plus);
end


%%
if strcmp(ControlMode,'AGC')
[pref0plus,f0plus]=control_law_glover('AGC',delta0plus,omega0plus,e0plus,m0plus,...
                    v0plus, theta0plus, pg0plus,qg0plus,y0plus);
                
elseif strcmp(ControlMode,'LQR')
    
    [pref0plus,f0plus]=control_law_glover('LQR',delta0plus,omega0plus,e0plus,m0plus,...
                    v0plus, theta0plus, pg0plus,qg0plus,[]);
end
                
                



[ deltaDot0plus, omegaDot0plus, eDot0plus, mDot0plus] = gFunctionVectorized( ...
      delta0plus, omega0plus, e0plus,m0plus,...
     vg0plus, thetag0plus,pg0plus,pref0plus,f0plus);
 
 
 z0plus=[delta0plus; omega0plus; e0plus; m0plus; v0plus; theta0plus; pg0plus; qg0plus;pref0plus;f0plus];
 
 
%% 14. Initial conditions for dynamical simulation  (LQR-control)
znew0plus=z0plus(1:prefIdx(1)-1); % initial z
if strcmp(ControlMode,'AGC')
znew0plus=[znew0plus;y0plus];
end

zDot0plus=zeros(size(znew0plus)); % initial zdots
zDot0plus(deltaIdx)=deltaDot0plus; 
zDot0plus(omegaIdx)=omegaDot0plus; 
zDot0plus(eIdx)=eDot0plus; 
zDot0plus(mIdx)=mDot0plus;

if strcmp(ControlMode,'AGC')
zDot0plus(yIdx)=yDot0plus;
end









options.Jacobian={[],Mass};


%%
disp('Running stabilityOPF dynamics'); 
[~,ZNEW]=ode15i(@(t,znew,zDot)...
    runDynamics(t,znew,zDot,NoiseVector), t, znew0plus, zDot0plus, options);
ZNEW = transpose(ZNEW);
disp('Finished stabilityOPF dynamics'); 

%%


[deltaVec, omegaVec, eVec, mVec,yVec,...
    thetaVec, vVec, pgVec, qgVec, ...
    prefVec,fVec, ...
    ploadVec, qloadVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec,yDotVec,ACEVec,...
    PmeasuredVec, OmegaMeasuredVec] =...
    retrieve_output_glover( t, ZNEW , NoiseVector);

[ sanityCheck1,sanityCheck2,sanityCheck3 , success] = sanityCheck(...
    deltaVec, omegaVec, eVec, mVec, ...
    thetaVec, vVec, pgVec, qgVec,...
    prefVec, fVec, ...
    ploadVec,qloadVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec);

[ trCost3] = calculateTrCostUsingIntegration(pgS, qgS, alpha,...
    deltaVec, omegaVec, eVec, mVec, prefVec, fVec, ...
    deltaS, omegaS, eS, mS, prefS, fS,...
    z0);


if exist('Results')~=7
    mkdir('Results'); 
end

cd('Results'); 
if exist(casefile)~=7
    mkdir(casefile);
end

cd(casefile);

if exist(stcontrol)~=7
    mkdir(stcontrol);
end

cd(stcontrol);

if exist(lfcontrol)~=7
    mkdir(lfcontrol);
end

savename=[casefile,'_',st_control,'_',lf_control];
save(savename); 





if exist('figures')~=7
    mkdir('figures');
end
% plots
figx0=0;
figy0=1;
width=8;
height=5;


freqyMin=min(min(omegaVec))./(2*pi);
freqyMax=max(max(omegaVec))./(2*pi);
freqyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'genfreq');
plot(t,omegaVec./(2*pi),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}\boldmath{\omega}$ (Hz)'); 
axis([0 Tfinal freqyMin-freqyOffSet freqyMax+freqyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title(' Generator Frequencies'); 
  print -dpdf freq.pdf
print -depsc2 freq






angleyMin=min(min(deltaVec));
angleyMax=max(max(deltaVec));
angleyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'genangle');
plot(t,deltaVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{\delta}$ (Rad)'); 
axis([0 Tfinal angleyMin-angleyOffSet angleyMax+angleyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator angles'); 
  print -dpdf angles.pdf
print -depsc2 angles




eyMin=min(min(eVec));
eyMax=max(max(eVec));
eyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'gene');
plot(t,eVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{e}$ (pu)'); 
axis([0 Tfinal eyMin-eyOffSet eyMax+eyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator angles'); 
  print -dpdf e.pdf
print -depsc2 e


myMin=min(min(mVec));
myMax=max(max(mVec));
myOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'gene');
plot(t,mVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{m}$ (pu)'); 
axis([0 Tfinal myMin-myOffSet myMax+myOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator angles'); 
  print -dpdf m.pdf
print -depsc2 m





    
vyMin=min(min(vVec));
vyMax=max(max(vVec));
vyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'Voltage mags');
plot(t,vVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{v}$ (pu)'); 
axis([0 Tfinal vyMin-vyOffSet vyMax+vyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator angles'); 
  print -dpdf VoltageMags.pdf
print -depsc2 VoltageMags



thetayMin=min(min(thetaVec));
thetayMax=max(max(thetaVec));
thetayOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'voltage angles');
plot(t,thetaVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{\theta}$ (Rad)'); 
axis([0 Tfinal thetayMin-thetayOffSet thetayMax+thetayOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator angles'); 
  print -dpdf VoltageAngles.pdf
print -depsc2 VoltageAngles





pgyMin=min(min(pgVec));
pgyMax=max(max(pgVec));
pgyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'voltage angles');
plot(t,pgVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{pg}$ (pu)'); 
axis([0 Tfinal pgyMin-pgyOffSet pgyMax+pgyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generated power'); 
  print -dpdf pg.pdf
print -depsc2 pg



qgyMin=min(min(qgVec));
qgyMax=max(max(qgVec));
qgyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'pg');
plot(t,qgVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{qg}$ (pu)'); 
axis([0 Tfinal qgyMin-qgyOffSet qgyMax+qgyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Reactive powers'); 
  print -dpdf qg.pdf
print -depsc2 qg


% 
% 
% qgyMin=min(min(qgVec));
% qgyMax=max(max(qgVec));
% qgyOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'pg');
% plot(t,deltaVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{pg}$ (pu)'); 
% axis([0 Tfinal qgyMin-qgyOffSet qgyMax+qgyOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Reactive powers'); 
%   print -dpdf qg.pdf
% print -depsc2 qg




prefyMin=min(min(prefVec));
prefyMax=max(max(prefVec));
prefyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'pg');
plot(t,prefVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{pref}$ (pu)'); 
axis([0 Tfinal prefyMin-prefyOffSet prefyMax+prefyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator reference'); 
  print -dpdf pref.pdf
print -depsc2 pref



fyMin=min(min(fVec));
fyMax=max(max(fVec));
fyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'pg');
plot(t,fVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{f}$ (pu)'); 
axis([0 Tfinal fyMin-fyOffSet fyMax+fyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator emf'); 
  print -dpdf f.pdf
print -depsc2 f


if strcmp(ControlMode,'AGC')
    
    
    yMin=min(min(yVec));
yMax=max(max(yVec));
yOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'y');
plot(t,yVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{y}$ '); 
axis([0 Tfinal yMin-yOffSet yMax+yOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator AGC signal'); 
  print -dpdf y.pdf
print -depsc2 y
    
    


    ACEMin=min(min(ACEVec));
ACEMax=max(max(ACEVec));
ACEOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'y');
plot(t,ACEVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$ACE$ '); 
axis([0 Tfinal ACEMin-ACEOffSet ACEMax+ACEOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Area ACE signal'); 
  print -dpdf ace.pdf
print -depsc2 ace
end
    

cd(currentdirectory);
end


