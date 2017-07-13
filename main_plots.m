%% Frequency plots
freqyMin=min(min([a.omegaVec, b.omegaVec]))./(2*pi);
freqyMax=max(max([a.omegaVec,b.omegaVec]))./(2*pi);
freqyOffSet=0;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:genfreq');
plot(t,a.ZNEW(omegaIdx,:).'/(2*pi),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}\boldmath{\omega}$ (Hz)'); 
axis([0 Tfinal freqyMin-freqyOffSet freqyMax+freqyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('LQR-OPF: Generator Frequencies'); 
  print -dpdf figures/test39LQROPFfreq.pdf
print -depsc2 figures/test39LQROPFfreq

% 
% % 
figure2=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure2, 'Name', 'OPF:genfreq');
plot(t,b.ZNEW(omegaIdx,:).'/(2*pi),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}$\boldmath$\omega$ (Hz)'); 
axis([0 Tfinal freqyMin-freqyOffSet freqyMax+freqyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('OPF: Generator Frequencies'); 
 print -dpdf figures/test39OPFfreq.pdf
print -depsc2 figures/test39OPFfreq


% 
% %% Control inputs
myMin=min(min([a.mVec, b.mVec]));
myMax=max(max([a.mVec,b.mVec]));



x0=0;
y0=1;
width=8;
height=5;
figure3=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure3, 'Name', 'LQR-OPF:mech');
plot(t,a.mVec.','lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$ \mathbf{m} $ (pu)'); 
 axis([0 Tfinal myMin-1 myMax+1]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('LQR-OPF: Generator Mechanical Power Input'); 
  print -dpdf figures/test39LQROPFdeltaM.pdf
print -depsc2 figures/test39LQROPFdeltaM





x0=0;
y0=1;
width=8;
height=5;
figure4=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure4, 'Name', 'OPF:mech');
plot(t,b.mVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$ \mathbf{m}$ (pu)'); 
  axis([0 Tfinal myMin-1 myMax+1]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
 title('OPF: Generator Mechanical Power Input'); 
  print -dpdf figures/test39OPFdeltaM.pdf
print -depsc2 figures/test39OPFdeltaM



fyMin=min(min([a.fVec, b.fVec]));
fyMax=max(max([a.fVec,b.fVec]));





x0=0;
y0=1;
width=8;
height=5;
figure5=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure5, 'Name', 'LQR-OPF:field');
plot(t,a.fVec.','lineWidth',2);
axis([0 Tfinal fyMin-1 fyMax+1]); 
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$ \mathbf{f} $ (pu)'); 
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('LQR-OPF: Generator Field Voltage'); 
  print -dpdf figures/test39LQROPFdeltaF.pdf
print -depsc2 figures/test39LQROPFdeltaF




x0=0;
y0=1;
width=8;
height=5;
figure6=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure6, 'Name', 'OPF:field');
plot(t,b.fVec.','lineWidth',2);
axis([0 Tfinal fyMin-1 fyMax+1]); 
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{f}$ (pu)'); 
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
 title('OPF: Generator Field Voltage'); 
  print -dpdf figures/test39OPFdeltaF.pdf
print -depsc2 figures/test39OPFdeltaF

 %% Voltage profile
 vyMin=min(min([a.vVec, b.vVec]));
vyMax=max(max([a.vVec,b.vVec]));
 
 
 x0=0;
y0=1;
width=8;
height=5;
figure7=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure7, 'Name', 'LQR-OPF:voltageprofile');
plot(t,a.ZNEW(vIdx,:).','lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{v} $ (pu)'); 
 axis([0 Tfinal vyMin-0.1 vyMax+0.1]); 
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('LQR-OPF: Time-varying nodal voltages'); 
  print -dpdf figures/test39LQROPFvoltage.pdf
print -depsc2 figures/test39LQROPFvoltage


 x0=0;
y0=1;
width=8;
height=5;
figure8=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure8, 'Name', 'OPF:voltageprofile');
plot(t,b.ZNEW(vIdx,:).','lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{v} $ (pu)'); 
  axis([0 Tfinal vyMin-0.1 vyMax+0.1]); 
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('OPF: Time-varying nodal voltages'); 
  print -dpdf figures/test39OPFvoltage.pdf
print -depsc2 figures/test39OPFvoltage
