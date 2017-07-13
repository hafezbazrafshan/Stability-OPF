clear all; 
close all;
clc;
alpha=0;
lfcontrol='LQR';
casefile='case39wmac_con';
lqropf=workflow(casefile,'LQR-OPF',lfcontrol,alpha);
opf=workflow(casefile,'OPF',lfcontrol,alpha);

%%
t=lqropf.t;
Tfinal=lqropf.Tfinal;
G=lqropf.G;
OMEGA_S=lqropf.OMEGA_S;



%% Frequency plots
freqyMax=max(max([lqropf.omegaVec,opf.omegaVec]))./(2*pi);
freqyMin=min(min([lqropf.omegaVec,opf.omegaVec]))./(2*pi);
freqyOffSet=0.05;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:genfreq');
plot(t,lqropf.omegaVec.'/(2*pi),'lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}\bf{\omega}$ (Hz)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal freqyMin-freqyOffSet freqyMax+freqyOffSet+0.1]);
yticks(fig_handle,floor(freqyMin*100./5)*0.05:0.05:ceil(freqyMax*100./5)*0.05);
 GenIdx=cellstr(num2str([1:G].'));
 GenText=cellstr(repmat('Gen. ', G,1));
 legendText=string(GenText)+' '+ string(GenIdx);
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='East';
 grid on;
subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
zoomIndex = (t<=5) & (t>=0);
plot(subfig_handle,t(zoomIndex),lqropf.omegaVec(:,zoomIndex).'/(2*pi)) % plot on new axes
axis(subfig_handle,'tight'); 
% subfig_handle.YTick=...
%     floor(min(min(lqropf.omegaVec))*100./5).*0.05: 0.05: ceil(max(max(lqropf.omegaVec))*100./5).*0.05;
subfig_handle.XTick= 0: 1: 5;


set(subfig_handle,'box','on');
set(subfig_handle,'fontSize',14); 
set(subfig_handle,'defaulttextinterpreter','latex');
subfig_handle.TickLabelInterpreter='latex';
 subfig_handle.XLabel.String='Time (sec)';
 subfig_handle.YLabel.String='$\frac{1}{2\pi}\bf{\omega}$ (Hz)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39LQROPFfreq.pdf
print -depsc2 test39LQROPFfreq
cd('..');



%% Frequency plot 2

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:genfreq');
plot(t,(opf.omegaVec.')/(2*pi),'lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}\bf{\omega}$ (Hz)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal freqyMin-freqyOffSet freqyMax+freqyOffSet+0.1]);
yticks(fig_handle,floor(freqyMin*100./5)*0.05:0.05:ceil(freqyMax*100./5)*0.05);
 GenIdx=cellstr(num2str([1:G].'));
 GenText=cellstr(repmat('Gen. ', G,1));
 legendText=string(GenText)+' '+ string(GenIdx);
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='East';
 grid on;

 
 subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
zoomIndex = (t<=5) & (t>=0);
plot(subfig_handle,t(zoomIndex),opf.omegaVec(:,zoomIndex).'/(2*pi)) % plot on new axes
axis(subfig_handle,'tight'); 
% subfig_handle.YTick=...
%     floor(min(min(opf.omegaVec))*100./5).*0.05: 0.05: ceil(max(max(opf.omegaVec))*100./5).*0.05;subfig_handle.XTick= 0: 1: 5;
subfig_handle.XTick= 0: 1: 5;

set(subfig_handle,'box','on');
set(subfig_handle,'fontSize',14); 
set(subfig_handle,'defaulttextinterpreter','latex');
subfig_handle.TickLabelInterpreter='latex';
 subfig_handle.XLabel.String='Time (sec)';
 subfig_handle.YLabel.String='$\frac{1}{2\pi}\bf{\omega}$ (Hz)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)
 
  cd('Figures'); 
  print -dpdf test39OPFfreq.pdf
print -depsc2 test39OPFfreq
cd('..');



% %% voltage plots
% vy=max(max([abs(lqropf.vVec-lqropf.vS),abs(opf.vVec-opf.vS)]));
% vyOffSet=0.05;
% 
% x0=0;
% y0=1;
% width=8;
% height=5;
% figure1=figure('Units','inches',...
% 'Position',[x0 y0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'LQR-OPF:busvolt');
% plot(t,(lqropf.vVec-lqropf.vS).','lineWidth',2);
% fig_handle=gca; 
% set(fig_handle,'box','on');
% set(fig_handle,'fontSize',20); 
% set(fig_handle,'defaulttextinterpreter','latex');
% grid on;
% xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\mathbf{v}-\mathbf{v}^\mathrm{eq}$ (pu)'); 
% fig_handle.TickLabelInterpreter='latex';
% axis(fig_handle,[0 Tfinal -vy-vyOffSet vy+vyOffSet]);
% yticks(fig_handle,floor(-vy*100./5)*0.05:0.05:ceil(vy*100./5)*0.05);
%  
%  
%  
%  if exist('Figures')~=7
%     mkdir('Figures'); 
%  end
% 
%  cd('Figures'); 
%   print -dpdf2 test39LQROPFvolt.pdf
% print -depsc2 test39LQROPFvolt
% cd('..');
% 
% 
% 
% %% voltage plots
% 
% 
% x0=0;
% y0=1;
% width=8;
% height=5;
% figure1=figure('Units','inches',...
% 'Position',[x0 y0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'LQR-OPF:busvolt');
% plot(t,(opf.vVec-opf.vS).','lineWidth',2);
% fig_handle=gca; 
% set(fig_handle,'box','on');
% set(fig_handle,'fontSize',20); 
% set(fig_handle,'defaulttextinterpreter','latex');
% grid on;
% xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\mathbf{v}-\mathbf{v}^\mathrm{eq}$ (pu)'); 
% fig_handle.TickLabelInterpreter='latex';
% axis(fig_handle,[0 Tfinal -vy-vyOffSet vy+vyOffSet]);
% yticks(fig_handle,floor(-vy*100./5)*0.05:0.05:ceil(vy*100./5)*0.05);
%  
%  
% 
%  
%  if exist('Figures')~=7
%     mkdir('Figures'); 
%  end
% 
%  cd('Figures'); 
%   print -dpdf test39OPFvolt.pdf
% print -depsc2 test39OPFvolt
% cd('..');


%% mech plots
my=max(max([abs(lqropf.mVec-lqropf.mS),abs(opf.mVec-opf.mS)]));
myOffSet=0.05;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:mech');
plot(t,(lqropf.mVec-lqropf.mS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal -my-myOffSet my+myOffSet+2]);
yticks(fig_handle,floor(-my):1:ceil(my));
 GenIdx=cellstr(num2str([1:G].'));
 GenText=cellstr(repmat('Gen. ', G,1));
 legendText=string(GenText)+' '+ string(GenIdx);
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='East';
 grid on;
 


  subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
zoomIndex = (t<=5) & (t>=0);
plot(subfig_handle,t(zoomIndex),(lqropf.mVec(:,zoomIndex)-lqropf.mS).') % plot on new axes
axis(subfig_handle,'tight'); 
subfig_handle.XTick= 0: 1: 5;

set(subfig_handle,'box','on');
set(subfig_handle,'fontSize',14); 
set(subfig_handle,'defaulttextinterpreter','latex');
subfig_handle.TickLabelInterpreter='latex';
 subfig_handle.XLabel.String='Time (sec)';
 subfig_handle.YLabel.String='$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)

 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39LQROPFmech.pdf
print -depsc2 test39LQROPFmech
cd('..');




%% mech plots
my=max(max([abs(lqropf.mVec-lqropf.mS),abs(opf.mVec-opf.mS)]));
myOffSet=0.05;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:mech');
plot(t,(opf.mVec-opf.mS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal -my-myOffSet my+myOffSet+2]);
yticks(fig_handle,floor(-my):1:ceil(my));
 GenIdx=cellstr(num2str([1:G].'));
 GenText=cellstr(repmat('Gen. ', G,1));
 legendText=string(GenText)+' '+ string(GenIdx);
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='East';
 grid on;
 
 
 
  subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
zoomIndex = (t<=5) & (t>=0);
plot(subfig_handle,t(zoomIndex),(opf.mVec(:,zoomIndex)-opf.mS).') % plot on new axes
axis(subfig_handle,'tight'); 
% subfig_handle.YTick= -0.1: 0.05: 0.1;
subfig_handle.XTick= 0: 1: 5;

set(subfig_handle,'box','on');
set(subfig_handle,'fontSize',14); 
set(subfig_handle,'defaulttextinterpreter','latex');
subfig_handle.TickLabelInterpreter='latex';
 subfig_handle.XLabel.String='Time (sec)';
 subfig_handle.YLabel.String='$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)


 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39OPFmech.pdf
print -depsc2 test39OPFmech
cd('..');



%% pg plots
pgy=max(max([abs(lqropf.pgVec-lqropf.pgS),abs(opf.pgVec-opf.pgS)]));
pgyOffSet=0.05;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:pg');
plot(t,(lqropf.pgVec-lqropf.pgS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{p}_{\mathrm{g}}-\mathbf{p}_{\mathrm{g}}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal -pgy-pgyOffSet pgy+pgyOffSet+2]);
yticks(fig_handle,floor(-pgy):1:ceil(pgy));
 GenIdx=cellstr(num2str([1:G].'));
 GenText=cellstr(repmat('Gen. ', G,1));
 legendText=string(GenText)+' '+ string(GenIdx);
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='East';
 grid on;
 


  subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
zoomIndex = (t<=5) & (t>=0);
plot(subfig_handle,t(zoomIndex),(lqropf.pgVec(:,zoomIndex)-lqropf.pgS).') % plot on new axes
axis(subfig_handle,'tight'); 
subfig_handle.XTick= 0: 1: 5;

set(subfig_handle,'box','on');
set(subfig_handle,'fontSize',14); 
set(subfig_handle,'defaulttextinterpreter','latex');
subfig_handle.TickLabelInterpreter='latex';
 subfig_handle.XLabel.String='Time (sec)';
 subfig_handle.YLabel.String='$\mathbf{p}_{\mathrm{g}}-\mathbf{p}_{\mathrm{g}}^\mathrm{eq}$ (pu)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)

 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39LQROPFpg.pdf
print -depsc2 test39LQROPFpg
cd('..');




%% pg plots
my=max(max([abs(lqropf.pgVec-lqropf.pgS),abs(opf.pgVec-opf.pgS)]));
myOffSet=0.05;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:mech');
plot(t,(opf.pgVec-opf.pgS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{p}_{\mathrm{g}}-\mathbf{p}_{\mathrm{g}}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal -pgy-pgyOffSet pgy+pgyOffSet+2]);
yticks(fig_handle,floor(-pgy):1:ceil(pgy));
 GenIdx=cellstr(num2str([1:G].'));
 GenText=cellstr(repmat('Gen. ', G,1));
 legendText=string(GenText)+' '+ string(GenIdx);
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='East';
 grid on;
 
 
 
  subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
zoomIndex = (t<=5) & (t>=0);
plot(subfig_handle,t(zoomIndex),(opf.pgVec(:,zoomIndex)-opf.pgS).') % plot on new axes
axis(subfig_handle,'tight'); 
% subfig_handle.YTick= -0.1: 0.05: 0.1;
subfig_handle.XTick= 0: 1: 5;

set(subfig_handle,'box','on');
set(subfig_handle,'fontSize',14); 
set(subfig_handle,'defaulttextinterpreter','latex');
subfig_handle.TickLabelInterpreter='latex';
 subfig_handle.XLabel.String='Time (sec)';
 subfig_handle.YLabel.String='$\mathbf{p}_{\mathrm{g}}-\mathbf{p}_{\mathrm{g}}^\mathrm{eq}$ (pu)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)


 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39OPFpg.pdf
print -depsc2 test39OPFpg
cd('..');

%% e plots
eyMax=max(max([lqropf.eVec-lqropf.eS,opf.eVec-opf.eS]));
eyMin=min(min([lqropf.eVec-lqropf.eS,opf.eVec-opf.eS]));

eyOffSet=0.01;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:e');
plot(t,(lqropf.eVec-lqropf.eS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{e}-\mathbf{e}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal eyMin-eyOffSet eyMax+eyOffSet]);
yticks(fig_handle,floor(eyMin*100/5)*0.05:0.05:ceil(eyMax*100/5)*0.05);
%  GenIdx=cellstr(num2str([1:G].'));
%  GenText=cellstr(repmat('Gen. ', G,1));
%  legendText=string(GenText)+' '+ string(GenIdx);
%  lgd=legend(fig_handle,legendText);
%  lgd.FontSize=12;
%  lgd.Interpreter='Latex';
%  lgd.Location='East';
 grid on;
 
 
 
%   subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
% zoomIndex = (t<=5) & (t>=0);
% plot(subfig_handle,t(zoomIndex),(lqropf.eVec(:,zoomIndex)-lqropf.eS).') % plot on new axes
% axis(subfig_handle,'tight'); 
% % subfig_handle.YTick= -0.1: 0.05: 0.1;
% subfig_handle.XTick= 0: 1: 5;
% 
% set(subfig_handle,'box','on');
% set(subfig_handle,'fontSize',14); 
% set(subfig_handle,'defaulttextinterpreter','latex');
% subfig_handle.TickLabelInterpreter='latex';
%  subfig_handle.XLabel.String='Time (sec)';
%  subfig_handle.YLabel.String='$\mathbf{e}-\mathbf{e}^\mathrm{eq}$ (pu)'; 
%  x=[0.21 0.16];
%  y=[0.55 0.49];
%  annotation('textarrow',x,y)


 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39LQROPFe.pdf
print -depsc2 test39LQROPFe
cd('..');



%% e plots
eyMax=max(max([lqropf.eVec-lqropf.eS,opf.eVec-opf.eS]));
eyMin=min(min([lqropf.eVec-lqropf.eS,opf.eVec-opf.eS]));

eyOffSet=0.01;


x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:mech');
plot(t,(opf.eVec-opf.eS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{e}-\mathbf{e}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal eyMin-eyOffSet eyMax+eyOffSet]);
yticks(fig_handle,floor(eyMin*100/5)*0.05:0.05:ceil(eyMax*100/5)*0.05);
%  GenIdx=cellstr(num2str([1:G].'));
%  GenText=cellstr(repmat('Gen. ', G,1));
%  legendText=string(GenText)+' '+ string(GenIdx);
%  lgd=legend(fig_handle,legendText);
%  lgd.FontSize=12;
%  lgd.Interpreter='Latex';
%  lgd.Location='East';
 grid on;
 
 
 
%   subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
% zoomIndex = (t<=5) & (t>=0);
% plot(subfig_handle,t(zoomIndex),(opf.eVec(:,zoomIndex)-opf.eS).') % plot on new axes
% axis(subfig_handle,'tight'); 
% % subfig_handle.YTick= -0.1: 0.05: 0.1;
% subfig_handle.XTick= 0: 1: 5;

% set(subfig_handle,'box','on');
% set(subfig_handle,'fontSize',14); 
% set(subfig_handle,'defaulttextinterpreter','latex');
% subfig_handle.TickLabelInterpreter='latex';
%  subfig_handle.XLabel.String='Time (sec)';
%  subfig_handle.YLabel.String='$\mathbf{e}-\mathbf{e}^\mathrm{eq}$ (pu)'; 
%  x=[0.21 0.16];
%  y=[0.55 0.49];
%  annotation('textarrow',x,y)


 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39OPFe.pdf
print -depsc2 test39OPFe
cd('..');

%% f plots
fyMax=max(max([lqropf.fVec-lqropf.fS,opf.fVec-opf.fS]));
fyMin=min(min([lqropf.fVec-lqropf.fS,opf.fVec-opf.fS]));

fyOffSet=0.01;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:f');
plot(t,(lqropf.fVec-lqropf.fS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{f}-\mathbf{f}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal fyMin-fyOffSet fyMax+fyOffSet]);
yticks(fig_handle,floor(fyMin*100/5)*0.05-0.05:0.1:ceil(fyMax*100/5)*0.05+0.05);
 GenIdx=cellstr(num2str([1:G].'));
 GenText=cellstr(repmat('Gen. ', G,1));
%  legendText=string(GenText)+' '+ string(GenIdx);
%  lgd=legend(fig_handle,legendText);
%  lgd.FontSize=12;
%  lgd.Interpreter='Latex';
%  lgd.Location='East';
 grid on;
 
 
 
%   subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
% zoomIndex = (t<=5) & (t>=0);
% plot(subfig_handle,t(zoomIndex),(lqropf.eVec(:,zoomIndex)-lqropf.eS).') % plot on new axes
% axis(subfig_handle,'tight'); 
% % subfig_handle.YTick= -0.1: 0.05: 0.1;
% subfig_handle.XTick= 0: 1: 5;
% 
% set(subfig_handle,'box','on');
% set(subfig_handle,'fontSize',14); 
% set(subfig_handle,'defaulttextinterpreter','latex');
% subfig_handle.TickLabelInterpreter='latex';
%  subfig_handle.XLabel.String='Time (sec)';
%  subfig_handle.YLabel.String='$\mathbf{e}-\mathbf{e}^\mathrm{eq}$ (pu)'; 
%  x=[0.21 0.16];
%  y=[0.55 0.49];
%  annotation('textarrow',x,y)


 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39LQROPFf.pdf
print -depsc2 test39LQROPFf
cd('..');



%% e plots
fyMax=max(max([lqropf.fVec-lqropf.fS,opf.fVec-opf.fS]));
fyMin=min(min([lqropf.fVec-lqropf.fS,opf.fVec-opf.fS]));

fyOffSet=0.01;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:field');
plot(t,(opf.fVec-opf.fS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{f}-\mathbf{f}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal fyMin-fyOffSet fyMax+fyOffSet]);
yticks(fig_handle,floor(fyMin*100/5)*0.05-0.05:0.1:ceil(fyMax*100/5)*0.05+0.05);
%  GenIdx=cellstr(num2str([1:G].'));
%  GenText=cellstr(repmat('Gen. ', G,1));
%  legendText=string(GenText)+' '+ string(GenIdx);
%  lgd=legend(fig_handle,legendText);
%  lgd.FontSize=12;
%  lgd.Interpreter='Latex';
%  lgd.Location='East';
 grid on;
 
 
 



 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39OPFf.pdf
print -depsc2 test39OPFf
cd('..');




%% delta plots
deltayMax=max(max(radians2degrees([lqropf.deltaVec-lqropf.deltaS,opf.deltaVec-opf.deltaS])));
deltayMin=min(min(radians2degrees([lqropf.deltaVec-lqropf.deltaS,opf.deltaVec-opf.deltaS])));

deltayOffSet=0.01;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:genangles');
plot(t,radians2degrees(lqropf.deltaVec-lqropf.deltaS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\delta-\delta^{\mathrm{eq}}$ (deg.)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal deltayMin-deltayOffSet deltayMax+deltayOffSet]);
yticks(fig_handle,floor(deltayMin./5)*5:5:ceil(deltayMax./5)*5);
%  GenIdx=cellstr(num2str([1:G].'));
%  GenText=cellstr(repmat('Gen. ', G,1));
%  legendText=string(GenText)+' '+ string(GenIdx);
%  lgd=legend(fig_handle,legendText);
%  lgd.FontSize=12;
%  lgd.Interpreter='Latex';
%  lgd.Location='East'
 grid on;
%  grid on;
% subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
% zoomIndex = (t<=5) & (t>=0);
% plot(subfig_handle,t(zoomIndex),radians2degrees(lqropf.deltaVec(:,zoomIndex)-lqropf.deltaS).'/(2*pi)) % plot on new axes
% axis(subfig_handle,'tight'); 
% subfig_handle.YTick= -0.05: 0.05: 0.05;
% subfig_handle.XTick= 0: 1: 5;
% 
% 
% set(subfig_handle,'box','on');
% set(subfig_handle,'fontSize',14); 
% set(subfig_handle,'defaulttextinterpreter','latex');
% subfig_handle.TickLabelInterpreter='latex';
%  subfig_handle.XLabel.String='Time (sec)';
%  subfig_handle.YLabel.String='$\delta-\delta^{\mathrm{eq}}$ (deg.)'; 
%  x=[0.21 0.16];
%  y=[0.55 0.49];
%  annotation('textarrow',x,y)

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39LQROPFdelta.pdf
print -depsc2 test39LQROPFdelta
cd('..');



%% delta plots
deltayMax=max(max(radians2degrees([lqropf.deltaVec-lqropf.deltaS,opf.deltaVec-opf.deltaS])));
deltayMin=min(min(radians2degrees([lqropf.deltaVec-lqropf.deltaS,opf.deltaVec-opf.deltaS])));

deltayOffSet=0.01;
x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:genangles');
plot(t,radians2degrees(opf.deltaVec-opf.deltaS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\delta-\delta^{\mathrm{eq}}$ (deg.)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal deltayMin-deltayOffSet deltayMax+deltayOffSet]);
yticks(fig_handle,floor(deltayMin./5)*5:5:ceil(deltayMax./5)*5);
%  GenIdx=cellstr(num2str([1:G].'));
%  GenText=cellstr(repmat('Gen. ', G,1));
%  legendText=string(GenText)+' '+ string(GenIdx);
%  lgd=legend(fig_handle,legendText);
%  lgd.FontSize=12;
%  lgd.Interpreter='Latex';
%  lgd.Location='East'
 grid on;
%  grid on;
% subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
% zoomIndex = (t<=5) & (t>=0);
% plot(subfig_handle,t(zoomIndex),radians2degrees(lqropf.deltaVec(:,zoomIndex)-lqropf.deltaS).'/(2*pi)) % plot on new axes
% axis(subfig_handle,'tight'); 
% subfig_handle.YTick= -0.05: 0.05: 0.05;
% subfig_handle.XTick= 0: 1: 5;
% 
% 
% set(subfig_handle,'box','on');
% set(subfig_handle,'fontSize',14); 
% set(subfig_handle,'defaulttextinterpreter','latex');
% subfig_handle.TickLabelInterpreter='latex';
%  subfig_handle.XLabel.String='Time (sec)';
%  subfig_handle.YLabel.String='$\delta-\delta^{\mathrm{eq}}$ (deg.)'; 
%  x=[0.21 0.16];
%  y=[0.55 0.49];
%  annotation('textarrow',x,y)

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39OPFdelta.pdf
print -depsc2 test39OPFdelta
cd('..');



%% pref plots
prefyMax=max(max([lqropf.prefVec-lqropf.prefS,opf.prefVec-opf.prefS]));
prefyMin=min(min([lqropf.prefVec-lqropf.prefS,opf.prefVec-opf.prefS]));

prefyOffSet=0.05;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:pref');
plot(t,(lqropf.prefVec-lqropf.prefS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{r}-\mathbf{r}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal prefyMin-prefyOffSet prefyMax+prefyOffSet+2]);
yticks(fig_handle,floor(prefyMin):1:ceil(prefyMax));
 GenIdx=cellstr(num2str([1:G].'));
 GenText=cellstr(repmat('Gen. ', G,1));
 legendText=string(GenText)+' '+ string(GenIdx);
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='East';
 grid on;
 
 
 
  subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
zoomIndex = (t<=5) & (t>=0);
plot(subfig_handle,t(zoomIndex),(lqropf.prefVec(:,zoomIndex)-lqropf.prefS).') % plot on new axes
axis(subfig_handle,'tight'); 
% subfig_handle.YTick= -0.1: 0.05: 0.1;
subfig_handle.XTick= 0: 1: 5;

set(subfig_handle,'box','on');
set(subfig_handle,'fontSize',14); 
set(subfig_handle,'defaulttextinterpreter','latex');
subfig_handle.TickLabelInterpreter='latex';
 subfig_handle.XLabel.String='Time (sec)';
 subfig_handle.YLabel.String='$\mathbf{r}-\mathbf{r}^\mathrm{eq}$ (pu)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)


 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39LQROPFpref.pdf
print -depsc2 test39LQROPFpref
cd('..');


%%


x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:pref');
plot(t,(opf.prefVec-opf.prefS).','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{r}-\mathbf{r}^\mathrm{eq}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal prefyMin-prefyOffSet prefyMax+prefyOffSet+2]);
yticks(fig_handle,floor(prefyMin):1:ceil(prefyMax));
 GenIdx=cellstr(num2str([1:G].'));
 GenText=cellstr(repmat('Gen. ', G,1));
 legendText=string(GenText)+' '+ string(GenIdx);
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='East';
 grid on;
 
 
 
  subfig_handle=axes('position',[0.26 0.6 0.3 .3]);
zoomIndex = (t<=5) & (t>=0);
plot(subfig_handle,t(zoomIndex),(opf.prefVec(:,zoomIndex)-opf.prefS).') % plot on new axes
axis(subfig_handle,'tight'); 
% subfig_handle.YTick= -0.1: 0.05: 0.1;
subfig_handle.XTick= 0: 1: 5;

set(subfig_handle,'box','on');
set(subfig_handle,'fontSize',14); 
set(subfig_handle,'defaulttextinterpreter','latex');
subfig_handle.TickLabelInterpreter='latex';
 subfig_handle.XLabel.String='Time (sec)';
 subfig_handle.YLabel.String='$\mathbf{r}-\mathbf{r}^\mathrm{eq}$ (pu)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)


 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39OPFpref.pdf
print -depsc2 test39OPFpref
cd('..');




%% f plots
vyMax=max(max([lqropf.vVec,opf.vVec]));
vyMin=min(min([lqropf.vVec,opf.vVec]))
vyOffSet=0.02;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:v');
plot(t,lqropf.vVec.','lineWidth',2);
% hold on;
% vMaxPlot=plot(t,repmat(1.06,length(t),1),'k--','lineWidth',2);
% hold on
% vMinPlot=plot(t,repmat(0.94,length(t),1),'k--','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{v}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal vyMin-vyOffSet vyMax+vyOffSet]);
yticks(fig_handle,floor(vyMin*100/5)*0.05:0.02:ceil(vyMax*100/5)*0.05);
% legLines=[vMaxPlot vMinPlot];
%  legendText={'$\mathbf{v}_{\max}$';'$\mathbf{v}_{\min}$'};
%  lgd=legend(legLines,legendText);
%  lgd.FontSize=20;
%  lgd.Interpreter='Latex';
%  lgd.Location='South';
%   lgd.Orientation='Horizontal';
 grid on;
 
 



 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39LQROPFvolt.pdf
print -depsc2 test39LQROPFvolt
cd('..');




%% f plots

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'OPF:v');
plot(t,opf.vVec.','lineWidth',2);
% hold on;
% vMaxPlot=plot(t,repmat(1.06,length(t),1),'k--','lineWidth',2);
% hold on
% vMinPlot=plot(t,repmat(0.94,length(t),1),'k--','lineWidth',2);
fig_handle=gca; 
set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{v}$ (pu)'); 
fig_handle.TickLabelInterpreter='latex';
axis(fig_handle,[0 Tfinal vyMin-vyOffSet vyMax+vyOffSet]);
yticks(fig_handle,floor(vyMin*100/5)*0.05:0.02:ceil(vyMax*100/5)*0.05);
% legLines=[vMaxPlot vMinPlot];
%  legendText={'$\mathbf{v}_{\max}$';'$\mathbf{v}_{\min}$'};
%  lgd=legend(legLines,legendText);
%  lgd.FontSize=20;
%  lgd.Interpreter='Latex';
%  lgd.Location='South';
%  lgd.Orientation='Horizontal';
 grid on;
 
 



 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf test39OPFvolt.pdf
print -depsc2 test39OPFvolt
cd('..');

