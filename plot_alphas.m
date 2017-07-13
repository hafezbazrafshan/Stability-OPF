load('Results/alpha_results');
x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'barcomparison');
lqropfCoordinates=alphaVec-0.03;
lqropfMat = [lqropfVec{1}.ssCost lqropfVec{1}.trCost3; ...
    lqropfVec{2}.ssCost, lqropfVec{2}.trCost3;...
    lqropfVec{3}.ssCost, lqropfVec{3}.trCost3;...
    lqropfVec{4}.ssCost, lqropfVec{4}.trCost3; 
    lqropfVec{5}.ssCost, lqropfVec{5}.trCost3];
lqropfBarPlot=bar(lqropfCoordinates,  lqropfMat,'stacked');
set(lqropfBarPlot,  'BarWidth',0.3);
lqropfBarPlot(1).FaceColor=[0 0 0.5]; % blue
lqropfBarPlot(2).FaceColor=	[100,206,210]/255; % cyan
 
set(gca,'nextplot','add') ;
opfCoordinates=alphaVec+0.03;
opfMat = [opfVec{1}.ssCost opfVec{1}.trCost3; ...
    opfVec{2}.ssCost, opfVec{2}.trCost3;...
    opfVec{3}.ssCost,opfVec{3}.trCost3;...
    opfVec{4}.ssCost, opfVec{4}.trCost3; 
   opfVec{5}.ssCost, opfVec{5}.trCost3];

opfBarPlot=bar(opfCoordinates,   opfMat,'stacked');
set(opfBarPlot,  'BarWidth',0.3);
opfBarPlot(1).FaceColor=[0 0 255]/255; % blue
opfBarPlot(2).FaceColor=[135,206,250]/255; %light sky blue

fig_handle=gca;

fig_handle.XTick=alphaVec;
fig_handle.YLim=[0  max( max( [sum(lqropfMat,2); sum(opfMat,2)]))+40000];
fig_handle.YTick=[0:10000:80000];
fig_handle.YTickLabel=[num2str(fig_handle.YTick.'/1000),repmat('k',length(fig_handle.YTick.'),1)];


set(gca,'nextplot','add') ;
hold on
linePlot=plot(alphaVec,sum(opfMat,2)-sum(lqropfMat,2));
linePlot.LineStyle='-';
linePlot.LineWidth=3;
linePlot.Color=[1 0.7 0];
linePlot.Marker='s';
linePlot.MarkerSize=10;
linePlot.MarkerFaceColor='w';


set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Coupling coefficient $\alpha$', 'FontWeight','bold');
 ylabel('Cost (\$)'); 
fig_handle.TickLabelInterpreter='latex';

yticks=fig_handle.YTick.';
legendText={'LQR-OPF steady-state costs'; 'LQR-OPF stability costs';'OPF steady-state costs'; 'OPF stability costs'; 'LQR-OPF savings'};
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='North';
 lgd.Orientation='Vertical';

 
  cd('Figures'); 
  print -dpdf test39LQROPFbarplot.pdf
print -depsc2 test39LQROPFbarplot
cd('..');