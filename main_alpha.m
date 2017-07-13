clear all;
clc;
alphaVec=[0;0.2;0.4;0.6;0.8];
lqropfVec=cell(size(alphaVec));
opfVec=cell(size(alphaVec));

casefile='case39wmac_con';
lfcontrol='LQR';
fileID=fopen('Results/LQR-comparisonAlpha.txt','w'); 
fprintf(fileID,'%-10s & %-10s & %-10s & %-20s & %-20s  & %-20s & %-20s \n',...
  'alpha3', 'Method', 'ss-cost', 'st-cost', 'total cost', 'max freq dev.', 'max volt. dev.');
for ii=1:length(alphaVec)
    alpha=alphaVec(ii);
    lqropf=workflow(casefile,'LQR-OPF',lfcontrol,alpha);
    lqropfVec{ii}=lqropf;
    fprintf(fileID, '%-10.2f  & %-10s & %-10.2f & %-20.2f & %-20.2f  & %-20.4f  %-20.2f\n', ...
     alpha,'lqr-opf', lqropf.ssCost, lqropf.trCost3, lqropf.ssCost+lqropf.trCost3, ...
    max(max(abs(lqropf.omegaVec-lqropf.OMEGA_S)))./(2*pi), max(max(abs(lqropf.vVec-lqropf.vS))));
     opf=workflow(casefile,'OPF',lfcontrol,alpha);
     opfVec{ii}=opf;
          fprintf(fileID, '%-10.2f & %-10s & %-10.2f & %-20.2f & %-20.2f  & %-20.4f  %-20.2f\n', ...
     alpha,'opf', opf.ssCost, opf.trCost3, opf.ssCost+opf.trCost3, ...
    max(max(abs(opf.omegaVec-opf.OMEGA_S)))./(2*pi), max(max(abs(opf.vVec-opf.vS))));
    
end

fclose(fileID);

