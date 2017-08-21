clear all;
clc;
alpha=0;
% casefiles={'case9wmac_con'; 'case14wmac_con';'case39wmac_con'};
casefiles={'case14wmac_con'}; 

lfcontrol='LQR';
fileID=fopen('Results/LQR-NoCoupling.txt','w'); 
fprintf(fileID,'%-20s & %-10s & %-10s & %-20s & %-20s  & %-20s & %-20s \n',...
    'Network', 'Method', 'ss-cost', 'st-cost', 'total cost', 'max freq dev.', 'max volt. dev.');
for case_index=1:length(casefiles)
    casefile=casefiles{case_index};
    lqropf=workflow(casefile,'LQR-OPF',lfcontrol,alpha);
    fprintf(fileID, '%-20s & %-10s & %-10.2f & %-20.2f & %-20.2f  & %-20.4f  %-20.2f\n', ...
    casefile, 'lqr-opf', lqropf.ssCost, lqropf.trCost3, lqropf.ssCost+lqropf.trCost3, ...
    max(max(abs(lqropf.omegaVec-lqropf.OMEGA_S)))./(2*pi), max(max(abs(lqropf.vVec-lqropf.vS))));
     opf=workflow(casefile,'OPF',lfcontrol,alpha);
          fprintf(fileID, '%-20s & %-10s & %-10.2f & %-20.2f & %-20.2f  & %-20.4f  %-20.2f\n', ...
    casefile, 'opf', opf.ssCost, opf.trCost3, opf.ssCost+opf.trCost3, ...
    max(max(abs(opf.omegaVec-opf.OMEGA_S)))./(2*pi), max(max(abs(opf.vVec-opf.vS))));
    
end

fclose(fileID);


