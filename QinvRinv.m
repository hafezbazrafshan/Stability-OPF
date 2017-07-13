function [ Qinv,Rinv ] = QinvRinv( pgS,qgS,alpha, Qinv, Rinv,networkS )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


global Sbase G omegaIdx deltaIdx eIdx mIdx freqR_vec







% Qinv(sub2ind([4*G,4*G], omegaIdx, omegaIdx)) =0.001*((-alpha*pgS + networkS.gen(:,9)./Sbase)./(networkS.gen(:,9)./Sbase));
% Qinv(sub2ind([4*G,4*G], omegaIdx, omegaIdx)) =0.01*((-alpha*pgS+networkS.gen(:,9)./Sbase)./networkS.gen(:,9)./Sbase);
% Qinv(sub2ind([4*G,4*G], omegaIdx, omegaIdx)) =0.0001;
% Qinv(sub2ind([4*G,4*G], deltaIdx, deltaIdx)) =1*( 0*pgS + networkS.gen(:,9)./Sbase);
Qinv(sub2ind([4*G,4*G], omegaIdx, omegaIdx))=1*(1-alpha*pgS./(networkS.gen(:,9)./Sbase));
Qinv(sub2ind([4*G,4*G], deltaIdx, deltaIdx)) =1*(1-alpha*pgS./(networkS.gen(:,9)./Sbase));
Qinv(sub2ind([4*G, 4*G], eIdx,eIdx)) =1*(1-alpha*qgS./(networkS.gen(:,4)./Sbase));
Qinv(sub2ind([4*G,4*G],mIdx,mIdx))=1*(1-alpha*pgS./(networkS.gen(:,9)./Sbase));
Rinv(sub2ind([2*G, 2*G], 1:G,1:G))  =1*(1-alpha*pgS./(networkS.gen(:,9)./Sbase));
Rinv(sub2ind([2*G, 2*G], G+1:2*G, G+1:2*G)) =1*(1-alpha*qgS./(networkS.gen(:,4)./Sbase));

end

