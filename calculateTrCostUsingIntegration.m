function [ trCost3] = calculateTrCostUsingIntegration(pgS, qgS, alpha,...
    deltaVec, omegaVec, eVec, mVec, prefVec, fVec, ...
    deltaS, omegaS, eS, mS, prefS, fS,...
    z0)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
global G n_samples n_pertSamples deltaIdx omegaIdx eIdx mIdx prefIdx fIdx networkS
global Tlqr

[Qinv,Rinv]=QinvRinv(pgS,qgS,alpha, zeros(4*G), zeros(2*G),networkS);
Q=inv(Qinv); 
R=inv(Rinv); 

xMat=[deltaVec;omegaVec;eVec;mVec];
uMat=[prefVec;fVec];
xs=[deltaS;omegaS;eS;mS];
us=[prefS;fS];
xsMat=repmat(xs,1, n_samples-n_pertSamples+1); 
usMat=repmat(us, 1, n_samples-n_pertSamples+1);
x0=z0([deltaIdx;omegaIdx; eIdx;mIdx]);
u0=z0([prefIdx,fIdx]); 
trCost3=Tlqr*calculateTrCostNumeric(xMat,uMat,Q,R, xsMat, usMat); 

end

