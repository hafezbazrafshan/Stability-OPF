function [ trCost] = calculateTrCostNumeric(x, u, Q, R , xs, us)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global fsample 
x=x-xs; 
u=u-us;

trCost=(1./2/fsample)*sum(diag(x.'*Q*x+u.'*R*u)); 

 
end