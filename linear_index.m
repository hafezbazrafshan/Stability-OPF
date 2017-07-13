function [ branchNumber] = linear_index(...
    fromNode, toNode , network)
%LINEAR_INDEX finds the linear index of a branch
%   branchNumber=linear_index(
% fromNode,toNode,network) finds the linear index
% of the branch associated with the edge
% (fromNode,toNode) for the input  structure 
% network
%
% Description of outputs:
% 1. branchNumber: scalar corresponding to the branch number
% 
% Description of inputs:
% fromNode: the node that link initiates from
% toNode: the node that link terminates at
% network: the network structure created by MATPOWER
% Updates: 
% Hafez Bazrafshan 8/29/2016 1:07 PM

branchNumber=find( and(network.branch(:,1)==...
    fromNode,network.branch(:,2)==toNode));

end

