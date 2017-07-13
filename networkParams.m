function [ N,G,L,Ymat, Gmat, Bmat,...
    node_set, gen_set, load_set,...
    Cg, yff_vec, yft_vec, ytf_vec, ytt_vec] = networkParams( network )
%NETWORKPARAMS builds network parameters such as bus admittance matrix
% [ Ymat, Gmat, Bmat,...
%     node_set, gen_set, load_set,...
%     Ysh , Cf, Ct,Cg] = networkParams( network )  builds network
%     parameters for the network. The network must be a MATPOWER case
%     struct with appropriate fields.   A similar function to this is
%     incorporated in MATPOWER.  This uses the "from--to model". See
%     networkParamsMNModel for the "mn--model". 
% Description of outputs:
%  1. Ymat: is the complex bus admittance matrix of size(N,N), where N is
%  the number of nodes in the network. 
% 2. Gmat: is the real part of Ymat. 
% 3. Bmat: is the imaginary part of Ymat. 
% 4. node_set: is the set of node labels for the network and has size(N,1).
% 5. gen_set: is the set of generator labels in the network and has size(G,1). 
% 6. load_set: is the set of load labels with no generators and has
% size(L,1). 
% 7. Ysh: is the matrix of node shunt admittances (diagonal) of size(N,N). 
% 8. Cf: is the from matrix of size(B,N), where B is the number of
% branches. 
% 9. Ct: is the to matrix of size(B,N), where B is the number of branches. 
% Cg: is the generator indicator matrix of size(N,G). 
%
% Description of inputs:
% 1. network: is a MATPOWER case struct with appropriate fields. 
% Updates:
% 1. Hafez Bazrafshan 9/9/2016 4:00 PM



node_set=network.bus(:,1); 
gen_set=network.gen(:,1); 
load_set=setdiff(node_set, gen_set); 

N=length(node_set); 
G=length(gen_set); 
L=length(load_set);

fromNodes=network.branch(:,1); 
toNodes=network.branch(:,2); 
B=length(fromNodes);  % number of branches

Cf=zeros(B,N); 
Ct=zeros(B,N); 

Cf(sub2ind([B N], (1:B).', fromNodes))=1;
Ct(sub2ind([B N], (1:B).', toNodes))=1; 

rmn_vec=network.branch(:,3); 
xmn_vec=network.branch(:,4);
zmn_vec=rmn_vec+i*xmn_vec;
ymn_vec=1./zmn_vec;
gmn_vec=real( ymn_vec); 
bmn_vec=imag(ymn_vec);
bcmn_vec=network.branch(:,5); 
taumn_vec=network.branch(:,9); 
taumn_vec(taumn_vec==0)=1;
alphamn_vec=network.branch(:,10);
gs_vec=network.bus(:,5);
bs_vec=network.bus(:,6);
ys_vec=(gs_vec+i*bs_vec)/network.baseMVA; % vector of nodal shunt admittances

yff_vec=(gmn_vec+i*(bmn_vec+bcmn_vec./2))./(taumn_vec.^2); 
yft_vec=-(gmn_vec+i*bmn_vec)./(taumn_vec.*exp(-i*degrees2radians(alphamn_vec))); 
ytf_vec=-(gmn_vec+i*bmn_vec)./(taumn_vec.*exp(i*degrees2radians(alphamn_vec))); 
ytt_vec=gmn_vec+i*(bmn_vec+bcmn_vec/2); 

Yff=diag(yff_vec); 
Yft=diag(yft_vec); 
Ytf=diag(ytf_vec); 
Ytt=diag(ytt_vec); 
Ysh=diag(ys_vec); 

Ymat=Cf.'*Yff*Cf+Cf.'*Yft*Ct+Ct.'*Ytf*Cf+Ct.'*Ytt*Ct+Ysh;
Gmat=real(Ymat); 
Bmat=imag(Ymat); 


Cg=zeros(N,G); 

Cg(sub2ind([N G], gen_set, (1:G).'))=1; 

end

