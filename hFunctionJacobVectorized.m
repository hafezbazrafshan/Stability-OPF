function [ hx, ha,hu ] = hFunctionJacobVectorized( z )
%HFUNCTIONJACOBVECTORIZED calculates the jacobian of g(x,a,u); 
% [ hx, ha,hu ] = hFunctionJacob( z,...
 %    xq_vec, xprime_vec ) calculates the jacobian of h. 
% 
% Description of Outputs: 
% 1. hx: the jacobian of h with respect to x, size(2*N+2*G,4*G); 
% 2. ha: the jacobian of h with respect to a, size(2*N+2*G, 2*N+2*G)
% 3. hu: the jacobian of h with respect to u, size(2*N+2*G, 2*G)
% 
% Description of Inputs: 
% 1. z: vector of z=(x,a,u) combining states, algebraic and control
% variables, size(2*N+7*G,1).
% 5. xq_vec:  vector of quadrature axis synchronous reactance (pu) size(G,1)
% 6. xprime_vec: direct axis transient reactance pu, size(G,1).
%
% See also hFunctionJacob

global xq_vec xprime_vec

global N G L  gen_set load_set Gmat Bmat Cg...
    deltaIdx omegaIdx eIdx mIdx ...
    thetaIdx vIdx pgIdx qgIdx prefIdx fIdx 

% quantities:
delta=z(deltaIdx);
e=z(eIdx);
V=z(vIdx);
theta=z(thetaIdx);
Vg=V(gen_set);
thetag=theta(gen_set);


% constants
Xprime=diag(xprime_vec); 
Xq=diag(xq_vec);


% Matrix variables
E=diag(e);
thetag=theta(gen_set);
Vg=V(gen_set);
V_g=diag(Vg);
Vmat=diag(V);
cosMat=diag( cos(theta));
sinMat=diag( sin(theta));


% Indices
VgIdx=vIdx(gen_set);
thetagIdx=thetaIdx(gen_set);


h1Idx=(1:G).';
h2Idx=(h1Idx(end)+1:h1Idx(end)+G).';
h3Idx=(h2Idx(end)+1:h2Idx(end)+G).';
h4Idx=(h3Idx(end)+1:h3Idx(end)+G).';
h5Idx=(h4Idx(end)+1:h4Idx(end)+L).';
h6Idx=(h5Idx(end)+1:h5Idx(end)+L).';






hz=zeros(2*N+2*G, 2*N+8*G); 


hz(h1Idx,deltaIdx)= inv(Xprime)*E*V_g*diag(cos(delta-thetag))+...
    inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*V_g*diag(cos(2*(delta-thetag)));
hz(h1Idx,eIdx)= inv(Xprime)*V_g*diag(sin(delta-thetag));
hz(h1Idx,pgIdx)=-eye(G); 
hz(h1Idx,VgIdx)= inv(Xprime)*E*diag(sin(delta-thetag))+...
    inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*diag(sin(2*(delta-thetag)));
hz(h1Idx,thetagIdx)=-inv(Xprime)*E*V_g*diag(cos(delta-thetag))-...
    inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*V_g*diag(cos(2*(delta-thetag)));

hz(h2Idx,deltaIdx)=-inv(Xprime)*E*V_g*diag(sin(delta-thetag))...
    -inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*V_g*diag(sin(2*(delta-thetag)));
hz(h2Idx,eIdx)=inv(Xprime)*V_g*diag(cos(delta-thetag));
hz(h2Idx,qgIdx)=-eye(G); 
hz(h2Idx,VgIdx)=inv(Xprime)*E*diag(cos(delta-thetag))- ...
    inv(Xq)*inv(Xprime)*(Xprime+Xq)*V_g+...
    inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*diag(cos(2*(delta-thetag)));
hz(h2Idx,thetagIdx)=inv(Xprime)*E*V_g*diag(sin(delta-thetag))+...
    inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*V_g*diag(sin(2*(delta-thetag)));


h35=zeros(N,2*N+8*G);
h35(:, pgIdx)= -Cg;
h35(:,vIdx)= diag(cosMat*Gmat* cosMat*V) +Vmat*cosMat*Gmat*cosMat...
    - diag(cosMat*Bmat*sinMat*V)- Vmat*cosMat*Bmat*sinMat...
    +diag( sinMat*Bmat*cosMat*V)+ Vmat*sinMat*Bmat*cosMat...
    +diag(sinMat*Gmat*sinMat*V)+ Vmat*sinMat*Gmat*sinMat;

h35(:,thetaIdx)=-diag(Vmat*Gmat*Vmat*cos(theta))*sinMat- cosMat*Vmat*Gmat*Vmat*sinMat...
    +diag(Vmat*Bmat*Vmat*sin(theta))*sinMat-cosMat*Vmat*Bmat*Vmat*cosMat...
    +diag(Vmat*Bmat*Vmat*cos(theta))*cosMat-sinMat*Vmat*Bmat*Vmat*sinMat...
    +diag(Vmat*Gmat*Vmat*sin(theta))*cosMat+sinMat*Vmat*Gmat*Vmat*cosMat;


hz(h3Idx,:)= h35(gen_set,:);
hz(h5Idx,:)=h35(load_set,:); 

h46=zeros(N,2*N+8*G);
h46(:, qgIdx)= -Cg;
h46(:,vIdx)= diag(sinMat*Gmat*cosMat*V) +Vmat*sinMat*Gmat*diag(cos(theta))...
    - diag(sinMat*Bmat*sinMat*V)- Vmat*sinMat*Bmat*diag(sin(theta))...
    - diag(cosMat*Bmat*cosMat*V)- Vmat*cosMat*Bmat*diag(cos(theta))...
    -diag(cosMat*Gmat*sinMat*V)-Vmat*cosMat*Gmat*diag(sin(theta));

h46(:,thetaIdx)=diag(Vmat*Gmat*Vmat*cos(theta))*cosMat- sinMat*Vmat*Gmat*Vmat*sinMat...
    -diag(Vmat*Bmat*Vmat*sin(theta))*cosMat-sinMat*Vmat*Bmat*Vmat*cosMat...
    +diag(Vmat*Bmat*Vmat*cos(theta))*sinMat+cosMat*Vmat*Bmat*Vmat*sinMat...
    +diag(Vmat*Gmat*Vmat*sin(theta))*sinMat-cosMat*Vmat*Gmat*Vmat*cosMat;
     
hz(h4Idx,:)=h46(gen_set,:);
hz(h6Idx,:)=h46(load_set,:);


hx=hz(:, [deltaIdx;omegaIdx;eIdx;mIdx]); 
ha=hz(:,[vIdx;thetaIdx;pgIdx;qgIdx]); 
hu=hz(:,[prefIdx;fIdx]); 






end

