function [ F ] = algebraic(x,...
    delta,e,...
    pd,qd)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here




global N G L gen_set load_set Gmat Bmat Cg xq_vec xprime_vec...

algIdxV=1:N;
algIdxTheta=N+1:2*N;


V=x(algIdxV); 
theta=x(algIdxTheta); 



% constants:
Xprime=diag(xprime_vec); 
Xq=diag(xq_vec);

% vector variables:
thetag=theta(gen_set);
Vg=V(gen_set);


% matrix variables
E=diag(e);
V_g=diag(Vg);
Vmat=diag(V);
cosMat=diag( cos(theta));
sinMat=diag( sin(theta));

% h1=-pg+ inv(Xprime)*E*V_g*sin(delta-thetag)+...
%    (1/2)*inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*V_g*sin(2*(delta-thetag));
%   

% h2=-qg+inv(Xprime)*E*V_g*cos(delta-thetag)-...
%     (1/2)*inv(Xq)*inv(Xprime)*(Xprime+Xq)*V_g*Vg+ ...
%     (1/2)*inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*V_g*cos(2*(delta-thetag));



h35=Vmat*cosMat*Gmat*Vmat*cos(theta)-Vmat*cosMat*Bmat*Vmat*sin(theta)...
  +Vmat*sinMat*Bmat*Vmat*cos(theta)+Vmat*sinMat*Gmat*Vmat*sin(theta)-Cg*( inv(Xprime)*E*V_g*sin(delta-thetag)+...
   (1/2)*inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*V_g*sin(2*(delta-thetag)));

h3=h35(gen_set);
h5=h35(load_set);

h46=Vmat*sinMat*Gmat*Vmat*cos(theta)-Vmat*sinMat*Bmat*Vmat*sin(theta)...
  -Vmat*cosMat*Bmat*Vmat*cos(theta)-  Vmat*cosMat*Gmat*Vmat*sin(theta)-Cg*(inv(Xprime)*E*V_g*cos(delta-thetag)-...
    (1/2)*inv(Xq)*inv(Xprime)*(Xprime+Xq)*V_g*Vg+ ...
    (1/2)*inv(Xq)*inv(Xprime)*(Xprime-Xq)*V_g*V_g*cos(2*(delta-thetag)));
h4=h46(gen_set);
h6=h46(load_set);

F=[h3+pd(gen_set);h4+qd(gen_set);h5+pd(load_set);h6+qd(load_set)];
end

