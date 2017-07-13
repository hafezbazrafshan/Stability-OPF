function [ hx, ha,hu ] = hFunctionJacob( z,...
     xq_vec, xprime_vec )
%HFUNCTIONJACOB calculates the jacobian of g(x,a,u); 
% [ hx, ha,hu ] = hFunctionJacob( z,...
 %    xq_vec, xprime_vec ) calculates the jacobian of h. 
% 
% Description of Outputs: 
% 1. hx: the jacobian of h with respect to x, size(2*N+2*G,3*G); 
% 2. ha: the jacobian of h with respect to a, size(2*N+2*G, 2*N+2*G)
% 3. hu: the jacobian of h with respect to u, size(2*N+2*G, 2*G)
% 
% Description of Inputs: 
% 1. z: vector of z=(x,a,u) combining states, algebraic and control
% variables, size(2*N+7*G,1).
% 5. xq_vec:  vector of quadrature axis synchronous reactance (pu) size(G,1)
% 6. xprime_vec: direct axis transient reactance pu, size(G,1).
%
% See also hFunctionJacobVectorized

global N G L  gen_set load_set Gmat Bmat ...
    deltaIdx omegaIdx eIdx ...
    thetaIdx vIdx pgIdx qgIdx mIdx fIdx 

h1Idx=1:G;
h2Idx=h1Idx(end)+1:h1Idx(end)+G;
h3Idx=h2Idx(end)+1:h2Idx(end)+G;
h4Idx=h3Idx(end)+1:h3Idx(end)+G;
h5Idx=h4Idx(end)+1:h4Idx(end)+L;
h6Idx=h5Idx(end)+1:h5Idx(end)+L;

V=z(vIdx);
theta=z(thetaIdx);
Vg=V(gen_set);
thetag=theta(gen_set);
VgIdx=vIdx(gen_set);
thetagIdx=thetaIdx(gen_set);
delta=z(deltaIdx);
e=z(eIdx);

hz=zeros(2*N+2*G, 2*N+7*G); 
hznRows=2*N+2*G; 
hznCols=2*N+7*G;

hz(sub2ind([hznRows hznCols], h1Idx, deltaIdx)) =  (1./xprime_vec).*(e.*Vg.*cos(delta-thetag))...
    + (1./xq_vec).*(1./xprime_vec).*(xprime_vec-xq_vec).*Vg.*Vg.*cos(2*(delta-thetag));
hz(sub2ind([hznRows hznCols], h1Idx,eIdx)) = (1./xprime_vec).*Vg.*sin(delta-thetag);
hz(sub2ind([hznRows hznCols], h1Idx,pgIdx))= -1; 
hz(sub2ind([hznRows hznCols], h1Idx, VgIdx)) = (1./xprime_vec).*e.*sin(delta-thetag) +...
    (1./xq_vec).*(1./xprime_vec).*(xprime_vec-xq_vec).*Vg.*sin(2*(delta-thetag));
hz(sub2ind([hznRows hznCols], h1Idx, thetagIdx)) = -(1./xprime_vec).*e.*Vg.*cos(delta-thetag)...
    -(1./xq_vec).*(1./xprime_vec).* (xprime_vec-xq_vec).*Vg.*Vg.*cos(2*(delta-thetag));





hz(sub2ind([hznRows hznCols], h2Idx, deltaIdx)) = - (1./xprime_vec).*(e.*Vg.*sin(delta-thetag))...
    - (1./xq_vec).*(1./xprime_vec).*(xprime_vec-xq_vec).*Vg.*Vg.*sin(2*(delta-thetag));
hz(sub2ind([hznRows hznCols], h2Idx,eIdx)) = (1./xprime_vec).*Vg.*cos(delta-thetag);
hz(sub2ind([hznRows hznCols], h2Idx,qgIdx))= -1; 
hz(sub2ind([hznRows hznCols], h2Idx, VgIdx)) = (1./xprime_vec).*e.*cos(delta-thetag) ...
    -  (1./xq_vec).*(1./xprime_vec).*(xprime_vec+xq_vec).*Vg+...
    (1./xq_vec).*(1./xprime_vec).*(xprime_vec-xq_vec).*Vg.*cos(2*(delta-thetag));
hz(sub2ind([hznRows hznCols], h2Idx, thetagIdx)) = (1./xprime_vec).*e.*Vg.*sin(delta-thetag)...
    +(1./xq_vec).*(1./xprime_vec).* (xprime_vec-xq_vec).*Vg.*Vg.*sin(2*(delta-thetag));


h3v=zeros(G,N); % only as big as V

for ii=1:G
    mIndex=gen_set(ii);
    for jj=1:N
        if mIndex==jj
            h3v(ii,jj)= Gmat(mIndex,:)*diag(V)*cos(theta(mIndex)-theta)+...
                Bmat(mIndex,:)*diag(V)*sin(theta(mIndex)-theta)+...
                V(mIndex)*Gmat(mIndex,mIndex);
        else
            h3v(ii, jj)= V(mIndex)*(  Gmat(mIndex,jj)*cos(theta(mIndex)-theta(jj))+...
            Bmat(mIndex,jj)*sin(theta(mIndex)-theta(jj)));
        
        end
    end
end


h3theta=zeros(G,N); % only as big as theta

for ii=1:G
    mIndex=gen_set(ii);
    for jj=1:N
        if mIndex==jj
            
            h3theta(ii,jj)=V(mIndex)* (-Gmat(mIndex,:)*diag(V)*sin(theta(mIndex)-theta)+...
                Bmat(mIndex,:)*diag(V)*cos(theta(mIndex)-theta))-Bmat(mIndex,mIndex)*V(mIndex)^2;
                
        
        else
            h3theta(ii, jj)= V(mIndex)*(  Gmat(mIndex,jj)*V(jj)*sin(theta(mIndex)-theta(jj))-...
            Bmat(mIndex,jj)*V(jj)*cos(theta(mIndex)-theta(jj)));
        
        end
    end
end

hz(h3Idx,vIdx)=h3v;
hz(h3Idx,thetaIdx)=h3theta;
hz(sub2ind([hznRows hznCols], h3Idx,pgIdx))=-1;







h4v=zeros(G,N); % only as big as V

for ii=1:G
    mIndex=gen_set(ii);
    for jj=1:N
        if mIndex==jj
            
            h4v(ii,jj)= -Bmat(mIndex,:)*diag(V)*cos(theta(mIndex)-theta)+...
                Gmat(mIndex,:)*diag(V)*sin(theta(mIndex)-theta)...
                -V(mIndex)*Bmat(mIndex,mIndex);
        
        else
            h4v(ii, jj)= V(mIndex)*( - Bmat(mIndex,jj)*cos(theta(mIndex)-theta(jj))+...
            Gmat(mIndex,jj)*sin(theta(mIndex)-theta(jj)));
        
        end
    end
end


h4theta=zeros(G,N); % only as big as theta

for ii=1:G
    mIndex=gen_set(ii);
    for jj=1:N
        if mIndex==jj
            
            h4theta(ii,jj)=V(mIndex)* (Bmat(mIndex,:)*diag(V)*sin(theta(mIndex)-theta)+...
                Gmat(mIndex,:)*diag(V)*cos(theta(mIndex)-theta))-Gmat(mIndex,mIndex)*V(mIndex)^2;
                
        
        else
            h4theta(ii, jj)= V(mIndex)*( - Bmat(mIndex,jj)*V(jj)*sin(theta(mIndex)-theta(jj))-...
            Gmat(mIndex,jj)*V(jj)*cos(theta(mIndex)-theta(jj)));
        
        end
    end
end

hz(h4Idx,vIdx)=h4v;
hz(h4Idx,thetaIdx)=h4theta;
hz(sub2ind([hznRows hznCols], h4Idx,qgIdx))=-1;




%% h5
h5v=zeros(L,N); % only as big as V

for ii=1:L
    mIndex=load_set(ii);
    for jj=1:N
        if mIndex==jj
            
            h5v(ii,jj)= Gmat(mIndex,:)*diag(V)*cos(theta(mIndex)-theta)+...
                Bmat(mIndex,:)*diag(V)*sin(theta(mIndex)-theta)+...
                V(mIndex)*Gmat(mIndex,mIndex);
        
        else
            h5v(ii, jj)= V(mIndex)*(  Gmat(mIndex,jj)*cos(theta(mIndex)-theta(jj))+...
            Bmat(mIndex,jj)*sin(theta(mIndex)-theta(jj)));
        
        end
    end
end


h5theta=zeros(L,N); % only as big as theta

for ii=1:L
    mIndex=load_set(ii);
    for jj=1:N
        if mIndex==jj
            
            h5theta(ii,jj)=V(mIndex)* (-Gmat(mIndex,:)*diag(V)*sin(theta(mIndex)-theta)+...
                Bmat(mIndex,:)*diag(V)*cos(theta(mIndex)-theta))-Bmat(mIndex,mIndex)*V(mIndex)^2;
                
        
        else
            h5theta(ii, jj)= V(mIndex)*(  Gmat(mIndex,jj)*V(jj)*sin(theta(mIndex)-theta(jj))-...
            Bmat(mIndex,jj)*V(jj)*cos(theta(mIndex)-theta(jj)));
        
        end
    end
end

hz(h5Idx,vIdx)=h5v;
hz(h5Idx,thetaIdx)=h5theta;



%% 

h6v=zeros(L,N); % only as big as V

for ii=1:L
    mIndex=load_set(ii);
    for jj=1:N
        if mIndex==jj
            
            h6v(ii,jj)= -Bmat(mIndex,:)*diag(V)*cos(theta(mIndex)-theta)+...
                Gmat(mIndex,:)*diag(V)*sin(theta(mIndex)-theta)...
                -V(mIndex)*Bmat(mIndex,mIndex);
        
        else
            h6v(ii, jj)= V(mIndex)*( - Bmat(mIndex,jj)*cos(theta(mIndex)-theta(jj))+...
            Gmat(mIndex,jj)*sin(theta(mIndex)-theta(jj)));
        
        end
    end
end


h6theta=zeros(L,N); % only as big as theta

for ii=1:L
    mIndex=load_set(ii);
    for jj=1:N
        if mIndex==jj
            
            h6theta(ii,jj)=V(mIndex)* (Bmat(mIndex,:)*diag(V)*sin(theta(mIndex)-theta)+...
                Gmat(mIndex,:)*diag(V)*cos(theta(mIndex)-theta))-Gmat(mIndex,mIndex)*V(mIndex)^2;
                
        
        else
            h6theta(ii, jj)= V(mIndex)*(  -Bmat(mIndex,jj)*V(jj)*sin(theta(mIndex)-theta(jj))-...
            Gmat(mIndex,jj)*V(jj)*cos(theta(mIndex)-theta(jj)));
        
        end
    end
end

hz(h6Idx,vIdx)=h6v;
hz(h6Idx,thetaIdx)=h6theta;


hx=hz(:, [deltaIdx,omegaIdx,eIdx]); 
ha=hz(:,[vIdx,thetaIdx,pgIdx,qgIdx]); 
hu=hz(:,[mIdx,fIdx]); 

end

