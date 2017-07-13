function F= genStateFunctions( x, X_QI , XPRIME_DI,  vg, thetag, pg, qg)
% genStateFunctions implicit generator state functions
% F=genStateFunctions( x, X_QI , XPRIME_DI,  theta, v, pg, qg)
% calculates the difference from zero of functions
% for generator states according to
% [CDC 2016 equations (2a), (2b)]
%  Given the load-flow solutions, the generated real 
% power output pgi  and the terminal voltage
% magnitude and phase, respectively vi and
% thetai are known for the generator at node 
% i, hence these equations provide a means to 
% compute for the internal angle deltai and
% quadrature axis internal electromotive force ei.
% x(1) = delta ; x(2) = e; 
%
% Description of outputs:
% 1. F of size(2,1), where F(1) and F(2) respectively
% show the difference from zero of equations (2a)
% and (2b).
%
% Description of inputs:
% 1. x of size(2,1) where x(1) is the generator steady-state 
% angle delta in radians,  x(2) is the generator 
% electromotive force e in pu volts.
% 2. X_QI: quadrature axis synchrnous reactance
% 3. XPRIME_DI: direct axis transient reactance
% 4. theta: generator terminal voltage angle in radians
% 5. v: generator terminal voltage magnitude in pu volts
% 6. pg: generator computed real power set point given by 
% power flow 
% 7. qg: generator computed reactive power set point given by 
% power flow
% Updates:
% 1. Hafez Bazrafshan 8/28/2016 1:05 PM
% 2. Hafez Bazrafshan 8/29/2016 2:58 PM

F(1) = -pg + (x(2).*vg./XPRIME_DI).*sin(x(1)- thetag) ...
    + ((XPRIME_DI-X_QI)./(2.*X_QI.*XPRIME_DI)).*(vg.^2).*sin(2*(x(1)-thetag)); 
F(2) = -qg + ( x(2).*vg./XPRIME_DI).*cos(x(1) - thetag)...
    -( ( (XPRIME_DI+X_QI)./(2*X_QI.*XPRIME_DI)).*(vg.^2) ...
    -  ( (XPRIME_DI- X_QI)./(2*X_QI.*XPRIME_DI)).*(vg.^2).*cos(2*(x(1)-thetag)));

end