function f = statica(K_stru,theta,rho,U,B,L,alpha,Cm)
% theta is in [deg], but should be converted in [rad] when you compute K_theta*theta
% to find Cm(theta) use the 'interp1' function
% fsolve finds theta that makes f=0, so write the stati equilibrium accordingly
Cm_i = interp1(alpha,Cm,theta);
kt = L*K_stru(3,3);

f = kt*(theta*pi/180) - 0.5*rho*(U.^2)*(B^2)*Cm_i;