function out1 = sharm_S_5_5(x,y,z)
%SHARM_S_5_5
%    OUT1 = SHARM_S_5_5(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:24

t2 = y.^2;
t3 = x.^2;
t4 = t2.^2;
out1 = t4.*y.*7.01560760020114e-1+t3.^2.*y.*3.50780380010057-t2.*t3.*y.*7.01560760020114;
