function out1 = sharm_S_5_4(x,y,z)
%SHARM_S_5_4
%    OUT1 = SHARM_S_5_4(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:23

t2 = y.^2;
t3 = x.^2;
out1 = t2.*x.*y.*z.*-8.874119674649425+t3.*x.*y.*z.*8.874119674649425;
