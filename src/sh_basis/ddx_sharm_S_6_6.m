function out1 = ddx_sharm_S_6_6(x,y,z)
%DDX_SHARM_S_6_6
%    OUT1 = DDX_SHARM_S_6_6(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:32

t2 = y.^2;
t3 = x.^2;
t4 = t2.^2;
out1 = t4.*y.*4.030159736288377+t3.^2.*y.*2.015079868144188e1-t2.*t3.*y.*4.030159736288377e1;
