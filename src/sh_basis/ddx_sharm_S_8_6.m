function out1 = ddx_sharm_S_8_6(x,y,z)
%DDX_SHARM_S_8_6
%    OUT1 = DDX_SHARM_S_8_6(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:51

t2 = y.^2;
t3 = t2.^2;
t4 = x.^2;
t5 = t4.^2;
t6 = z.^2;
out1 = t2.*t3.*y.*-2.746090971726902+t2.*t5.*y.*3.203772800348052e1+t3.*t4.*y.*1.922263680208831e1+t3.*t6.*y.*3.844527360417662e1-t4.*t5.*y.*1.922263680208831e1+t5.*t6.*y.*1.922263680208831e2-t2.*t4.*t6.*y.*3.844527360417662e2;
