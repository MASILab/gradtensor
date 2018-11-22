function out1 = ddx_sharm_S_6_4(x,y,z)
%DDX_SHARM_S_6_4
%    OUT1 = DDX_SHARM_S_6_4(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:30

t2 = y.^2;
t3 = x.^2;
t4 = t2.^2;
t5 = z.^2;
out1 = t4.*y.*1.984313483298443-t3.^2.*y.*9.921567416492214-t2.*t5.*y.*1.984313483298443e1+t3.*t5.*y.*5.952940449895329e1;
