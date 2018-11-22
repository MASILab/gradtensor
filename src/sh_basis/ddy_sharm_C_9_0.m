function out1 = ddy_sharm_C_9_0(x,y,z)
%DDY_SHARM_C_9_0
%    OUT1 = DDY_SHARM_C_9_0(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:54

t2 = y.^2;
t3 = z.^2;
t4 = t3.^2;
t5 = t2.^2;
t6 = x.^2;
t7 = t6.^2;
out1 = t2.*t4.*y.*z.*1.89e2+t2.*t5.*y.*z.*1.96875e1-t3.*t4.*y.*z.*3.6e1-t3.*t5.*y.*z.*1.575e2+t2.*t7.*y.*z.*5.90625e1-t3.*t7.*y.*z.*1.575e2+t4.*t6.*y.*z.*1.89e2+t5.*t6.*y.*z.*5.90625e1+t6.*t7.*y.*z.*1.96875e1-t2.*t3.*t6.*y.*z.*3.15e2;
