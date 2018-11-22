function out1 = ddy_sharm_C_6_1(x,y,z)
%DDY_SHARM_C_6_1
%    OUT1 = DDY_SHARM_C_6_1(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:26

t2 = z.^2;
t3 = y.^2;
t4 = x.^2;
out1 = t2.*x.*y.*z.*-2.29128784747792e1+t3.*x.*y.*z.*1.14564392373896e1+t4.*x.*y.*z.*1.14564392373896e1;
