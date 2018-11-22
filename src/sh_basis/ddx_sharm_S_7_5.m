function out1 = ddx_sharm_S_7_5(x,y,z)
%DDX_SHARM_S_7_5
%    OUT1 = DDX_SHARM_S_7_5(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:40

t2 = x.^2;
t3 = y.^2;
t4 = t3.^2;
t5 = t2.^2;
t6 = z.^2;
out1 = t4.*x.*y.*8.549259836383499-t5.*x.*y.*1.42487663939725e1+t2.*t3.*x.*y.*9.499177595981665+t2.*t6.*x.*y.*1.1399013115178e2-t3.*t6.*x.*y.*1.1399013115178e2;
