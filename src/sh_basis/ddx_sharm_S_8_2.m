function out1 = ddx_sharm_S_8_2(x,y,z)
%DDX_SHARM_S_8_2
%    OUT1 = DDX_SHARM_S_8_2(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:46

t2 = y.^2;
t3 = t2.^2;
t4 = x.^2;
t5 = z.^2;
t6 = t4.^2;
t7 = t5.^2;
out1 = t2.*t3.*y.*-7.843687748756958e-1-t3.*t4.*y.*7.059318973881262-t2.*t6.*y.*1.176553162313544e1+t3.*t5.*y.*2.353106324627087e1-t2.*t7.*y.*6.274950199005567e1-t4.*t6.*y.*5.490581424129871-t4.*t7.*y.*1.88248505970167e2+t5.*t6.*y.*1.176553162313544e2+t5.*t7.*y.*2.509980079602227e1+t2.*t4.*t5.*y.*1.411863794776252e2;
