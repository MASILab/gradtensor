function out1 = ddx_sharm_S_9_3(x,y,z)
%DDX_SHARM_S_9_3
%    OUT1 = DDX_SHARM_S_9_3(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:58

t2 = x.^2;
t3 = y.^2;
t4 = t3.^2;
t5 = t2.^2;
t6 = z.^2;
t7 = t6.^2;
out1 = t2.*t4.*x.*y.*-9.011711130523436-t2.*t5.*x.*y.*9.011711130523436-t3.*t5.*x.*y.*1.802342226104687e1-t2.*t7.*x.*y.*5.407026678314062e2-t3.*t7.*x.*y.*1.802342226104687e2+t4.*t6.*x.*y.*2.703513339157031e1+t5.*t6.*x.*y.*2.433162005241328e2+t6.*t7.*x.*y.*1.44187378088375e2+t2.*t3.*t6.*x.*y.*2.703513339157031e2;
