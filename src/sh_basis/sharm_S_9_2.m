function out1 = sharm_S_9_2(x,y,z)
%SHARM_S_9_2
%    OUT1 = SHARM_S_9_2(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:57

t2 = z.^2;
t3 = t2.^2;
t4 = y.^2;
t5 = t4.^2;
t6 = x.^2;
t7 = t6.^2;
out1 = t2.*t3.*x.*y.*z.*3.146426544510455e1+t2.*t5.*x.*y.*z.*6.88280806611662e1-t3.*t4.*x.*y.*z.*1.101249290578659e2+t2.*t7.*x.*y.*z.*6.88280806611662e1-t3.*t6.*x.*y.*z.*1.101249290578659e2-t4.*t5.*x.*y.*z.*6.88280806611662-t4.*t7.*x.*y.*z.*2.064842419834986e1-t5.*t6.*x.*y.*z.*2.064842419834986e1-t6.*t7.*x.*y.*z.*6.88280806611662+t2.*t4.*t6.*x.*y.*z.*1.376561613223324e2;
