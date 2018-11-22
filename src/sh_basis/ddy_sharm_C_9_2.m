function out1 = ddy_sharm_C_9_2(x,y,z)
%DDY_SHARM_C_9_2
%    OUT1 = DDY_SHARM_C_9_2(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:56

t2 = y.^2;
t3 = z.^2;
t4 = t3.^2;
t5 = t2.^2;
t6 = x.^2;
t7 = t6.^2;
out1 = t2.*t4.*y.*z.*2.202498581157318e2+t2.*t5.*y.*z.*2.753123226446648e1-t3.*t4.*y.*z.*3.146426544510455e1-t3.*t5.*y.*z.*2.064842419834986e2+t3.*t7.*y.*z.*6.88280806611662e1+t5.*t6.*y.*z.*4.129684839669972e1-t6.*t7.*y.*z.*1.376561613223324e1-t2.*t3.*t6.*y.*z.*1.376561613223324e2;
