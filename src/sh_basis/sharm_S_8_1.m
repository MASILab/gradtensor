function out1 = sharm_S_8_1(x,y,z)
%SHARM_S_8_1
%    OUT1 = SHARM_S_8_1(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:44

t2 = y.^2;
t3 = z.^2;
t4 = t3.^2;
t5 = t2.^2;
t6 = x.^2;
t7 = t6.^2;
out1 = t2.*t4.*y.*z.*-3.15e1-t2.*t5.*y.*z.*3.28125+t3.*t4.*y.*z.*6.0+t3.*t5.*y.*z.*2.625e1-t2.*t7.*y.*z.*9.84375+t3.*t7.*y.*z.*2.625e1-t4.*t6.*y.*z.*3.15e1-t5.*t6.*y.*z.*9.84375-t6.*t7.*y.*z.*3.28125+t2.*t3.*t6.*y.*z.*5.25e1;
