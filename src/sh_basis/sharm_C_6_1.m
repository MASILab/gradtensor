function out1 = sharm_C_6_1(x,y,z)
%SHARM_C_6_1
%    OUT1 = SHARM_C_6_1(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:25

t2 = x.^2;
t3 = z.^2;
t4 = t3.^2;
t5 = t2.^2;
t6 = y.^2;
out1 = t6.^2.*x.*z.*2.8641098093474+t4.*x.*z.*4.58257569495584+t5.*x.*z.*2.8641098093474-t2.*t3.*x.*z.*1.14564392373896e1+t2.*t6.*x.*z.*5.7282196186948-t3.*t6.*x.*z.*1.14564392373896e1;
