function out1 = sharm_C_8_0(x,y,z)
%SHARM_C_8_0
%    OUT1 = SHARM_C_8_0(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:43

t2 = y.^2;
t3 = t2.^2;
t4 = x.^2;
t5 = t4.^2;
t6 = z.^2;
t7 = t6.^2;
out1 = t3.*t5.*1.640625+t3.*t7.*2.625e1+t5.*t7.*2.625e1+t3.^2.*2.734375e-1+t5.^2.*2.734375e-1+t7.^2+t2.*t3.*t4.*1.09375-t2.*t3.*t6.*8.75+t2.*t4.*t5.*1.09375+t2.*t4.*t7.*5.25e1-t2.*t5.*t6.*2.625e1-t3.*t4.*t6.*2.625e1-t2.*t6.*t7.*1.4e1-t4.*t5.*t6.*8.75-t4.*t6.*t7.*1.4e1;
