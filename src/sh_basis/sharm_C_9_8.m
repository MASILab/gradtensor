function out1 = sharm_C_9_8(x,y,z)
%SHARM_C_9_8
%    OUT1 = SHARM_C_9_8(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:54:04

t2 = x.^2;
t3 = t2.^2;
t4 = y.^2;
t5 = t4.^2;
out1 = t3.^2.*z.*2.583977731709147+t5.^2.*z.*2.583977731709147+t3.*t5.*z.*1.808784412196403e2-t2.*t3.*t4.*z.*7.235137648785613e1-t2.*t4.*t5.*z.*7.235137648785613e1;
