function out1 = ddz_sharm_S_9_8(x,y,z)
%DDZ_SHARM_S_9_8
%    OUT1 = DDZ_SHARM_S_9_8(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:54:05

t2 = x.^2;
t3 = y.^2;
t4 = t3.^2;
t5 = t2.^2;
out1 = t2.*t4.*x.*y.*1.447027529757123e2+t2.*t5.*x.*y.*2.067182185367318e1-t3.*t4.*x.*y.*2.067182185367318e1-t3.*t5.*x.*y.*1.447027529757123e2;
