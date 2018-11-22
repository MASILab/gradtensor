function out1 = ddy_sharm_S_8_4(x,y,z)
%DDY_SHARM_S_8_4
%    OUT1 = DDY_SHARM_S_8_4(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:49

t2 = x.^2;
t3 = y.^2;
t4 = t2.^2;
t5 = z.^2;
t6 = t3.^2;
t7 = t5.^2;
out1 = t2.*t4.*x.*1.645305822636023+t3.*t4.*x.*4.935917467908069-t2.*t6.*x.*8.226529113180115+t2.*t7.*x.*6.581223290544092e1-t3.*t6.*x.*1.151714075845216e1-t4.*t5.*x.*3.948733974326455e1-t3.*t7.*x.*1.974366987163227e2+t5.*t6.*x.*1.974366987163227e2;
