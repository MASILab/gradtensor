function out1 = ddx_sharm_S_8_4(x,y,z)
%DDX_SHARM_S_8_4
%    OUT1 = DDX_SHARM_S_8_4(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:49

t2 = y.^2;
t3 = t2.^2;
t4 = x.^2;
t5 = z.^2;
t6 = t4.^2;
t7 = t5.^2;
out1 = t2.*t3.*y.*-1.645305822636023-t3.*t4.*y.*4.935917467908069+t2.*t6.*y.*8.226529113180115+t3.*t5.*y.*3.948733974326455e1-t2.*t7.*y.*6.581223290544092e1+t4.*t6.*y.*1.151714075845216e1+t4.*t7.*y.*1.974366987163227e2-t5.*t6.*y.*1.974366987163227e2;
