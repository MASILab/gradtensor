function out1 = ddx_sharm_S_9_4(x,y,z)
%DDX_SHARM_S_9_4
%    OUT1 = DDX_SHARM_S_9_4(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:54:00

t2 = y.^2;
t3 = z.^2;
t4 = t3.^2;
t5 = t2.^2;
t6 = x.^2;
t7 = t6.^2;
out1 = t2.*t4.*y.*z.*-1.06119036934944e2-t2.*t5.*y.*z.*1.326487961686799e1+t3.*t5.*y.*z.*1.06119036934944e2+t2.*t7.*y.*z.*6.632439808433997e1-t3.*t7.*y.*z.*5.305951846747198e2+t4.*t6.*y.*z.*3.183571108048319e2-t5.*t6.*y.*z.*3.979463885060398e1+t6.*t7.*y.*z.*9.285415731807596e1;
