function out1 = ddy_sharm_S_9_5(x,y,z)
%DDY_SHARM_S_9_5
%    OUT1 = DDY_SHARM_S_9_5(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:54:01

t2 = y.^2;
t3 = t2.^2;
t4 = x.^2;
t5 = t4.^2;
t6 = z.^2;
t7 = t6.^2;
out1 = t3.*t5.*-2.774548633055023e1+t3.*t7.*1.109819453222009e2+t5.*t7.*1.109819453222009e2+t3.^2.*3.567276813927887+t5.^2.*1.98182045218216-t2.*t3.*t4.*2.219638906444019e1-t2.*t3.*t6.*7.768736172554066e1-t2.*t4.*t7.*6.658916719332056e2+t2.*t5.*t6.*1.664729179833014e2+t3.*t4.*t6.*4.994187539499042e2-t4.*t5.*t6.*5.549097266110047e1;
