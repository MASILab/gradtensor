function out1 = ddz_sharm_C_9_0(x,y,z)
%DDZ_SHARM_C_9_0
%    OUT1 = DDZ_SHARM_C_9_0(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:54

t2 = y.^2;
t3 = t2.^2;
t4 = x.^2;
t5 = t4.^2;
t6 = z.^2;
t7 = t6.^2;
out1 = t3.*t5.*1.4765625e1+t3.*t7.*2.3625e2+t5.*t7.*2.3625e2+t3.^2.*2.4609375+t5.^2.*2.4609375+t7.^2.*9.0+t2.*t3.*t4.*9.84375-t2.*t3.*t6.*7.875e1+t2.*t4.*t5.*9.84375+t2.*t4.*t7.*4.725e2-t2.*t5.*t6.*2.3625e2-t3.*t4.*t6.*2.3625e2-t2.*t6.*t7.*1.26e2-t4.*t5.*t6.*7.875e1-t4.*t6.*t7.*1.26e2;
