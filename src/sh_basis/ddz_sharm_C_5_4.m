function out1 = ddz_sharm_C_5_4(x,y,z)
%DDZ_SHARM_C_5_4
%    OUT1 = DDZ_SHARM_C_5_4(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:23

t2 = x.^2;
t3 = y.^2;
out1 = t2.*t3.*-1.331117951197414e1+t2.^2.*2.218529918662356+t3.^2.*2.218529918662356;
