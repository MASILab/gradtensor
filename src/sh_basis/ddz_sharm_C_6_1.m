function out1 = ddz_sharm_C_6_1(x,y,z)
%DDZ_SHARM_C_6_1
%    OUT1 = DDZ_SHARM_C_6_1(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:26

t2 = x.^2;
t3 = y.^2;
t4 = z.^2;
t5 = t2.^2;
out1 = t5.*x.*2.8641098093474+t3.^2.*x.*2.8641098093474+t4.^2.*x.*2.29128784747792e1+t2.*t3.*x.*5.7282196186948-t2.*t4.*x.*3.43693177121688e1-t3.*t4.*x.*3.43693177121688e1;
