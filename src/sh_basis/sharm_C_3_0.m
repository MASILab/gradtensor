function out1 = sharm_C_3_0(x,y,z)
%SHARM_C_3_0
%    OUT1 = SHARM_C_3_0(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:10

t2 = z.^2;
out1 = t2.*z-x.^2.*z.*1.5-y.^2.*z.*1.5;
