function out1 = ddy_sharm_S_5_2(x,y,z)
%DDY_SHARM_S_5_2
%    OUT1 = DDY_SHARM_S_5_2(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:21

t2 = z.^2;
t3 = x.^2;
out1 = x.*y.^2.*z.*-1.53704261489394e1+t2.*x.*z.*1.02469507659596e1-t3.*x.*z.*5.123475382979799;
