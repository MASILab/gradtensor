function out1 = ddy_sharm_C_4_2(x,y,z)
%DDY_SHARM_C_4_2
%    OUT1 = DDY_SHARM_C_4_2(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:15

t2 = y.^2;
out1 = t2.*y.*2.23606797749979-y.*z.^2.*6.708203932499369;
