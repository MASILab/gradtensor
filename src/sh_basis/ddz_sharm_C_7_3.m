function out1 = ddz_sharm_C_7_3(x,y,z)
%DDZ_SHARM_C_7_3
%    OUT1 = DDZ_SHARM_C_7_3(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:37

t2 = x.^2;
t3 = z.^2;
t4 = t2.^2;
t5 = y.^2;
out1 = t5.^2.*x.*z.*5.15539765682532e1-t4.*x.*z.*1.71846588560844e1+t2.*t3.*x.*z.*4.58257569495584e1+t2.*t5.*x.*z.*3.43693177121688e1-t3.*t5.*x.*z.*1.374772708486752e2;
