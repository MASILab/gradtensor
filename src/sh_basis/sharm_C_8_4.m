function out1 = sharm_C_8_4(x,y,z)
%SHARM_C_8_4
%    OUT1 = SHARM_C_8_4(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:48

t2 = y.^2;
t3 = t2.^2;
t4 = x.^2;
t5 = t4.^2;
t6 = z.^2;
t7 = t6.^2;
out1 = t3.*t5.*-4.113264556590057+t3.*t7.*1.645305822636023e1+t5.*t7.*1.645305822636023e1+t3.^2.*4.113264556590057e-1+t5.^2.*4.113264556590057e-1-t2.*t3.*t4.*1.645305822636023-t2.*t3.*t6.*9.871834935816137-t2.*t4.*t5.*1.645305822636023-t2.*t4.*t7.*9.871834935816137e1+t2.*t5.*t6.*4.935917467908069e1+t3.*t4.*t6.*4.935917467908069e1-t4.*t5.*t6.*9.871834935816137;
