function out1 = sharm_C_8_6(x,y,z)
%SHARM_C_8_6
%    OUT1 = SHARM_C_8_6(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:50

t2 = y.^2;
t3 = t2.^2;
t4 = x.^2;
t5 = t4.^2;
t6 = z.^2;
out1 = t3.^2.*4.576818286211503e-1-t5.^2.*4.576818286211503e-1-t2.*t3.*t4.*6.407545600696104-t2.*t3.*t6.*6.407545600696104+t2.*t4.*t5.*6.407545600696104-t2.*t5.*t6.*9.611318401044156e1+t3.*t4.*t6.*9.611318401044156e1+t4.*t5.*t6.*6.407545600696104;
