function out1 = sharm_C_7_7(x,y,z)
%SHARM_C_7_7
%    OUT1 = SHARM_C_7_7(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:42

t2 = x.^2;
t3 = y.^2;
t4 = t2.^2;
t5 = t3.^2;
out1 = t2.*t4.*x.*6.472598492877493e-1+t2.*t5.*x.*2.265409472507123e1-t3.*t4.*x.*1.359245683504274e1-t3.*t5.*x.*4.530818945014245;
