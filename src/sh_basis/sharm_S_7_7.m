function out1 = sharm_S_7_7(x,y,z)
%SHARM_S_7_7
%    OUT1 = SHARM_S_7_7(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:42

t2 = y.^2;
t3 = t2.^2;
t4 = x.^2;
t5 = t4.^2;
out1 = t2.*t3.*y.*-6.472598492877493e-1-t2.*t5.*y.*2.265409472507123e1+t3.*t4.*y.*1.359245683504274e1+t4.*t5.*y.*4.530818945014245;
