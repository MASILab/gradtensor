function out1 = sharm_S_5_3(x,y,z)
%SHARM_S_5_3
%    OUT1 = SHARM_S_5_3(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:53:22

t2 = y.^2;
t3 = x.^2;
t4 = t2.^2;
t5 = z.^2;
out1 = t4.*y.*5.229125165837972e-1-t3.^2.*y.*1.568737549751392-t2.*t3.*y.*1.045825033167594-t2.*t5.*y.*4.183300132670378+t3.*t5.*y.*1.254990039801113e1;
