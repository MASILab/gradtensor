function out1 = sharm_C_9_5(x,y,z)
%SHARM_C_9_5
%    OUT1 = SHARM_C_9_5(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Jan-2017 14:54:00

t2 = x.^2;
t3 = t2.^2;
t4 = y.^2;
t5 = z.^2;
t6 = t4.^2;
t7 = t3.^2;
t8 = t5.^2;
out1 = t7.*x.*3.963640904364319e-1+t6.^2.*x.*1.98182045218216-t3.*t6.*x.*5.549097266110047+t3.*t8.*x.*2.219638906444019e1+t6.*t8.*x.*1.109819453222009e2-t2.*t3.*t4.*x.*3.170912723491455-t2.*t3.*t5.*x.*1.109819453222009e1+t3.*t4.*t5.*x.*9.988375078998085e1+t2.*t5.*t6.*x.*5.549097266110047e1-t2.*t4.*t8.*x.*2.219638906444019e2-t4.*t5.*t6.*x.*5.549097266110047e1;
