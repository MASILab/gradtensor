function D = sim_dt(n1_v, D0, dPerpNorm1) %trb, dirs_m)
% Usage: s_v = simDualTensors(f1, n1_v, n2_v, D0, dPerpNorm1, dPerpNorm2, trb, dirs_m)
%
% simDualTensors calculates the normalized DW signal, s_v, for each direction
% in the nDirs by 3 matrix dirs_m and diffusion weighting trb (in the same
% unists as D0).
% The signal is calculated for a voxel with two volume-averaged fibers.
% Fiber 1 has volume fraction f1 (fiber 2 has volume fraction f2 = 1 - f1). 
% Fiber 1 has its major diffusion axis parallel to the vector n1_v,
% and has perpendicular diffusivity dPerpNorm1*D0 (analogous statements
% hold for fiber 2). Diffusion is assumed to be axially symmetric in both 
% fibers.
%
% simDualTensors uses the diffusion weighting directions defined 
% in the file harDirs46.txt (mathematica output). simDualTensors is a child
% of harSim4.
%
% Last modified: 2009/01/24 (AWA)
%
% simDualFibers calls the functions:
%   op_v = perpVector(ip_v)
% Preliminaries:
%dirsFile_s = ['C:\Documents and Settings\Adam Anderson\My Documents\', ...
   % 'mathematica\harDirs46.txt'];

% For script execution:
% f1 = 0.5;
% n1_v = [1, 0, 0];
% rot = 45;       % Degrees rotation from fiber 1. 
% Note resolution is ~35 degrees for maxOrder = 6
%                    ~60 degrees for maxOrder = 4.
% n2_v = [cos(rot*pi/180), sin(rot*pi/180), 0];
% D0 = 0.8e-5;
% dPerpNorm1 = 0.3;
% dPerpNorm2 = 0.4;

% Define flags:
%tensorTestFlag = 0;         % 1 to verify tensor construction.

% Force fiber direction vectors to be normalized column vectors
% (assume they are 1D arrays):
    n1_v = n1_v / norm(n1_v);
    if (size(n1_v, 2) > 1)
        n1_v = n1_v.';
    end
    % Create second (m_v) and third (o_v) eigenvectors (these are 
    % perpendicular to n_v, but otherwise arbitrary):
    m1_v = perpVector(n1_v);    % A vector perpendicular to n1_v.
    o1_v = cross(n1_v, m1_v);   % Third orthogonal eigenvector.

    % Define true diffusion properties. Assume axisymmetric diffusion (with the
    % same trace) in both fibers:
    dPerp1 = dPerpNorm1 * D0;
    dPara1 = 3*D0 - 2*dPerp1;

    D = dPara1*(n1_v * n1_v.') + dPerp1*(m1_v * m1_v.') + ...
        dPerp1*(o1_v * o1_v.');
end
%{
% Define relative magnitudes of crossed fibers:
f2 = 1 - f1;

% Test tensor construction:
if (tensorTestFlag == 1)
    [v1_m, d1_m] = eig(trueTensor1_m);
    disp(['Fiber 1 eigenvalues are ', num2str(diag(d1_m).'/D0), ' *D0'])
    disp(['        corresponding eigenvectors are '])
    v1_m
    [v2_m, d2_m] = eig(trueTensor2_m);
    disp(['Fiber 2 eigenvalues are ', num2str(diag(d2_m).'/D0), ' *D0'])
    disp(['        corresponding eigenvectors are '])
    v2_m
end

% Now find the theta and phi coordinates of each measurement. These 
% are column vectors, rows correspond to measurements:
phi_v = atan2(dirs_m(:, 2), dirs_m(:, 1));	% Y comes first in atan2.
theta_v = acos(dirs_m(:, 3));
​
% Find normalized signal for each measurement:
nDirs = size(dirs_m, 1);
s0 = 1;
s_v = [];
for dir = 1:nDirs
    n_v = dirs_m(dir, :).';
    s_v = [s_v; s0 * f1 * exp(-trb * n_v.' * trueTensor1_m * n_v) + ...
	    s0 * f2 * exp(-trb * n_v.' * trueTensor2_m * n_v)]	    
end
%}
