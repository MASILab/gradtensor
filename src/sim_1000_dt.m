
MD = 0.00008;
% Trace fixed
%TrD = 0.0021; %derek jones
%MD = TrD/3;

% FA change
FA = 0.9;
%l1 = (MD) * (1 + (2*FA) / sqrt( 3 - (2 * FA^2)));
%lprep = (MD) * (1 - FA / sqrt( 3 - (2 * FA^2)));
lprep = MD * (1 - FA/ (3-2*FA^2)^(1/2));

% for PVE
% n = 100;
% for i = 1:n
     PEV = [1, 0, 0];
%     
%     while norm(PEV) < 0.001
%         x = randn();
%         y = randn();
%         z = randn();
%         PEV = [x y z];
%     end
     PEV = PEV / norm(PEV);
     disp(PEV)
    D = sim_dt(PEV, MD,  lprep);
    FA_new = compute_FA(D);
    plotDTI(D,0.002);
% end



