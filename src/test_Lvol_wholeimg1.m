function L_pred = test_Lvol_wholeimg1(DT_true,g,b,S_corpt, b0)
    DT_true = squeeze(DT_true);
    d = DT_true;
    % DT_true fit is from mrtrix; the convension is D11, D22, D33, D12, D13, D23
    Dxx = d(1);
    Dxy = d(4);
    Dxz = d(5);
    Dyy = d(2);
    Dyz = d(6);
    Dzz = d(3);


    gx = g(1,:)';
    gy = g(2,:)';
    gz = g(3,:)';

    qii = [Dxx.*gx.^2  2.*Dxx.*gx.*gy  2.*Dxx.*gx.*gz  2.*Dxy.*gx.^2  2.*Dxy.*gx.*gy  2.*Dxy.*gx.*gz  2.*Dxz.*gx.^2  ...
        2.*Dxz.*gx.*gy  2.*Dxz.*gx.*gz  Dxx.*gy.^2  2.*Dxx.*gy.*gz  2.*Dxy.*gx.*gy  2.*Dxy.*gy.^2  2.*Dxy.*gy.*gz  ...
        2.*Dxz.*gx.*gy  2.*Dxz.*gy.^2  2.*Dxz.*gy.*gz  Dxx.*gz.^2  2.*Dxy.*gx.*gz  2.*Dxy.*gy.*gz  2.*Dxy.*gz.^2  ...
        2.*Dxz.*gx.*gz  2.*Dxz.*gy.*gz  2.*Dxz.*gz.^2  Dyy.*gx.^2  2.*Dyy.*gx.*gy  2.*Dyy.*gx.*gz  2.*Dyz.*gx.^2  ...
        2.*Dyz.*gx.*gy  2.*Dyz.*gx.*gz  Dyy.*gy.^2  2.*Dyy.*gy.*gz  2.*Dyz.*gx.*gy  2.*Dyz.*gy.^2  2.*Dyz.*gy.*gz  ...
        Dyy.*gz.^2  2.*Dyz.*gx.*gz  2.*Dyz.*gy.*gz  2.*Dyz.*gz.^2  Dzz.*gx.^2  2.*Dzz.*gx.*gy  2.*Dzz.*gx.*gz  ...
        Dzz.*gy.^2  2.*Dzz.*gy.*gz  Dzz.*gz.^2];

    S_corpt = squeeze(S_corpt);
    S_corpt_b0 = S_corpt / b0;
    

    
    L0 = diag([1 1 1]);
    opts.TolFun=1e-9;
    opts.TolX = 1e-9;
%     opts.Display='iter';
    opts.PlotFcns=@optimplotfval;
    L_pred = fminsearch(@(L) myLoss(S_corpt_b0,b,qii,L),L0(:),opts);
    disp(L_pred)
    