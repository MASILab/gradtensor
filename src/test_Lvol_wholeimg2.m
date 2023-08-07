function L_pred = test_Lvol_wholeimg2(n,DT_trues,gs,bs,S_corpts, b0s)
    
    DT_trues = squeeze(DT_trues);
    S_corpts = squeeze(S_corpts);
    b0s = squeeze(b0s);
    bs = squeeze(bs);
    
    qii_stack = [];
    S_corpt_stack = [];
    b_stack = [];

    for i = 1:n
        d = DT_trues(:,i); 
        % DT_true fit is from mrtrix; the convension is D11, D22, D33, D12, D13, D23
        Dxx = d(1);
        Dxy = d(4);
        Dxz = d(5);
        Dyy = d(2);
        Dyz = d(6);
        Dzz = d(3);
        
        g = gs(:,:,i);
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
     
        qii_stack = cat(1, qii_stack, qii);
        
        b0 = b0s(i);
        S_corpt = S_corpts(:,i) / b0;
        S_corpt_stack = cat(1,S_corpt_stack,S_corpt);
        
        b = bs(:,i) ;
        b_stack = cat(1,b_stack,b);
        
        
    end
    
L0 = [  -0.000368297971855   0.002735312206693   0.002301347115275;
   0.001433235048658  -0.000134046615614   0.000836823353303;
   0.001040932606284   0.001257188013506   0.000090831513272 ];
%    1.0e-03 *
% 
%   -0.046257276345828  -0.002302050044001   0.021828850702372
%   -0.035685580669770  -0.082761974887268  -0.084670300247762
%   -0.006414668555071  -0.025708968131075   0.115919844861839

%     L0 = diag([1 1 1]);
    opts.TolFun=1e-9;
    opts.TolX = 1e-9;
%     opts.Display='iter';
    opts.PlotFcns=@optimplotfval;
    L_pred = fminsearch(@(L) myLoss2(S_corpt_stack,b_stack,qii_stack,L),L0(:),opts);
    
end