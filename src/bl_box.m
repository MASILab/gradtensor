colors = [0.25, 0.45, 0.75; ... % Blue
          0.75, 0.00, 0.00; ... % Red
          0.25, 0.33, 0.40; ... % Dark Gray
          1.00, 0.75, 0.00; ... % Yellow
          0.90, 0.90, 0.90; ... % Light Gray
          0.29,    0, 0.57; ... % Purple
             0, 0.57, 0.57; ... % Teal
           0.1,  0.1,  0.1; ... % Black
           0.6,  0.6,  0.6];    % Gray


FAdiff_sm = FA_corr_sm_fy_x - FA_true_x;

FAdiff = FA_sim_corpt_x - FA_true_x;

FAdiff_bx = FA_corr_bx_x - FA_true_x;

figure; 
boxplot([abs(FAdiff), abs(FAdiff_sm), abs(FAdiff_bx)],  'Labels',{'No correction','BL method','Full correction method'}, 'color', colors([1 2 3], :),'Symbol', '.','OutlierSize', 8); 
title('Abs FA Difference');
pause(1);

MDdiff_sm = MD_corr_sm_fy_x - MD_true_x;

MDdiff = MD_sim_corpt_x - MD_true_x;

MDdiff_bx = MD_corr_bx_x - MD_true_x;

figure; 
boxplot([abs(MDdiff), abs(MDdiff_sm), abs(MDdiff_bx)] ,  'Labels',{'No correction','BL method','Full correction method'}, 'color', colors([1 2 3], :),'Symbol', '.','OutlierSize', 8); 
title('Abs MD Difference');

