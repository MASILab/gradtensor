figure;
FAdiff = FA_sim_corpt_x - FA_true_x;
FAdiff_bx = FA_corr_bx_x - FA_true_x;
FAdiff_sm = FA_corr_sm_x - FA_true_x;
FAdiff_sm_fy = FA_corr_sm_fy_x - FA_true_x;
h = boxplot([abs(FAdiff), abs(FAdiff_bx), abs(FAdiff_sm), abs(FAdiff_sm_fy)], 'Labels',{'No correction','Full correction method','Diffusivity correction method', 'Anisotropy correction method'}, 'color', colors([1 2 3 7], :),'Symbol', '.','OutlierSize', 8);
%h = violin1([FAdiff, FAdiff_bx, FAdiff_sm]);% 'Labels',{'No correction','Full correction method','Approximation method'}, 'color', colors([1 2 3], :),'Symbol', '.','OutlierSize', 12);

title('FA Comparsion for correction methods ');
ylabel('FA difference');
pretty_boxes();
