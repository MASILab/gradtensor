figure;
MDdiff = MD_sim_corpt_x - MD_true_x;
MDdiff_bx = MD_corr_bx_x - MD_true_x;
MDdiff_sm = MD_corr_sm_x - MD_true_x;
MDdiff_sm_fy = MD_corr_sm_fy_x - MD_true_x;
h = boxplot([MDdiff, MDdiff_bx, MDdiff_sm, MDdiff_sm_fy], 'Labels',{'No correction','Full correction method','Diffusivity correction method', 'Anisotropy correction method'}, 'color', colors([1 2 3 7], :),'Symbol', '.','OutlierSize', 8);
%h = violin1([FAdiff, FAdiff_bx, FAdiff_sm]);% 'Labels',{'No correction','Full correction method','Approximation method'}, 'color', colors([1 2 3], :),'Symbol', '.','OutlierSize', 12);

title('MD Comparsion for correction methods ');
ylabel('MD difference');
pretty_boxes();
