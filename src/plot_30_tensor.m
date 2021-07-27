figure;
subplot(4,4,1);
plot(1:20, 1:20);
subplot(4,4,2);
plot(1:20, 1:20);

subplot(4,4,1);
plotDTI(DT_30(:,:,[42],[38],35),0.002);
subplot(4,4,2);
plotDTI(DT_rot_lest_posB(:,:,[42],[38],35),0.002);
subplot(4,4,3);
plotDTI(DT_rot_lest_corr_bx_posB(:,:,[42],[38],35),0.002);
subplot(4,4,4);
plotDTI((DT_rot_lest_corr_bx_posB(:,:,[42],[38],35) - DT_rot_lest_posB(:,:,[42],[38],35)),0.002);
subplot(4,4,7);
plotDTI(DT_rot_lest_corr_posB(:,:,[42],[38],35),0.002);
subplot(4,4,8);
plotDTI((DT_rot_lest_corr_posB(:,:,[42],[38],35) - DT_rot_lest_posB(:,:,[42],[38],35)),0.002);

subplot(4,4,9);
plotDTI(DT_45(:,:,[42],[38],35),0.002);
subplot(4,4,10);
plotDTI(DT_rot45_lest_posB(:,:,[42],[38],35),0.002);
subplot(4,4,11);
plotDTI(DT_rot45_lest_corr_bx_posB(:,:,[42],[38],35),0.002);
subplot(4,4,12);
plotDTI((DT_rot45_lest_corr_bx_posB(:,:,[42],[38],35) - DT_rot45_lest_posB(:,:,[42],[38],35)),0.002);
subplot(4,4,15);
plotDTI(DT_rot45_lest_corr_posB(:,:,[42],[38],35),0.002);
subplot(4,4,16);
plotDTI((DT_rot45_lest_corr_posB(:,:,[42],[38],35) - DT_rot45_lest_posB(:,:,[42],[38],35)),0.002);

