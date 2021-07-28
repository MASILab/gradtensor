figure;
axis equal;
subplot(1,4,1);
D = sim_dt(0,0,0);
plotDTI(D,0.002);
subplot(1,4,2);
D = sim_dt(30,30,30);
plotDTI(D,0.002);
subplot(1,4,3)
D = sim_dt(60,60,60);
plotDTI(D,0.002);
subplot(1,4,4)
D = sim_dt(90,90,90);
plotDTI(D,0.002);
