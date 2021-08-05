figure;
theta_allCopy(theta_allCopy < 18) = 1;
theta_allCopy(theta_allCopy > 18*1 & theta_allCopy <= 18*2)  = 2;
theta_allCopy(theta_allCopy > 18*2 & theta_allCopy <= 18*3)  = 3;
theta_allCopy(theta_allCopy > 18*3 & theta_allCopy <= 18*4)  = 4;
theta_allCopy(theta_allCopy > 18*4 & theta_allCopy <= 18*5)  = 5;
theta_allCopy(theta_allCopy > 18*5 & theta_allCopy <= 18*6)  = 6;
theta_allCopy(theta_allCopy > 18*6 & theta_allCopy <= 18*7)  = 7;
theta_allCopy(theta_allCopy > 18*7 & theta_allCopy <= 18*8)  = 8;
theta_allCopy(theta_allCopy > 18*8 & theta_allCopy <= 18*9)  = 9;
theta_allCopy(theta_allCopy > 18*9 & theta_allCopy <= 18*10)  = 10;

boxplot(FA_true_x,baba,theta_allCopy);
title('theta affect on FA');xlabel('FA true');ylabel('FA diff');
legend('0-18theta','18-36theta','36-54theta','54-72theta','72-90theta','90-108theta','108-126theta','126-144theta','144-162theta','162-180theta');

figure;
mo = 0.2;
FA_allCopy(FA_allCopy <= 0.2) = 0;
FA_allCopy(FA_allCopy > mo*1 & FA_allCopy <= mo*2)  = 2;
FA_allCopy(FA_allCopy > mo*2 & FA_allCopy <= mo*3)  = 3;
FA_allCopy(FA_allCopy > mo*3 & FA_allCopy <= mo*4)  = 4;
FA_allCopy(FA_allCopy > mo*4 & FA_allCopy <= mo*5)  = 5;
cmap = [1 0 0; 0 0 0 ; 0 1 0 ; 0 1 1 ; 1 0 1 ];
gscatter(FA_sim_corpt_x,baba,FA_allCopy,cmap, '');
title('FA affect on FA');xlabel('FA corpt');ylabel('FA diff');
legend('0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1');

figure;
Ldet_allCopy(Ldet_allCopy < 0.8) = 0;
Ldet_allCopy(Ldet_allCopy > 0.8 & Ldet_allCopy <= 0.9)  = 2;
Ldet_allCopy(Ldet_allCopy > 0.9 & Ldet_allCopy <= 1.0524)  = 3;
%Ldet_allCopy(Ldet_allCopy > mo*3 & Ldet_allCopy <= mo*4)  = 4;
%Ldet_allCopy(Ldet_allCopy > mo*4 & FA_allCopy <= 1.0524)  = 5;
gscatter(FA_sim_corpt_x,baba,Ldet_allCopy, cmap);
title('det L affect on FA');xlabel('FA corpt');ylabel('FA diff');
%legend('0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.9-1');
legend('0-0.8','0.8-0.9','0.9-1')%,'0.6-0.8','0.8-1')










