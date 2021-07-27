function sim_Dt(alpha, theta, lamda_per, lamda_parl)
	diffusivity = [lamda_per 0 0 ; 0 lamda_per 0 ; 0 0 lamda_parl];
	D = [sin(alpha)*cos(theta) sin(alpha)*sin(theta) cos(alpha)] .* diffusivity .* [sin(alpha)*cos(theta) sin(alpha)*sin(theta) cos(alpha)];
    plotDTI(D,0.002)
