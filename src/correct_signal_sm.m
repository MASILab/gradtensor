function new_dwi_signal = correct_signal_sm(S,g,b,L_mat) 
	b0 = 1;
	%L_mat = [ 1.0091 -0.0027 0.0070; -0.0031 1.0025 -0.0018 ; -0.0161 0.0041 1.0042 ];
	for i = 1:length(b)
    		%g = bvec(:,i);
    		%b = bval(i);
            inv_l2 = 1 / det(L_mat(:,:)) ^ 2 ;
    		%new_dwi_signal(i) = b0 *  ((S(i) / b0) ^ (det(L_mat(:,:))^-2));
            new_dwi_signal(i) = (b0^(1-inv_l2)) * (S(i)^inv_l2);
                                %b0_ovl *((S(i)/b0 )^ (det(L_mat(:,:))^-2))
                                %b0_vol(x,y,z)*((dwmri_vols(x,y,z,i)/b0_vol(x,y,z))^ (det(L_mat(:,:))^-2))
	end
end
