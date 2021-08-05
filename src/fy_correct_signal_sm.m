function new_dwi_signal = fy_correct_signal_sm(S,g,b,L_mat) 
	b0 = 1;
	%L_mat = [ 1.0091 -0.0027 0.0070; -0.0031 1.0025 -0.0018 ; -0.0161 0.0041 1.0042 ];
	for i = 1:length(b)
        %g = g(:,i);
        %b = bval(i); 
        % scaling by lenght of adjusted bvec
        sm_abvec = L_mat * g(:,i); 
        len2 = sum(sm_abvec.^2);
        inv_l2 = 1/len2;
        
    	%new_dwi_signal(i) = b0 *  ((S(i) / b0) ^ (det(L_mat(:,:))^-2));
        %new_dwi_signal(i) = b0  *  ((S(i) / b0) ^ (len2));
        new_dwi_signal(i) = (b0^(1-inv_l2)) * (S(i)^inv_l2);
	end
end