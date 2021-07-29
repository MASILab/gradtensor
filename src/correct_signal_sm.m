function new_dwi_signal = correct_signal_sm(S,b) 
	b0 = 1;
	%[S, SS, ag, ab, g, b] = sim_signal;
	L_mat = [ 1.0091 -0.0027 0.0070; -0.0031 1.0025 -0.0018 ; -0.0161 0.0041 1.0042 ];
	for i = 1:length(b)
    		%g = bvec(:,i);
    		%b = bval(i);
    		new_dwi_signal(i) = b0 *  (S(i) ^ (det(L_mat(:,:))^-2));
	end
end
