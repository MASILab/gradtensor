function new_dwi_signal = bl_correct_signal_sm(S_corpt,g,b,abvec,abval,L_mat) 
	b0 = 1;
	for v = 1:length(b)
         og = g(:,v);
         ob = b(v);

	 % adjus bvec by L*bvec and then compute the length change
          gg = L_mat * og;
          len2 = sum(gg.^2);

         % compute the new signal with length change
        new_dwi_signal(v) = b0 * exp( (log(S_corpt(v)/b0)) / len2 );
      
	end
end
