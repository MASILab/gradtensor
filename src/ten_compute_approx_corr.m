function new_dwi_signal = ten_compute_approx_corr(S_corpt,g,b,L_mat) 
        % ten_compute_approx_corr Compute the approximate corrected signal in tensor simulation
        %
        % Inputs
        %    S_corpt          Corrput signal
        %    g                B vector
        %    b                B value
	% Outputs
	%    new_dwi_signal   Corrected signal
	b0 = 1;

	% For every volume compute the corrected signal
	for v = 1:length(b)
         og = g(:,v);
         ob = b(v);

	 % Adjus bvec by L*bvec and then compute the square of length change
         og(1) = -og(1);
         gg = L_mat * og;
         len2 = sum(gg.^2);

         % Compute the new signal by scaling the corput signal with square of length change
        new_dwi_signal(v) = b0 * exp( (log(S_corpt(v)/b0)) / len2 );

	% To compute ADC signal
        %bv_b0 = 0;
        %ADC_simple_corr = log(new_dwi_signal(v) / (b0)) * (1 / (bv_b0 - ob));
        %fprintf('ADC simple corr %f for volume %i\n', [ADC_simple_corr, v]);
      
	end
end
