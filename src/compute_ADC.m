function compute_ADC(D_corr_bx,g,b)
	% compute_ADC Compute ADC
        % D_corr_bx 	Corrected diffusion tensor
	% g		orginial bvec
	% b		orginial bval
        b0 = 1;
        for v = 1:length(b)
         og = g(:,v);
         ob = b(v)
	 Sbx(v) =  b0*exp(-1*ob*og'*D_corr_bx(:,:)*og);
	 bv_b0 = 0;
	 ADC_Sbx = log(Sbx(v) / (b0)) * (1 / (bv_b0 - ob));
        fprintf('ADC full corr %f for volume %i\n', [ADC_Sbx, v]);
        end
end
