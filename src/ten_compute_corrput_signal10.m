function [S_corpt, SS, abvec, abval] = ten_compute_corrput_signal(DT_mat, L_mat, g, b)
    % ten_compute_corrput_signal Compute signal with gradient non-linear fields induced
    % Inputs
    % D_mat     diffusion tensor
    % L		LR matrix
    % Output
    % S_corpt   esimated signal with LR induced
    % SS	esimated signal
    % abvec	adjected bvec (LR induced)
    % abval	adjected bval (LR induced)
    % g		original bvec
    % b		original bval

    % Original bvec and bval
    %g = [-0.0017637 0.79489 0.63422 -0.14367 -0.40313 -0.88413 -0.0020413 -0.79716 -0.6349 0.14195 0.40312 0.88295 -0.0019587 0.79384 0.63729 -0.14742 -0.39953 -0.88496 -0.0011682 -0.79556 -0.63546 0.14314 0.40294 0.88313;0.0015011 -0.40864 0.63071 -0.88244 0.79923 -0.13693 -0.015948 0.40028 -0.63819 0.87423 -0.80556 0.13019 0.0088742 -0.40719 0.63178 -0.87857 0.80395 -0.12925 -0.0081889 0.40645 -0.63243 0.87942 -0.80123 0.13217;1 0.44851 0.44717 0.44795 0.44578 0.44673 0.99987 0.45201 0.43546 0.46431 0.43425 0.45105 0.99996 0.45168 0.44127 0.4543 0.44049 0.44737 0.99997 0.44931 0.44297 0.45402 0.44235 0.45012];
    %b = [1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000];
    b0 = 1;
    abvec = zeros(size(g));
    abval = zeros(size(b));
    
    SNR = 10;
    noise_std = 1 / SNR;
    ndw_vol = length(b);
    randn('state',sum(100*clock));
    real_noise = noise_std * randn(1,ndw_vol);
    img_noise = 1i * noise_std * randn(1,ndw_vol);
    % For all volumes
    for v = 1:length(b)
        og = g(:,v);
        ob = b(v);
	
	% Estimate signal
        SS(v) =  abs((b0*exp(-1*ob*og'*DT_mat(:,:)*og)) + real_noise(v) + img_noise(v));
        bv_b0 = 0;
        %ADC_SS = log(SS(v) / (b0)) * (1 / (bv_b0 - ob));
        %fprintf('ADC no corpt %f for volume %i\n', [ADC_SS, v]);

        %here b and g need to be original uncorrected
        %bvals, and bvecs
        % Most simply, the adjusted bvec is simply L * bvec. Here we are
        % operating in the image space.
        og(1) = -og(1);
        ab = L_mat * og;

        % The bvecs were length 1 before adjustment, so now compute the length
        % change and adjust bvals accordingly. NOTE: adjust bval by the square
        % of the length, because the b value has a G^2 term but the vector
        % length is for G.
        len2 = sum(ab.^2);
        adjbval = ob .* len2;

        % Re-normalize bvecs to length 1 to compensate for the b value
        % adjustment we just made. Skip cases where b=0.
        len = sqrt(sum(ab.^2));
        lenkeeps = len~=0;
        ab(:,lenkeeps) = ab(:,lenkeeps) ./ repmat(len(lenkeeps),3,1);
        adjbvec = ab;

        adjbvec(1) = -adjbvec(1);

        % Estimate LR induced signal
        S_corpt(v) = abs((b0*exp(-1*adjbval*adjbvec'*DT_mat(:,:)*adjbvec)) + real_noise(v) + img_noise(v));
        %ADC_Scorpt = log(S_corpt(v) / (b0)) * (1 / (bv_b0 - ob));
        %fprintf('ADC corpt %f for volume %i\n', [ADC_Scorpt, v]);

        abvec(:,v) = adjbvec; 
        abval(v) = adjbval;
    end
end
	
