function [abvec, abval] = ten_compute_emp_noise_LR(g, b)
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

    b0 = 1;
    abvec = zeros(size(g));
    abval = zeros(size(b));
    SNR = 100;
    noise_std = 0.3333 / SNR;
    ndw_vol = length(b);
    randn('state',sum(100*clock));
    real_noise = noise_std * randn(1,ndw_vol);
    img_noise = 1i * noise_std * randn(1,ndw_vol);
    L_mat = [1 0 0 ; 0 1 0 ; 0 0 1];
    % For all volumes
    for v = 1:length(b)
        og = g(:,v);
        ob = b(v);
	
        % Estimate signal
        L_noise =  abs(L_mat + real_noise(v) + img_noise(v));

        % Most simply, the adjusted bvec is simply L * bvec. Here we are
        % operating in the image space.
        og(1) = -og(1);
        ab = L_noise * og;

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

        abvec(:,v) = adjbvec; 
        abval(v) = adjbval;

    end
end
	
