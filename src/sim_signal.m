function [S, SS, ag, ab, g, b] = sim_signal()
    DT_mat = sim_dt(0,0,0);
    disp(DT_mat)
    g = [-0.000737 0.79383 0.63632 -0.14534 -0.39912 -0.88473; 0.00026791 -0.41093 0.62867 -0.88214 0.80005 -0.13373; 1 0.44829 0.44707 0.448 0.44791 0.44652];
    b = [1000 1000 1000 1000 1000 1000];
    L_mat = [ 1.0091 -0.0027 0.0070; -0.0031 1.0025 -0.0018 ; -0.0161 0.0041 1.0042 ];
    %L_mat = [1 0 0 ; 0 1 0 ; 0 0 1];
    %L_mat = [ 0.8091 -0.0027 0.1070; 0.0031 -1.3025 -0.1018 ; -0.0161 -0.0041 -1.4042 ];
    b0 = 1;

    ag = zeros(size(g));
    ab = zeros(size(b));
    for v = 1:length(b)
        og = g(:,v);
        ob = b(v);


        SS(v) =  b0*exp(-1*ob*og'*DT_mat(:,:)*og);

        %here b and g need to be original uncorrected
        %bvals, and bvecs


        % Most simply, the adjusted bvec is simply L * bvec. Here we are
        % operating in the image space.
        og(1) = -og(1);
        %disp(L_mat);
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

        S(v) = b0*exp(-1*adjbval*adjbvec'*DT_mat(:,:)*adjbvec);

        ag(:,v) = adjbvec; 
        ab(v) = adjbval;
    end
end
	
