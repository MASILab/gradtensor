function new_dwi_signal = bl_correct_signal_sm(S_corpt,g,b,abvec,abval,L_mat) 
	b0 = 1;
	%L_mat = [ 1.0091 -0.0027 0.0070; -0.0031 1.0025 -0.0018 ; -0.0161 0.0041 1.0042 ];
%     g(1,:) = -g(1,:);
%          
%     gg = L_mat * g;
%          %norm_gg = norm(gg);
%     len2 = sum(gg.^2);
    
	for v = 1:length(b)
         og = g(:,v);
         ob = b(v);
% 
%         %here b and g need to be original uncorrected
%         %bvals, and bvecs
%         % Most simply, the adjusted bvec is simply L * bvec. Here we are
%         % operating in the image space.
          %og(1) = og(1);
%          
          gg = L_mat * og;
%          %norm_gg = norm(gg);
          len2 = sum(gg.^2);
          %len2 = norm(gg);
        %new_dwi_signal(v) = S(v) * exp(-1*ggabvec(:,v)'*eye(3)*abvec(:,v));
        new_dwi_signal(v) = b0 * exp( (log(S_corpt(v)/b0)) / len2 );
        
	end
end
