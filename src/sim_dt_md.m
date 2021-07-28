%l1 = 30;
%lp = 10;
%l2 = lp;
%l3 = lp;
%lh = (l1 + l2 + l3)/3;

%FA = (3/2)^2 * (sqrt((l1-lh)^2 + (l2-lh)^2 + (l3-lh)^2) ./ sqrt(l1^2 + l2^2 + l3^2));

function D = sim_dt(phi,theta,psi) %sharpiness and MD
    MD = 0.0007; %average diffusion
    l1 = 0.0017; %2.8; %randomly picked
    lp = ((3 * MD) - l1)/2;

    %phi = 90;
    %theta = 90;
    %psi = 90;
    R = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]*[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)] * [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]; 

    D = R * [ l1 0 0; 0 lp 0; 0 0 lp] * R';
end
