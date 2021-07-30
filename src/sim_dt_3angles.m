% from derek jones paper

function DT = sim_dt(phi,theta,psi)
    TrD = 0.0021; 
    FA = 0.3;

    lp = (TrD/3) * (1 + (2*FA) / sqrt( 3 - (2 * FA^2)));

    l2 = (TrD/3) * (1 - FA / sqrt( 3 - (2 * FA^2)));

    l3 = l2;

    D = [lp 0 0 ; 0 l2 0 ; 0 0 l3];

    R = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]*[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)] * [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]; 

    DT = R * D * R';
end
