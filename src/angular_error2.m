function diff = angular_error(PEa, PEb)
	chord = (PEa(:,1) - PEb(:,1).^2 + PEa(:,2) - PEb(:,2).^2 + PEa(:,3) - PEb(:,3).^2 );
	chord = sqrt(chord);
	ang = 2 * real(asin(chord / 2));
	ang(ang > (pi/2)) = pi - ang(ang > (pi/2)); 
	diff = rad2deg(ang);

	normalize = @(v) v ./ sqrt(sum(v.^2));
	diff = acosd(dot(normalize(PEa), normalize(PEb)));
end
