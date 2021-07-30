function op_v = perpVector(ip_v)
% Usage: op_v = perpVector(ip_v)
%
% perpVector finds a vector perpendicular to ip_v and 
% returns it in op_v.​
% Find a vector not parallel to ip_v:
    minIndex = min(find(ip_v == min(ip_v)));    % Index of smallest component of ip_v.
    np_v = zeros(size(ip_v));
    np_v(minIndex) = 1;

    % Find vector perpendicular to ip_v and np_v:
    op_v = cross(ip_v, np_v);

    % Normalize it:
    op_v = op_v / norm(op_v);
end

