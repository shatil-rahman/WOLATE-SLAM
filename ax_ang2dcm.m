function [ DCM ] = ax_ang2dcm( phi)
%AX_ANG2DCM Summary of this function goes here
%   Detailed explanation goes here
phi = reshape(phi,3,1);
angle = norm(phi);

phi_cross = [  0, -phi(3), phi(2);
               phi(3), 0, -phi(1);
               -phi(2), phi(1), 0];

DCM = cos(angle)*eye(3) + (1-cos(angle))*(phi./angle)*(phi./angle)' - sin(angle)*phi_cross./angle;

if (abs(det(DCM)-1)>100*eps)
    error('Determinant not 1, invalid DCM created')
end

if (norm(DCM*DCM' - eye(3)) > 100*eps)
    error('not orthogonal, Invalid DCM created')
end

end

