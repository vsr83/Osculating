function [a, ecc_norm, incl, Omega, omega, E, M, f] = osculating(r, v, mu)
% OSCULATING - Compute Osculating Orbit Parameters.
%   
% INPUTS: 
%   R              Position as a 3x1 vector (in meters).
%   V              Velocity as a 3x1 vector (in meters).
%   mu             The standard gravitational parameter (in m^3/s^2).
%                  If not filled, parameter for Earth is assumed.
% OUTPUTS:
%   E              The Eccentric Anomaly (in degrees).

incl_min = 1e-7;

if (nargin < 3)
    % Standard gravitational parameter for Earth (m^3/s^2).
    mu = 3.986004418e14;
end 

% Angular momentum per unit mass.
k = cross(r, v);

% Eccentricity vector.
ecc = cross(v, k)/mu - r/norm(r);
ecc_norm = norm(ecc);

% Inclination
incl = acosd(k(3) / norm(k));

% Energy integral.
h = 0.5 * norm(v) * norm(v) - mu / norm(r);

% Semi-major axis.
a = -mu / (2*h);
b = a * sqrt(1 - norm(ecc)^2);

% Longitude of the ascending node.
Omega = atan2d(k(1), -k(2));

% Argument of periapsis.
% When inclination is close to zero, we wish to compute the argument of periapsis
% on the XY-plane and avoid division by zero:
if incl < incl_min
    omega = atan2d(ecc(2), ecc(1)) - Omega;
else
    % We wish to avoid division by zero and thus use the formula, which has larger
    % absolute value for the denominator:
    if abs(sind(Omega)) < abs(cosd(Omega))
        asc_y = ecc(3) / sind(incl);
        asc_x = (1 / cosd(Omega)) * (ecc(1) + sind(Omega)*cosd(incl)*ecc(3) / sind(incl));

        omega = atan2d(asc_y, asc_x);
    else
        asc_y = ecc(3) / sind(incl);
        asc_x = (1 / sind(Omega)) * (ecc(2) - cosd(Omega)*cosd(incl)*ecc(3) / sind(incl));

        omega = atan2d(asc_y, asc_x)
    end
end

% Eccentric anomaly
R_Omega = [cosd(Omega), -sind(Omega), 0; sind(Omega), cosd(Omega), 0; 0, 0, 1];
R_incl = [1, 0, 0; 0, cosd(incl), -sind(incl); 0, sind(incl), cosd(incl)];
R_omega = [cosd(omega), -sind(omega), 0; sind(omega), cosd(omega), 0; 0, 0, 1];

r_orbital = R_omega' * R_incl' * R_Omega' * r;
E = atan2d(r_orbital(2)/b, r_orbital(1)/a + norm(ecc));

% Mean anomaly
M = E - ecc_norm*sind(E);

% Natural anomaly
xu = (cosd(E) - ecc_norm) / (1 - ecc_norm*cosd(E));
yu = sqrt(1 - ecc_norm * ecc_norm) * sind(E) / (1 - ecc_norm * cosd(E));

f = atan2d(yu, xu);

end