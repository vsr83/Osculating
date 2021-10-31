function [r, v] = cartesian(a, ecc_norm, incl, Omega, omega, M, mu)
% CARTESIAN - Compute Cartesian Orbit State Vector from Kepler Parameters
%   
% INPUTS: 
%   a              The Semi-Major Axis (in degrees).
%   ecc_norm       The scalar Eccentricity (0 < ecc_norm < 1).
%   incl           The Inclination (in degrees).
%   Omega          The Longitude of Ascending Node (in degrees).
%   omega          The Argument of Periapsis (in degrees).
%   M              The Mean Anomaly (in degrees).
%   mu             The standard gravitational parameter (in m^3/s^2).
%                  If not filled, parameter for Earth is assumed.
% OUTPUTS:
%   R              Position as a 3x1 vector (in meters).
%   V              Velocity as a 3x1 vector (in meters).

    if (nargin < 7)
        % Standard gravitational parameter for Earth (m^3/s^2).
        mu = 3.986004418e14;
    end 

    % The Eccentric Anomaly
    E = solve_kepler(M, ecc_norm, 1e-9, 20);
    % Semi-Minor Axis
    b = a * sqrt(1 - norm(ecc_norm)^2);

    % Coordinates on the Orbital Plane.
    r_orbital = [a * (cosd(E) - ecc_norm); b * sind(E); 0];

    dEdt = (sqrt(mu) / sqrt(a)^3) / (1 - ecc_norm*cosd(E));
    v_orbital = [-a * dEdt * sind(E); b * dEdt * cosd(E); 0];

    R_Omega = [cosd(Omega), -sind(Omega), 0; sind(Omega), cosd(Omega), 0; 0, 0, 1];
    R_incl = [1, 0, 0; 0, cosd(incl), -sind(incl); 0, sind(incl), cosd(incl)];
    R_omega = [cosd(omega), -sind(omega), 0; sind(omega), cosd(omega), 0; 0, 0, 1];

    % Cartesian coordinates.
    r = R_Omega * R_incl * R_omega * r_orbital;
    v = R_Omega * R_incl * R_omega * v_orbital;
end