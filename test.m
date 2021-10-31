% Load ISS Orbit State Vectors.
load iss_osv.txt

n = 150;
%iss_osv = iss_osv(1:n, :);

% Convert to meters.
iss_osv = iss_osv * 1000;
[num_osv, num_col] = size(iss_osv);

a_list = zeros(num_osv, 1);
ecc_list = zeros(num_osv, 1);
incl_list = zeros(num_osv, 1);
Omega_list = zeros(num_osv, 1);
omega_list = zeros(num_osv, 1);
E_list = zeros(num_osv, 1);
M_list = zeros(num_osv, 1);
f_list = zeros(num_osv, 1);
err_pos = zeros(num_osv, 1);
err_vel = zeros(num_osv, 1);

for ind_osv = 1:num_osv
    r_osv = iss_osv(ind_osv, 1:3)';
    v_osv = iss_osv(ind_osv, 4:6)';

    [a, ecc, incl, Omega, omega, E, M, f] = osculating(r_osv, v_osv);

    a_list(ind_osv) = a;
    ecc_list(ind_osv) = ecc;
    incl_list(ind_osv) = incl;
    Omega_list(ind_osv) = Omega;
    omega_list(ind_osv) = omega;
    E_list(ind_osv) = E;
    M_list(ind_osv) = M;
    f_list(ind_osv) = f;

    [r_check, v_check] = cartesian(a, ecc, incl, Omega, omega, M);

    err_pos(ind_osv) = norm(r_osv - r_check);
    err_vel(ind_osv) = norm(v_osv - v_check);
end

disp 'Drawing.'

figure(1);
clf
subplot(4, 3, 1)
plot(1:n, iss_osv(1:n, 1), 'r'); 
hold on
plot(1:n, iss_osv(1:n, 2), 'g'); 
plot(1:n, iss_osv(1:n, 3), 'b');
legend 'x' 'y' 'z' 
xlabel 'OSV'
ylabel 'Position (m)'
xlim([0 n])
subplot(4, 3, 4)
plot(1:n, iss_osv(1:n, 4), 'r'); 
hold on
plot(1:n, iss_osv(1:n, 5), 'g'); 
plot(1:n, iss_osv(1:n, 6), 'b');
legend 'v_x' 'v_y' 'v_z' 
xlabel 'OSV'
ylabel 'Velocity (m/s)'
xlim([0 n])
subplot(4, 3, 7)
plot(1:n, err_pos(1:n), 'r'); 
xlabel 'OSV'
ylabel 'Position Error (m)'
xlim([0 n])
subplot(4, 3, 10)
plot(1:n, err_vel(1:n), 'r'); 
xlabel 'OSV'
ylabel 'Velocity Error (m/s)'
xlim([0 n])

subplot(4, 3, 2)
plot(1:n, a_list(1:n))
xlabel 'OSV'
ylabel 'Semi-Major Axis (m)'
xlim([0 n])
title 'ISS Orbit Parameters'
subplot(4, 3, 5)
plot(1:n, ecc_list(1:n))
xlabel 'OSV'
ylabel 'Eccentricity'
xlim([0 n])
subplot(4, 3, 8)
plot(1:n, incl_list(1:n))
xlabel 'OSV'
ylabel 'Inclination (deg)'
xlim([0 n])
subplot(4, 3, 11)
plot(1:n, Omega_list(1:n))
xlabel 'OSV'
ylabel 'Longitude of Ascending Node (deg)'
xlim([0 n])

subplot(4, 3, 3)
plot(1:n, omega_list(1:n))
xlabel 'OSV'
ylabel 'Argument of Periapsis (deg)'
xlim([0 n])
subplot(4, 3, 6)
plot(1:n, E_list(1:n))
xlabel 'OSV'
ylabel 'Eccentric Anomaly (deg)'
xlim([0 n])
subplot(4, 3, 9)
plot(1:n, M_list(1:n))
xlim([0 n])
xlabel 'OSV'
ylabel 'Mean Anomaly (deg)'

subplot(4, 3, 12)
plot(1:n, f_list(1:n) + Omega_list(1:n) + omega_list(1:n))
ylim([-0 360])
xlim([0 n])
xlabel 'OSV'
ylabel 'Mean Longitude (deg)'

figure(2)
clf
for ind_osv = 1:200:num_osv 
    a = a_list(ind_osv);
    ecc = ecc_list(ind_osv);
    incl = incl_list(ind_osv);
    Omega = Omega_list(ind_osv);
    omega = omega_list(ind_osv);
    
    r_draw = zeros(361, 3);

    disp(ind_osv)

    for ind_M = 1:361
        M = ind_M;
        [r_out, v_out] = cartesian(a, ecc, incl, Omega, omega, M);
        r_draw(ind_M, :) = r_out';
    end
    plot3(r_draw(:, 1), r_draw(:, 2), r_draw(:, 3));
    hold on;
    axis equal;
end
xlabel 'x (m)'
ylabel 'y (m)'
zlabel 'z (m)'
title 'Reconstructed Orbits in Inertial Coordinates'