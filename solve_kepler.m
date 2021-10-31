function E = solve_kepler(M, ecc, tolerance, max_iterations)
% SOLVE_KEPLER - Solve the Kepler equation
%   
% INPUTS: 
%   M              The Mean Anomaly (in degrees).
%   ECC            The orbit eccentricity (0 < ECC < 1).
%   TOLERANCE      The tolerance.
%   MAX_ITERATIONS The maximum number of iterations.
% OUTPUTS:
%   E              The Eccentric Anomaly (in degrees).

    iteration_count = 0;
    error = tolerance + 1.0;

    E = M;

    while (error > tolerance)
        iteration_count = iteration_count + 1;

        if (iteration_count > max_iterations)
            error("Convergence failed");
        end

        % Newton-Raphson iteration:
        E = E - (E - ecc * sind(E) - M) / (1.0 - ecc * cosd(E));

        % An issue here might be that the iteration might converge to 
        % modulo angle of the target angle so we avoid using angles.
        error = abs(sind(M) - sind(E - ecc * sind(E)));
    end
end