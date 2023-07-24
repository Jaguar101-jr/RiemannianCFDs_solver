
%% This code is developed by HUSSEIN ABDUELHALEIM HUSSEIN MUHAMMED (https://github.com/Jaguar101-jr)
%% As a tutorial to implement Riemann solver for CFDs equ.
%% B.Sc.H. University of Khartoum (Sudan) 2015, M.Sc. from UPC-Qingdao (China) 2022.
%% The code implements a basic numerical solver **matlab function for the one-dimensional compressible Euler equations.
%% it demonstrates the solution to a simple Riemann problem using the method of characteristics.

function Riemanniansolver()
    % Define the initial conditions
    rho_L = 1.0;    % Density on the left side of the interface
    rho_R = 0.125;  % Density on the right side of the interface
    u_L = 0.0;      % Velocity on the left side of the interface
    u_R = 0.0;      % Velocity on the right side of the interface
    p_L = 1.0;      % Pressure on the left side of the interface
    p_R = 0.1;      % Pressure on the right side of the interface

    % Define material properties
    gamma = 1.4;    % Ratio of specific heats

    % Time to compute the solution
    t_end = 0.25;

    % Define the grid
    nx = 100;           % Number of grid points
    x_start = -1.0;     % Left boundary
    x_end = 1.0;        % Right boundary
    dx = (x_end - x_start) / (nx - 1); % Grid spacing

    % Initialize the solution arrays
    x = linspace(x_start, x_end, nx);
    rho = zeros(1, nx);
    u = zeros(1, nx);
    p = zeros(1, nx);

    % Set the initial conditions
    for i = 1:nx
        if x(i) <= 0
            rho(i) = rho_L;
            u(i) = u_L;
            p(i) = p_L;
        else
            rho(i) = rho_R;
            u(i) = u_R;
            p(i) = p_R;
        end
    end

    % Time integration using the method of characteristics
    t = 0;
    CFL = 0.5; % Courant-Friedrichs-Lewy number (adjust as needed)
    while t < t_end
        dt = CFL * dx / max(abs(u) + sqrt(gamma * p ./ rho));
        if t + dt > t_end
            dt = t_end - t;
        end

        % Update the solution using the method of characteristics
        for i = 2:nx-1
            rho(i) = rho(i) - dt/dx * (rho(i+1) * u(i+1) - rho(i-1) * u(i-1));
            u(i) = u(i) - dt/dx * (rho(i+1) * u(i+1)^2 + p(i+1) - rho(i-1) * u(i-1)^2 - p(i-1));
            p(i) = p(i) - dt/dx * (u(i+1) * (p(i+1) + rho(i+1) * u(i+1)^2) - u(i-1) * (p(i-1) + rho(i-1) * u(i-1)^2));
        end

        % Apply boundary conditions (reflection at the boundaries)
        rho(1) = rho(2);
        u(1) = -u(2);
        p(1) = p(2);

        rho(nx) = rho(nx-1);
        u(nx) = -u(nx-1);
        p(nx) = p(nx-1);

        t = t + dt;
    end

    % Plot the final solution
    figure;
    plot(x, rho, 'b-', 'LineWidth', 2);
    hold on;
    plot(x, u, 'r-', 'LineWidth', 2);
    plot(x, p, 'g-', 'LineWidth', 2);
    xlabel('displacement X-axis');
    ylabel('Density, Velocity, Pressure');
    legend('Density', 'Velocity', 'Pressure');
    title('CFDs Riemann Solver');
    grid on;
end
