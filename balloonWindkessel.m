function B = balloonWindkessel(S, t)

    % -------------------------------------------------------------
    % Balloon-Windkessel hemodynamic model - Equation 23
    %
    % INPUT:
    % S = synaptic activity from Equation 22
    % t = corresponding time vector
    %
    % OUTPUT:
    % B = simulated BOLD signal
    % -------------------------------------------------------------

    % Parameters from Table 3
    kappa = 0.65;
    gamma = 0.41;
    tau_H = 0.98;
    alpha = 0.32;
    rho   = 0.34;
    V0    = 0.02;

    % Initial conditions:
    % s = 0, f = 1, v = 1, q = 1
    y0 = [0; 1; 1; 1];

    % Define the ODE system
    model = @(tCurrent, y) balloonWindkesselODE( ...
        tCurrent, y, t, S, ...
        kappa, gamma, tau_H, alpha, rho);

    % Solve the differential equations
    [~, ySol] = ode45(model, t, y0);

    % Extract blood volume and deoxyhemoglobin
    v = ySol(:,3);
    q = ySol(:,4);

    % Calculate BOLD signal
    B = V0 .* ( ...
          7*rho .* (1 - q) ...
        + 2 .* (1 - q ./ v) ...
        + (2*rho - 0.2) .* (1 - v) );

end


function dydt = balloonWindkesselODE( ...
    tCurrent, y, tInput, S, ...
    kappa, gamma, tau_H, alpha, rho)

    % State variables
    s = y(1);   % vasodilatory signal
    f = y(2);   % normalized blood inflow
    v = y(3);   % normalized blood volume
    q = y(4);   % normalized deoxyhemoglobin content

    % Find S at the current time requested by ode45
    neuralInput = interp1(tInput, S, tCurrent, 'linear', 0);

    % Equation 23a
    dsdt = neuralInput ...
           - kappa*s ...
           - gamma*(f - 1);

    % Equation 23b
    dfdt = s;

    % Equation 23c
    dvdt = (f - v^(1/alpha)) / tau_H;

    % Equation 23d
    dqdt = (1/tau_H) * ( ...
          (f/rho) * (1 - (1-rho)^(1/f)) ...
          - q * v^(1/alpha - 1) );

    % Return derivatives to ode45
    dydt = [dsdt;
            dfdt;
            dvdt;
            dqdt];

end
