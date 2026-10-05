function B = balloonWindkessel(S, t, label)

    % -------------------------------------------------------------
    % Balloon-Windkessel hemodynamic model - Equation 23
    %
    % INPUT:
    % S     = synaptic activity from Equation 22
    % t     = corresponding time vector
    % label = optional string shown on the progress bar (default: 'Balloon-Windkessel')
    %
    % OUTPUT:
    % B = simulated BOLD signal
    % -------------------------------------------------------------

    if nargin < 3 || isempty(label)
        label = 'Balloon-Windkessel';
    end

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

    % ---- progress bar setup -------------------------------------
    t0 = t(1);
    tend = t(end);
    h = waitbar(0, sprintf('%s: 0%%', label), 'Name', label);
    cleanupObj = onCleanup(@() close(h));   % guarantees the bar closes even on error/Ctrl-C

    n_calls = 0;
    update_every = 200;   % throttle rendering -- updating every single accepted
                          % ODE step would noticeably slow down ode45 itself

    opts = odeset('OutputFcn', @progressFcn);
    % ---------------------------------------------------------------

    % Define the ODE system
    model = @(tCurrent, y) balloonWindkesselODE( ...
        tCurrent, y, t, S, ...
        kappa, gamma, tau_H, alpha, rho);

    % Solve the differential equations
    [~, ySol] = ode45(model, t, y0, opts);

    % Extract blood volume and deoxyhemoglobin
    v = ySol(:,3);
    q = ySol(:,4);

    % Calculate BOLD signal
    B = V0 .* ( ...
          7*rho .* (1 - q) ...
        + 2 .* (1 - q ./ v) ...
        + (2*rho - 0.2) .* (1 - v) );

    % ---- nested OutputFcn: shares workspace with the function above ----
    function status = progressFcn(tCurrent, ~, flag)
        status = 0;
        if isempty(flag)
            n_calls = n_calls + 1;
            frac = (tCurrent(end) - t0) / (tend - t0);
            if mod(n_calls, update_every) == 0 || frac >= 1
                waitbar(min(max(frac,0),1), h, sprintf('%s: %.0f%%', label, 100*frac));
            end
        end
    end

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
