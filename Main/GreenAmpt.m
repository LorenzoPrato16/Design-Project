function [GreenAmpt_output] = GreenAmpt(station, timescale, K_s, Psi, theta_i, theta_s)

    % Moisture deficit [-]
    theta_d = theta_s - theta_i;

    % Load appropriate dataset based on station
    if station == "Airolo"
        opts = detectImportOptions('Precipitazioni_airolo2024.txt');
        opts = setvartype(opts, 'time', 'char');
        T = readtable('Precipitazioni_airolo2024.txt', opts);
        T.time = datetime(T.time, 'InputFormat', 'yyyyMMddHH');
    elseif station == "Biasca"
        opts = detectImportOptions('precipitazioni_Biasca2024.txt');
        opts.Delimiter = ';';
        opts = setvartype(opts, 'time', 'char');
        T = readtable('precipitazioni_Biasca2024.txt', opts);
        T.time = datetime(T.time, 'InputFormat', 'yyyyMMddHH');
    end

    % Handle continuous full year case
    if timescale == "Year 2024"
        t = hours(T.time - T.time(1)); % Convert time to hours relative to the start
        precip = T.rre150h0 / 10;  % mm to cm

        % Initialize vectors to store results
        F_sol = zeros(size(t)); % Cumulative infiltration [cm]
        f_star = zeros(size(t)); % Potential infiltration rate [cm/h]
        f_sol = zeros(size(t)); % Actual infiltration rate [cm/h]
        z = zeros(size(t)); % Wetting front depth [cm]

        % Time loop to compute infiltration at each time step
        for i = 1:length(t)
            F = K_s * t(i);  % Initial guess for cumulative infiltration
            tol = 1e-6; % Tolerance for convergence
            max_iter = 100; % Max number of iteration allowed
            iter = 0; % Counter 

             % Define Green-Ampt cumulative infiltration function and its derivative
            F_fun = @(F) F - Psi * theta_d * log(1 + F / (Psi * theta_d)) - K_s * t(i);
            F_prime = @(F) 1 - (Psi * theta_d) / (F + Psi * theta_d);

            % Newton-Raphson iteration to solve for cumulative infiltration F
            while iter < max_iter
                f_val = F_fun(F);
                f_der = F_prime(F);
                F_new = F - f_val / f_der;
                if abs(F_new - F) < tol
                    break;
                end
                F = F_new;
                iter = iter + 1;
            end

            F_sol(i) = F;
            f_star(i) = K_s * (1 + (Psi * theta_d) / F); % Compute potential infiltration rate
            f_sol(i) = min(f_star(i), precip(i)); % Determine actual infiltration rate based on precipitation
            z(i) = F / theta_d; % Compute wetting front depth
        end

        % Store results
        GreenAmpt_output = table(t, f_sol, F_sol, f_star, z, ...
            'VariableNames', {'t[h]', 'f(t)', 'F(t)', 'f*(t)', 'z(t)'});

    % Handle peak/volume month selection
    elseif timescale == "June 2024" || timescale == "March 2024"

         % Initialisation of a cell array to store data for each month
        events = cell(1, 12);
        months_names = {'January', 'February', 'March', 'April', 'May', 'June', ...
                        'July', 'August', 'September', 'October', 'November', 'December'};
        months = month(T.time);
        for m = 1:12
            idx = months == m;
            events{m} = T(idx, {'stn', 'time', 'rre150h0'});
        end

        max_precipitation = zeros(1, 12);
        total_precipitation = zeros(1, 12);

        % Calculate maximum hourly and total precipitation for each month
        for m = 1:12
            month_data = events{m};
            max_precipitation(m) = max(month_data.rre150h0);
            total_precipitation(m) = sum(month_data.rre150h0, 'omitnan');
        end

        % Identify the month with the highest peak precipitation
        [~, month_with_max_precip] = max(max_precipitation);
        % Identify the month with the highest total precipitation
        [~, month_with_max_rain] = max(total_precipitation);

        if timescale == "June 2024"
            selected_month = month_with_max_precip;
        else
            selected_month = month_with_max_rain;
        end

        % Extract data for the selected month
        data = events{selected_month};
        t = hours(data.time - data.time(1));
        precip = data.rre150h0 / 10;

        % Initialize output vectors for infiltration results
        F_sol = zeros(size(t));
        f_star = zeros(size(t));
        f_sol = zeros(size(t));
        z = zeros(size(t));

        % Compute Green-Ampt infiltration for each time step
        for i = 1:length(t)
            F = K_s * t(i);
            tol = 1e-6;
            max_iter = 100;
            iter = 0;
            F_fun = @(F) F - Psi * theta_d * log(1 + F / (Psi * theta_d)) - K_s * t(i);
            F_prime = @(F) 1 - (Psi * theta_d) / (F + Psi * theta_d);

            % Newton-Raphson iteration to solve for F
            while iter < max_iter
                f_val = F_fun(F);
                f_der = F_prime(F);
 
                F_new = F - f_val / f_der;
                if abs(F_new - F) < tol
                    break;
                end
                F = F_new;
                iter = iter + 1;
            end
        
            F_sol(i) = F;
            f_star(i) = K_s * (1 + (Psi * theta_d) / F);
            f_sol(i) = min(f_star(i), precip(i));
            z(i) = F / theta_d;
        end
        
        % Store results of the infiltration model
        GreenAmpt_output = table(t, f_sol, F_sol, f_star, z, ...
            'VariableNames', {'t[h]', 'f(t)', 'F(t)', 'f*(t)', 'z(t)'});

    end
end
