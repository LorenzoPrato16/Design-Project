function [GreenAmpt_output] = GreenAmpt(station, timescale, K_s, Psi, theta_i, theta_s)

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
        t = hours(T.time - T.time(1));
        precip = T.rre150h0 / 10;  % mm to cm

        F_sol = zeros(size(t));
        f_star = zeros(size(t));
        f_sol = zeros(size(t));
        z = zeros(size(t));

        for i = 1:length(t)
            F = K_s * t(i);
            tol = 1e-6;
            max_iter = 100;
            iter = 0;
            F_fun = @(F) F - Psi * theta_d * log(1 + F / (Psi * theta_d)) - K_s * t(i);
            F_prime = @(F) 1 - (Psi * theta_d) / (F + Psi * theta_d);

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

        GreenAmpt_output = table(t, f_sol, F_sol, f_star, z, ...
            'VariableNames', {'t[h]', 'f(t)', 'F(t)', 'f*(t)', 'z(t)'});

    % Handle peak/volume month selection
    elseif timescale == "June 2024" || timescale == "March 2024"
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
        for m = 1:12
            month_data = events{m};
            max_precipitation(m) = max(month_data.rre150h0);
            total_precipitation(m) = sum(month_data.rre150h0, 'omitnan');
        end

        [~, month_with_max_precip] = max(max_precipitation);
        [~, month_with_max_rain] = max(total_precipitation);

        if timescale == "June 2024"
            selected_month = month_with_max_precip;
        else
            selected_month = month_with_max_rain;
        end

        data = events{selected_month};
        t = hours(data.time - data.time(1));
        precip = data.rre150h0 / 10;

        F_sol = zeros(size(t));
        f_star = zeros(size(t));
        f_sol = zeros(size(t));
        z = zeros(size(t));

        for i = 1:length(t)
            F = K_s * t(i);
            tol = 1e-6;
            max_iter = 100;
            iter = 0;
            F_fun = @(F) F - Psi * theta_d * log(1 + F / (Psi * theta_d)) - K_s * t(i);
            F_prime = @(F) 1 - (Psi * theta_d) / (F + Psi * theta_d);

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

        GreenAmpt_output = table(t, f_sol, F_sol, f_star, z, ...
            'VariableNames', {'t[h]', 'f(t)', 'F(t)', 'f*(t)', 'z(t)'});

    else
        % Default fallback for unknown or custom timescale
        t = hours(T.time - T.time(1));
        precip = T.rre150h0 / 10;

        F_sol = zeros(size(t));
        f_star = zeros(size(t));
        f_sol = zeros(size(t));
        z = zeros(size(t));

        for i = 1:length(t)
            F = K_s * t(i);
            tol = 1e-6;
            max_iter = 100;
            iter = 0;
            F_fun = @(F) F - Psi * theta_d * log(1 + F / (Psi * theta_d)) - K_s * t(i);
            F_prime = @(F) 1 - (Psi * theta_d) / (F + Psi * theta_d);

            while iter < max_iter
                f_val = F_fun(F);
                f_der = F_prime(F);
                if abs(f_der) < eps
                    warning('Derivative near zero at step %d', i);
                    break;
                end
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

        GreenAmpt_output = table(t, f_sol, F_sol, f_star, z, ...
            'VariableNames', {'t[h]', 'f(t)', 'F(t)', 'f*(t)', 'z(t)'});
    end
end
