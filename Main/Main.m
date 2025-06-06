clear all

% Definition of Green and Ampt parameters:
K_s = 0.7;           % [cm/h] hydraulic conductivty at saturation
Psi = 11.6;          % [cm] suction
theta_i = 5.6/100;   % [-] initial water content
theta_s = 15.5/100;  % [-] saturated water content (porosity)

C0=200;              % Initial concentration of pollutant [mg/kg] [200 for arsenic, 500 for nitrate]
H = 20;              % Total deposit height [m]
p=2;                 % Pollutant type [1=nitrate, 2=Arsenic]
tr=1;                % Timestep [h]
rho= 2670;           % Materials density [kg/m3]

station = "Airolo";         %  Change between "Biasca" and "Airolo"
timescale= "March 2024";     % Change between "June 2024", "March 2024" or "Year 2024"

GreenAmpt_output = GreenAmpt(station, timescale, K_s, Psi, theta_i, theta_s);

[Caq_n_filtered,time_filtered, H, p, tr, station, timescale] = Diffusion(H,p,tr,station,timescale, rho, GreenAmpt_output,C0);



% Pollutant-specific labels
if p == 1
    pollutant_name = 'Nitrate';
    color = [0.1 0.5 0.8];  % bluish
elseif p == 2
    pollutant_name = 'Arsenic';
    color = [0.8 0.2 0.2];  % reddish
end

% Plot leachate concentration
figure('Color', 'w')
hold on
scatter(time_filtered, Caq_n_filtered, 40, 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k')
plot(time_filtered, Caq_n_filtered, '-', 'LineWidth', 2, 'Color', color)

xlabel('Time [hours]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel(['$C_{\mathrm{aq}}$ [mg/L]'], 'Interpreter', 'latex', 'FontSize', 14)
title([pollutant_name ' concentration in Leachate '], 'Interpreter', 'latex', 'FontSize', 16 )

grid on
box on
set(gca, 'FontSize', 12, 'LineWidth', 1.2)
xlim([min(time_filtered) max(time_filtered)])
ylim([0 max(Caq_n_filtered)*1.1])

% Display a box with H and tr
info_text = {sprintf('H (Total Depth) = %d m', H), sprintf('Time Step = %d h', tr), sprintf('Station: %s', station),sprintf("Timescale: %s", timescale)};
annotation('textbox', [0.45 0.8 0.2 0.1], 'String', info_text, ...
    'FitBoxToText', 'on');
