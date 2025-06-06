function [Caq_n_filtered,time_filtered, H, p, tr, station, timescale] = Diffusion(H,p, tr, station, timescale, rho, GreenAmpt_output,C0)

%% Pollutant properties
 % Observed diffusivity [Nitrate Arsenic]
 Dobs=[1.9e-9 1.97e-9]; %[m2/s]

 if p==1
     Dobs=Dobs(1);
 else
     Dobs=Dobs(2);
 end


 %Soil-water partitioning coefficient kd [Nitrate Arsenic]
 kd= [0.1 6.7]; %[L/kg]

 if p==1
     kd=kd(1);
 else
     kd=kd(2);
 end

%Definition of deposit surface 
S=100;




%Extracting results

zt = GreenAmpt_output.("z(t)")/100;      % cumulative infiltration [m]
f_t = GreenAmpt_output.("f(t)");      % infiltration rate [cm/h]
time = GreenAmpt_output.("t[h]");         % time [h]

%Removing Nan values
time(isnan(time)) = 0;
f_t(isnan(f_t)) = 0;
zt(isnan(zt)) = 0;


%Interpolation functions for infiltration depth z[m]
z_interp = @(tt) interp1(time, zt, tt, 'linear', 'extrap');


%% Discretization in layers
t0 = time(1); %beginning of timestep

depth = []; %This will store the depth of each layer

while true
    t1 = t0 + tr; %end of timestep
    z0 = z_interp(t0);
    z1 = z_interp(t1);
    h_i = z1 - z0; %Depth of layer i

    if z0 >= H  %If the depth of the deposit is depassed, stop discretization
        break;
    end
    if z1 > H %depth of last layer
        h_i = H - z0;
        t1 = interp1(zt, time, H); % adjust to exact time water reaches H
    end

    depth(end+1) = h_i;
    

    if z1 >= H %Stop discretization if this condition is met
        break;
    end

    t0 = t1;
end

n = length(depth);  % number of layers



%Remotion of initial zero values
nonzero_start_index = find(f_t > 0, 1, 'first');

% Adjust time and infiltration rate data in order to avoid initial zeros
time_trimmed = time(nonzero_start_index:end);
f_t_trimmed = f_t(nonzero_start_index:end);

% Update precipitation duration based on trimmed time
precip_duration = time_trimmed(end) - time_trimmed(1);
m = round(precip_duration / tr) + 1;  % timesteps during precipitation

% Total simulation time including percolation through all layers
m_sim = m + n - 1;

% Interpolation + Averaging Based on tr
if tr <= 1
    % Interpolate directly using linear interpolation
    t_sim = time_trimmed(1):tr:(time_trimmed(1) + (m_sim - 1) * tr);
    f_raw = interp1(time_trimmed, f_t_trimmed, t_sim, 'linear', 'extrap');
else
    % Defining full simulation time vector
    t_sim = time_trimmed(1):tr:(time_trimmed(1) + (m_sim - 1) * tr);

    %Interpolate f at 0.5h timestep (high resolution)
    t_highres = time_trimmed(1):0.5:(time_trimmed(end));
    f_highres = interp1(time_trimmed, f_t_trimmed, t_highres, 'linear', 'extrap');

    %Aggregate interpolated value into tr blocks
    f_raw = zeros(size(t_sim));
    for i = 1:length(t_sim)
        t_start = t_sim(i);
        t_end = t_start + tr;

        %Get all high-res points within the interval
        idx = t_highres >= t_start & t_highres < t_end;
        if any(idx)
            f_raw(i) = mean(f_highres(idx)); %compute the mean
        else
            f_raw(i) = 0;  %if no data in the selected window
        end
    end
end

%Avoid zero values
f_raw(f_raw < 10^(-6)) = 0.01;

%Final infiltration vector
f = f_raw;

%Extending infiltration rate after precipitation ends
precip_end_time = time_trimmed(end);
precip_end_idx = find(t_sim > precip_end_time, 1, 'first');
if isempty(precip_end_idx)
    precip_end_idx = length(t_sim);
end

% Find last non-zero value before end of precipitation
last_nonzero_idx = find(f_raw(1:precip_end_idx - 1) > 0, 1, 'last');
last_nonzero_value = f_raw(last_nonzero_idx);

% Extend constant infiltration rate after precipitation ends
f(precip_end_idx:end) = last_nonzero_value;

%% Diffusion model

%Create matrix to store results
M_mass=zeros(n,m+n-1); %Cumulative mass released [mg/kg]
Ws=zeros(n,m+n-1); %Water to solid ratio [L/kg]
Caq=zeros(n,m+n-1); %Aqueous concentration of pollutant [mg/L]
Ceq=zeros(n,m+n-1); %Aqueous equilibrium concentration [mg/L]
C0=C0.*ones(n,1);

%Start the simulation
for j=1:m+n-1
    if j < n %when precipitation didn't reach the bottom layer yet. 
        for i = 1:j
            if i==1 %first layer, no initial aqueous concentration
                V=S*depth(i); %volume
                M_mass(i,j)= 2 *C0(i) * S/V * (Dobs * tr *3600 /pi)^(1/2); %mass released 
                Ws(i,j)=10*f(j)*tr/(rho*depth(i)); %Water to solid ratio
                Caq(i,j)=M_mass(i,j)/Ws(i,j); %Aqueous concentration
                Ceq(i,j)= (C0(i)-M_mass(i,j))/kd; %Equilibrium concentration
                if Ceq(i,j)<10^(-6) 
                    Ceq(i,j)=0;
                end
               
                if Caq(i,j)>Ceq(i,j) %If aqueous concentration is higher than the equilibrium one, they will be equal
                    Caq(i,j)=Ceq(i,j);
                end
                C0(i)=C0(i)-M_mass(i,j); %Updating leachable content
                if C0(i)<10^(-6)
                    C0(i)=0;
                end
            else    %layers i>1
                
                M_mass_tot = 0;
                for k = 1:i-1
                    M_mass_tot = M_mass_tot + M_mass(k, j - (i - k)); %cumulative mass release in previous timesteps
                end

                if M_mass_tot>C0(i)
                    M_mass_tot=C0(i);
                end
                V=S*depth(i);
                M_mass(i,j)= 2*(C0(i)-M_mass_tot) *S/V * (Dobs *tr *3600 /pi)^(1/2); %includes concentration gradient between solid and aqueous phase
                Ceq(i,j)= (C0(i)-M_mass(i,j)) / kd;
                if Ceq(i,j)<10^(-6)
                    Ceq(i,j)=0;
                end

                Ws(i,j)=10*f(j)*tr/(rho*depth(i));
                Caq_prev = Caq(i-1, j-1);
                Caq(i,j)=Caq_prev+M_mass(i,j)/Ws(i,j);


                if (Caq(i,j)) > Ceq(i,j)
                    Caq(i,j)= Ceq(i,j);

                end
                C0(i)=C0(i)-M_mass(i,j);
                if C0(i)<10^(-6)
                    C0(i)=0;
                end
            end
        end
    elseif j <= m %Timesteps within precipitation period. Water reaches bottom layer, recording of leached concentration starts
         for i=1:n
            if i==1
                V=S*depth(i);
                M_mass(i,j)= 2 *C0(i) * S/V * (Dobs * tr *3600 /pi)^(1/2);
                Ws(i,j)=10*f(j)*tr/(rho*depth(i));
                Caq(i,j)=M_mass(i,j)/Ws(i,j);
                Ceq(i,j)= (C0(i)-M_mass(i,j))/kd;
                if Ceq(i,j)<10^(-6)
                    Ceq(i,j)=0;
                end
                if Caq(i,j)>Ceq(i,j)
                    Caq(i,j)=Ceq(i,j);
                end
                C0(i)=C0(i)-M_mass(i,j);
                if C0(i)<10^(-6)
                    C0(i)=0;
                end
            else   
                
                M_mass_tot = 0;
                for k = 1:i-1
                    M_mass_tot = M_mass_tot + M_mass(k, j - (i - k));
                end

                if M_mass_tot>C0(i)
                    M_mass_tot=C0(i);
                end
                V=S*depth(i);
                M_mass(i,j)= 2*(C0(i)-M_mass_tot) *S/V * (Dobs *tr *3600 /pi)^(1/2);
                Ceq(i,j)= (C0(i)-M_mass(i,j)) / kd;
                if Ceq(i,j)<10^(-6)
                    Ceq(i,j)=0;
                end

                Ws(i,j)=10*f(j)*tr/(rho*depth(i));
                
                Caq_prev = Caq(i-1, j-1);
                Caq(i,j)=Caq_prev+M_mass(i,j)/Ws(i,j);

                if (Caq(i,j)) > Ceq(i,j)
                    Caq(i,j)= Ceq(i,j);
                    
                end
                C0(i)=C0(i)-M_mass(i,j);
                if C0(i)<10^(-6)
                    C0(i)=0;
                end
            end
        end

    else
        for i = j-m+1:n %timesteps over precipitation period. Water stops flowing in top layers
           
                M_mass_tot = 0;
                for k = 1:i-1
                    M_mass_tot = M_mass_tot + M_mass(k, j - (i - k));
                end

                if M_mass_tot>C0(i)
                    M_mass_tot=C0(i);
                end
                V=S*depth(i);
                M_mass(i,j)= 2*(C0(i)-M_mass_tot) *S/V * (Dobs *tr*3600 /pi)^(1/2);
                Ceq(i,j)= (C0(i)-M_mass(i,j)) / kd;
                if Ceq(i,j)<10^(-6)
                    Ceq(i,j)=0;
                end

                Ws(i,j)=10*f(j)*tr/(rho*depth(i));
                
                Caq_prev = Caq(i-1, j-1);
                Caq(i,j)=Caq_prev+M_mass(i,j)/Ws(i,j);

                if (Caq(i,j)) > Ceq(i,j)
                    Caq(i,j)= Ceq(i,j);
                    
                end
                C0(i)=C0(i)-M_mass(i,j);
                if C0(i)<10^(-6)
                    C0(i)=0;
                end
        end
    end
   
end




% Time vector in hours
time_plot = (1:(m+n-1)) * tr;

% Extract concentration for bottom layer
Caq_n = Caq(n, :);

% Filter out initial zero values
nonzero_idx = Caq_n > 0;
Caq_n_filtered = Caq_n(nonzero_idx);
time_filtered = time_plot(nonzero_idx);




end