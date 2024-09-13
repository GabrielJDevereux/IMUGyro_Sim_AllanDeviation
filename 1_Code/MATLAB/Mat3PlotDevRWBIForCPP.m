% Load CSV Allan Deviation Data
data = readtable('allan_deviation_output.csv');

% Averaging time (tau)
tau = data{:, 1}; 
% Deviation (sigma)
sigma = data{:, 2};

% Plot deviation
figure;
loglog(tau, sigma, 'b', 'LineWidth', 2);
hold on;
grid on;
xlabel('Averaging Time \tau [s]');
ylabel('Allan Deviation \sigma(\tau)');
title('Allan Deviation computed in C++');

%% Angle Random Walk (ARW)
slope = -0.5;  
logtau = log10(tau);  
logadev = log10(sigma);  
dlogadev = diff(logadev) ./ diff(logtau); 

% Point where the slope is closest to -0.5 / Gaussian White Noise characteristic slope
[~, ARW_idx] = min(abs(dlogadev - slope));


b = logadev(ARW_idx) - slope * logtau(ARW_idx);

% Determine the ARW coefficient
logN = b;  
N = 10^logN;  

% Compute ARW in degrees / root hour
tau_1 = 1; % = 1 second
y_tau1 = sigma(tau == tau_1); % Find the deviation at tau = 1
if isempty(y_tau1)
    y_tau1 = sigma(ARW_idx); % Use the closest tau if exact is not found
    disp('tau is not found')
end
disp(y_tau1);
ARW_deg_per_rt_hr = y_tau1 * 60; % Convert to degrees per root hour

% Plot ARW reference line
lineN = N ./ sqrt(tau); 
loglog(tau, lineN, '--r', 'LineWidth', 1.5);  
text(tau_1, N, 'N', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

%% Bias Instability (BI)
slope = 0;  
[~, BI_idx] = min(abs(dlogadev - slope)); 

b = logadev(BI_idx) - slope * logtau(BI_idx);

% Determine the BI coefficient from the line
scfB = 0.664;  
logB = b - log10(scfB);  
B = 10^logB;  

% Compute BI in degrees per hour
BI_tau = tau(BI_idx);
BI_y_value = sigma(BI_idx);
BI_deg_per_hr = (BI_y_value / 0.664) * 3600; % Convert to degrees per hour

% Plot BI reference line
lineB = B * scfB * ones(size(tau));  
loglog(tau, lineB, '--g', 'LineWidth', 1.5);  
text(BI_tau, scfB*B, '0.664B', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Rate Random Walk (RRW)
slope = 0.5;  
[~, RRW_idx] = min(abs(dlogadev - slope)); 

b = logadev(RRW_idx) - slope * logtau(RRW_idx);

% Determine the RRW coefficient from the line
logK = slope * log10(3) + b;
K = 10^logK; 

% Plot RRW reference line
lineK = K .* sqrt(tau / 3); 
loglog(tau, lineK, '--m', 'LineWidth', 1.5);  
text(3, K, 'K', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Display ARW and BI values 
text(max(tau)*0.1, max(sigma)*0.9, ...
    {sprintf('ARW: %.4f deg/rt-hr', ARW_deg_per_rt_hr), ...
     sprintf('BI: %.4f deg/hr', BI_deg_per_hr)}, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 10, 'FontWeight', 'bold');

% Finalize Legend 
legend('\sigma(\tau)', '\sigma_N (\tau^{-0.5})', '\sigma_B (\tau^{0})', '\sigma_K (\tau^{0.5})', 'Location', 'SouthWest');
hold off;



