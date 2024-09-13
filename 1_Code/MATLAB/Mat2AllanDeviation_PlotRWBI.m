% Load csv file
data = readtable('gyro_data.csv');  %assumed header

% Extract x values
gyro_x = data{:, 'Gyro_X'}; 

% Euler integration to convert angular rate 
t0 = 0.01; 
theta = cumsum(gyro_x) * t0; 

% Set parameters for Allan Deviation calculation
maxNumM = 100;  
L = length(theta);  
maxM = 32768; 

m = [0.01; logspace(log10(1), log10(maxM), maxNumM).'];  
m = unique(m); 

% Calculate the time interval
tau = m * t0;

% Allocate array 
avar = zeros(numel(m), 1);

% Calculate Allan variance for different values of m
for i = 1:numel(m)
    mi = round(m(i));  

    % Ensure mi is valid
    if mi > 0 && (2 * mi) <= L
        % Calculate the sum of squared differences for the current m
        sum_sq_diff = sum((theta(1+2*mi:L) - 2*theta(1+mi:L-mi) + theta(1:L-2*mi)).^2);

        % Store for variance calculation
        avar(i,:) = sum_sq_diff;

    end
end

% Normalize
avar = avar ./ (2 * tau.^2 .* (L - 2*m));

% Compute deviation (square root of variance)
adev = sqrt(avar);

% Debugging output 
for i = 1:numel(tau)
    if tau(i) >= 239 && tau(i) <= 400  % Adjust the range to print for desired tau values
        fprintf('Final Allan deviation at tau = %.2f is %.15f\n', tau(i), adev(i));
    end
end

% Plot the Allan deviation
figure
loglog(tau, adev, 'b', 'LineWidth', 2)
grid on
title('Allan Deviation computed in MATLAB')
xlabel('\tau (s)')
ylabel('\sigma(\tau) [rad]')
hold on

%% Angle Random Walk (ARW)
slope_arw = -0.5;  % slope
logtau = log10(tau); 
logadev = log10(adev);  
dlogadev = diff(logadev) ./ diff(logtau);  

% Find where the slope is closest to -0.5 
[~, ARW_idx] = min(abs(dlogadev - slope_arw));

% Find the y-intercept
b = logadev(ARW_idx) - slope_arw * logtau(ARW_idx);

% Angle random walk coefficient
logN = b;  
N = 10^logN; 

% Tau = 1 second
tau_1 = 1;

% Find the closest value to tau = 1 in the tau array
[~, idx] = min(abs(tau - tau_1));

% Get the Allan deviation at the closest tau value
y_tau1 = adev(idx);

% Is the closest tau is significantly different from tau = 1
if abs(tau(idx) - tau_1) > 1e-6  d
    disp('Exact tau = 1 is not found, using closest value');
end

disp(['Closest tau: ', num2str(tau(idx))]);
disp(['Allan deviation at closest tau: ', num2str(y_tau1)]);

% Compute ARW in degrees per root hour
ARW_deg_per_rt_hr = y_tau1 * 60; 


% Plot ARW reference line
lineN = N ./ sqrt(tau);  
loglog(tau, lineN, '--r', 'LineWidth', 1.5);  
text(tau_1, N, 'N', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');


%% Bias Instability (BI) 
slope_bi = 0;  % Slope 
logtau = log10(tau);  
logadev = log10(adev);  
dlogadev = diff(logadev) ./ diff(logtau); 

% Find where the slope is closest to 0
[~, BI_idx] = min(abs(dlogadev - slope_bi));

b = logadev(BI_idx) - slope_bi * logtau(BI_idx);

% Calculate BI from the y-intercept
scfB = 0.664; 
logB = b - log10(scfB);  
B = 10^logB;  

% Compute BI in degrees per hour
BI_tau = tau(BI_idx);  
BI_y_value = adev(BI_idx);  
BI_deg_per_hr = (BI_y_value / 0.664) * 3600;  % Convert to degrees per hour

% Plot Bias Instability reference line
lineB = B * scfB * ones(size(tau));  
loglog(tau, lineB, '--g', 'LineWidth', 1.5); 

% Display Bias Instability value
text(BI_tau, scfB * B, sprintf('BI: %.4f deg/hr', B * 3600), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');


%% Rate Random Walk (RRW)
slope_rrw = 0.5;  % slope
logtau = log10(tau); 
logadev = log10(adev);  
dlogadev = diff(logadev) ./ diff(logtau); 

% Find where the slope is closest to 0.5
[~, RRW_idx] = min(abs(dlogadev - slope_rrw));

b = logadev(RRW_idx) - slope_rrw * logtau(RRW_idx);

% Determine the RRW coefficient
logK = slope_rrw * log10(3) + b;  
K = 10^logK; 

% Plot RRW reference line
lineK = K .* sqrt(tau / 3); 
loglog(tau, lineK, '--m', 'LineWidth', 1.5); 
text(3, K, 'K', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');


%% Display ARW and BI 
annotation('textbox', [0.15, 0.75, 0.2, 0.1], 'String', {sprintf('ARW: %.4f deg/rt-hr', ARW_deg_per_rt_hr), sprintf('BI: %.4f deg/hr', BI_deg_per_hr)}, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 10, 'FontWeight', 'bold', 'EdgeColor', 'none');


%% Finalize the plot 
legend('\sigma(\tau)', 'ARW', 'BI', 'RRW', 'Location', 'SouthWest');
hold off;




