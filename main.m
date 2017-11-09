clear all;
clc;
close all;

% error to the power (n)^-n

v_inf = 130 * (5280/3600); % 130 mph
af_r  = 4412;
af_t  = 2412;
b     = 40;
c_r   = 6;
c_t   = 2;
geo_r = 7;
geo_t = 4;
N = 1000;
[a0_r, aero_r] = NACA_lift_slope(af_r, c_r);
[a0_t, aero_t] = NACA_lift_slope(af_t, c_t);

% find true e, cL, cDi
[e_true, c_L_true, c_Di_true] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, 2000);

% computer percent error for N values below
N_range = [2, 5, 7, 10, 15, 20, 25];
errs = zeros(1, length(N_range));
for N = N_range
  [e, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);
  errs(N == N_range) = 100*abs(c_L_true - c_L)/c_L_true;
end

% best fit line and plot
f = fit(N_range', errs', 'cubicinterp');
figure; hold on; grid on;
scatter(N_range, errs);
x = linspace(min(N_range), max(N_range), 100);
y = feval(f, x);
plot(x, y);
title('Number of Vortexes vs Percent Error');
xlabel('Number of Vortexes');
ylabel('Percent Error');

% minimum N for error less than the values below
for e = [5, 1, 1/10]
  idx = find(y < e, 1);
  x_e = x(idx);
  scatter(x_e, e, 'ro');
  fprintf('%.1f%% error: %d points\n', e, ceil(x_e));
end
print('n_vs_err', '-dpng');

% problem 2 (or 3, depending on how you count)
AR_range = [4, 6, 8, 10];
taper_range = linspace(0, 1, 100);
e_vals = zeros(length(AR_range), length(taper_range));
for AR = AR_range
  for t = taper_range
    % b = 20; % ft
    % c_r = 2*b/(AR*(1+t));
    c_r = 6;
    c_t = t*c_r;
    b = AR*((1+t)*c_r)/2;

    % lift slope = 2pi, zero twist, alpha_L0 = 0, alpha = 5, N = 50
    [e, ~, ~] = PLLT(b, 2*pi, 2*pi, c_t, c_r, 0, 0, 5, 5, 50);
    e_vals(AR == AR_range, t == taper_range) = e;
  end
end

figure; hold on;
xlabel('Taper Ratio C_t/C_r');
ylabel('Span efficiency');
for AR = AR_range
  plot(taper_range, e_vals(AR == AR_range, :), ...
       'DisplayName', sprintf('AR = %d', AR));
end
legend('show', 'location', 'south')
print('e_vs_ar', '-dpng');
