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
