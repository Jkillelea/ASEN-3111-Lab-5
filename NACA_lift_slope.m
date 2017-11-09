function [dcl_da, alpha_L0] = NACA_lift_slope(airfoil_num, c)
  % returns lift slope, units 1/degrees and zero-lift angle of attack
  N = 50;

  [m, p, t]   = NACA_from_4_digit(airfoil_num);
  alpha_range = deg2rad(-5:1:10);
  cl_vals     = zeros(1, length(alpha_range));
  [x, y]      = NACA_Airfoil(m, p, t, c, N);

  for i = 1:length(alpha_range)
    alpha  = alpha_range(i);
    [~, ~, cl] = Vortex_Panel(x, y, alpha);
    cl_vals(i) = cl;
  end

  a        = polyfit(alpha_range, cl_vals, 1); % linear fit
  dcl_da   = a(1); % cl slope
  alpha_L0 = -a(2)/a(1);
end
