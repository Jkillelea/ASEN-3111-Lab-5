function [m, p, t] = NACA_from_4_digit(num)
  % Ex: for a NACA 4412, t = 12, p = 4, m = 4
  t = mod(num, 100);
  p = mod(round(num/100), 10);
  m = mod(round(num/1000), 10);
end
