function [x_a, y_a] = NACA_Airfoil(m, p, t, c, N)
  % Creates a NACA 4-digit series airfoil based on the parameters
  % m -> maximum camber
  % p -> location of max camber
  % t -> thickness
  % c -> chord length
  % N -> number of panels to be used

  m = m/100;
  p = p/10;
  t = t/100;

  % x = linspace(0, c, N);
  dx = c/(N/2);
  x = 0:dx:c;

  % coefficients
  c1 =  0.2969;
  c2 = -0.1260;
  c3 = -0.3516;
  c4 =  0.2843;
  c5 = -0.1036;

  % half thickness (from the mean camber line)
  yt = (t/0.2)*c * (c1.*sqrt(x./c) + c2.*(x./c)    ...
                  + c3.*(x./c).^2  + c4.*(x./c).^3 ...
                  + c5.*(x./c).^4);

  [y_c, y_c_p] = yc(x, p, c, m);

  zeta = atan(y_c_p);

  % upper and lower surfaces
  xu = (x - yt.*sin(zeta))';
  xl = (x + yt.*sin(zeta))';

  yu = (y_c + yt.*cos(zeta))';
  yl = (y_c - yt.*cos(zeta))';

  % format and return
  x_a = [flip(xl); xu(2:end)];
  y_a = [flip(yl); yu(2:end)];
end

% mean camber line
function [y_c, y_c_p] = yc(x, p, c, m)
  if (0 <= x) & (x <= p*c)
    y_c = (m/p^2).*x.*(2*p - x./c);
    y_c_p = (m*(2*p - x/c))/p^2 - (m*x)/(c*p^2);
  else
    y_c = (m/(1-p)^2).*(c - x).*(1 + x./c - 2*p);
    y_c_p = (m*(c - x))/(c*(p - 1)^2) - (m*(x/c - 2*p + 1))/(p - 1)^2;
  end
end
