function [gamma, cp, cl] = Vortex_Panel(xb, yb, alpha)
  m = length(xb)-1;
  x      = zeros(1, m);
  y      = zeros(1, m);
  s      = zeros(1, m);
  sine   = zeros(1, m);
  cosine = zeros(1, m);
  theta  = zeros(1, m);
  v      = zeros(1, m);
  cp     = zeros(1, m);
  gamma  = zeros(1, m);
  rhs    = zeros(1, m);
  cn1    = zeros(m);
  cn2    = zeros(m);
  ct1    = zeros(m);
  ct2    = zeros(m);
  an     = zeros(m+1);
  at     = zeros(m, m+1);

  mp1 = m + 1;

  for i = 1:m
    ip1 = i + 1;
    x(i)      = 0.5*(xb(i) + xb(ip1));
    y(i)      = 0.5*(yb(i) + yb(ip1));
    s(i)      = sqrt( (xb(ip1) - xb(i))^2 + (yb(ip1) - yb(i))^2 );
    theta(i)  = atan2( (yb(ip1) - yb(i)), (xb(ip1) - xb(i)) );
    sine(i)   = sin(theta(i));
    cosine(i) = cos(theta(i));
    rhs(i)    = sin(theta(i) - alpha);
  end

  for i = 1:m
    for j = 1:m
      if i == j
        cn1(i, j) = -1.0;
        cn2(i, j) = 1.0;
        ct1(i, j) = 0.5*pi;
        ct2(i, j) = 0.5*pi;
      else
        a = -(x(i) - xb(j))*cosine(j) - (y(i) - yb(j))*sine(j);
        b = (x(i) - xb(j))^2 + (y(i) - yb(j))^2;
        c = sin(theta(i) - theta(j));
        d = cos(theta(i) - theta(j));
        e = (x(i) - xb(j))*sine(j) - (y(i) - yb(j))*cosine(j);
        f = log(1 + (s(j)^2 + 2*a*s(j))/b );
        g = atan2(e*s(j), b+a*s(j));

        p = (x(i) - xb(j)) * sin(theta(i)-2.*theta(j)) ...
				  + (y(i) - yb(j)) * cos(theta(i)-2.*theta(j));

        q = (x(i) - xb(j)) * cos(theta(i)-2.*theta(j)) ...
          - (y(i) - yb(j)) * sin(theta(i)-2.*theta(j));

        cn2(i, j) = d + 0.5*q*f/s(j) - (a*c + d*e)*g/s(j);
        cn1(i, j) = 0.5*d*f + c*g - cn2(i, j);
        ct2(i, j) = c + 0.5*p*f/s(j) + (a*d - c*e)*g/s(j);
        ct1(i, j) = 0.5*c*f - d*g - ct2(i, j);
      end % if
    end % loop j
  end % loop i

	for i = 1:m
		an(i, 1)   = cn1(i, 1);
		an(i, mp1) = cn2(i, m);
		at(i, 1)   = ct1(i, 1);
		at(i, mp1) = ct2(i, m);

		for j = 2:m
      an(i, j) = cn1(i, j) + cn2(i, j-1);
      at(i, j) = ct1(i, j) + ct2(i, j-1);
		end
	end

	an(mp1, 1)   = 1.0;
	an(mp1, mp1) = 1.0;

	for j = 2:m
		an(mp1, j) = 0.0;
	end

	rhs(mp1) = 0.0;

  [~, ~, gamma, ~] = cramer(an, rhs, gamma, mp1, m);

  for i = 1:m
    v(i) = cos(theta(i) - alpha);
    for j = 1:mp1
      v(i) = v(i) + at(i, j)*gamma(j);
      cp(i) = 1 - v(i)^2;
    end
  end

  gamma_net = sum(v.*s);
  cl = 2*gamma_net;

end % function airfoil

function [c, a, x, n] = cramer(c, a, x, n, m)
  cc    = zeros(m+1);
  denom = det(c); % call MATLAB builtin for speed instead of the translated code

  for k = 1:n
    for i = 1:n
      for j = 1:n
        cc(i, j) = c(i, j);
      end % j
    end % i

    for i = 1:n
      cc(i, k) = a(i);
    end % i
    x(k) = det(cc) / denom;
  end % k
end % function cramer
