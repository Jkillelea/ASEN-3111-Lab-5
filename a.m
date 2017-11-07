function a0 = a(t, a0_r, a0_t)
  a0 = a0_r + (a0_t - a0_r)*cos(t);
end
