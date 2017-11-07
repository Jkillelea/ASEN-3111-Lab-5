function a = aero(t, aero_r, aero_t)
  a = aero_r + (aero_t - aero_r)*cos(t);
end
