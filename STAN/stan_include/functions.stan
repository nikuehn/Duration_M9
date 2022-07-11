functions {
  real logistic_hinge(real x, real x0, real a, real b0, real b1, real delta) { 
    real xdiff = x - x0;
    return a + b0 * xdiff + (b1 - b0) * delta * log1p_exp(xdiff / delta);
  }
}
