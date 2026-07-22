// Inputs
double Te = x1;
double w  = x2;

// Time-based derivative
static double w_prev = 0.0;
static double t_prev = 0.0;
static double J_est = 0.0;

double dt = t - t_prev;
double dw_dt = 0.0;

if (dt > 0)
    dw_dt = (w - w_prev) / dt;

// Compute J
if (fabs(dw_dt) > 1e-6)
{
    double J_inst = Te / dw_dt;

    // Smooth J (not derivative)
    J_est = 0.99 * J_est + 0.01 * J_inst;
}

// Update memory
w_prev = w;
t_prev = t;

// Output
y1 = J_est;
