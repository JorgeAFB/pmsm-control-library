// Inputs
double Te    = x1;
double w     = x2;
double Tload = x3;

// Memory
static double w_prev = 0.0;
static double Te_prev = 0.0;
static double t_prev = 0.0;
static double B_est = 0.0;

// Time
double dt = t - t_prev;

// Derivative
double dw_dt = 0.0;
if (dt > 0)
    dw_dt = (w - w_prev) / dt;

// --- Stability conditions ---
int steady_speed  = fabs(dw_dt) < 5e-2;
int steady_torque = fabs(Te - Te_prev) < 0.05;
int valid_speed   = fabs(w) > 1.0;   // avoid low-speed noise

if (steady_speed && steady_torque && valid_speed)
{
    double B_inst = (Te - Tload) / w;

    // Reject unrealistic values
    if (B_inst > 0 && B_inst < 1.0)
    {
        B_est = 0.99 * B_est + 0.01 * B_inst;
    }
}

// Update memory
w_prev = w;
Te_prev = Te;
t_prev = t;

// Output
y1 = B_est;
