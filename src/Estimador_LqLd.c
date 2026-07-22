static double iq_prev = 0.0;
static double diq_dt_f = 0.0;
static double Lq_avg = 0.0;
static int valid_count = 0;

double vq3, iq3, id3, we3, Rs3, en3;
double diq_dt;
double Lq_raw;
double valid;

// thresholds
double eps_diq = 1e-3;
double eps_id  = 0.2;     // d-axis current should be near zero
double eps_we  = 1e-2;
double eps_vq  = 1.0;     // minimum applied vq
double iq_max_est = 3.0;  // limit to early transient

// Inputs
vq3 = x1;
iq3 = x2;
id3 = x3;
we3 = x4;
Rs3 = x5;
en3 = x6;

// raw derivative
diq_dt = (iq3 - iq_prev) / delt;

// simple low-pass filter on derivative
diq_dt_f = diq_dt_f + 0.1 * (diq_dt - diq_dt_f);

// default outputs
Lq_raw = 0.0;
valid = 0.0;

// Valid only during pulse transient and locked-rotor condition
if (en3 > 0.5) {
    if ((fabs(vq3) > eps_vq) &&
        (fabs(diq_dt_f) > eps_diq) &&
        (fabs(id3) < eps_id) &&
        (fabs(we3) < eps_we) &&
        (fabs(iq3) < iq_max_est)) {

        Lq_raw = (vq3 - Rs3 * iq3) / diq_dt_f;

        // reject non-physical values (adjust bounds if needed)
        if ((Lq_raw > 0.001) && (Lq_raw < 0.2)) {
            valid = 1.0;

            valid_count++;
            if (valid_count == 1) {
                Lq_avg = Lq_raw;
            } else {
                Lq_avg = Lq_avg + (Lq_raw - Lq_avg) / valid_count;
            }
        }
    }
}

// Save previous current
iq_prev = iq3;

// Outputs
y1 = diq_dt_f;
y2 = Lq_raw;
y3 = Lq_avg;
y4 = valid;
