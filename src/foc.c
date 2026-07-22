static double wMc_rpm_last;
static double wMc_acc_pu;

double wMc_rpm;
double vab, vbc, va, vb, vc, valpha, vbeta, vd, vq;

double ia, ib, ic,  ialpha, ibeta, id, iq;

double vc_d, vc_q;
double vc_alpha, vc_beta;
double vc_a, vc_b, vc_c;
double vdc;          
double ma_tst, mb_tst, mc_tst;

double thetaMc_pu;
double wm;
double thetam_pu;
double thetam;
double psi_f;
double Tem_tst1;
double pelec;

// --- Voltage limiting / anti-windup ---
double vc_d_raw, vc_q_raw;
double Vmag, Vmax, S;
double dv_d, dv_q;
double kaw = 400.0;       // using PSIM value
double eps = 1e-6;        // like 1u block 

// --- Static variables (Integrators must be static to persist) ---
static double xd_i = 0.0; 
static double xq_i = 0.0;

// variables for PI Control
double id_ref, iq_ref, ed, eq, kp, ki;
double vPI_d, vPI_q;

// Constants (matching  motor parameters)
double m_Rs_tst = 1.5; 
double m_Ls_tst = 0.060;

double tmp;
double tmpcos;
double tmpsin;

#define pi 3.1415926535897932385
#define _2pi 6.2831853071795864770

//------------------INPUTS-----------------------------------
wMc_rpm = x1;
vab = x2;
vbc = x3;
ia = x4;
ib = x5;
psi_f = x6;
vdc     = x7;
id_ref  = x8;
iq_ref  = x9;

//-------------------------------------------------------------------
wMc_acc_pu = wMc_acc_pu + (wMc_rpm+wMc_rpm_last) * delt / (2*60.0);
         if ( wMc_acc_pu > 1 ) wMc_acc_pu = wMc_acc_pu - 1;
else if ( wMc_acc_pu < 0 ) wMc_acc_pu = wMc_acc_pu + 1;
wMc_rpm_last = wMc_rpm;

thetaMc_pu = wMc_acc_pu;

//-------------------------------------------------------------------
wm = wMc_rpm * _2pi / 60 * m_p;
tmp = thetaMc_pu*m_p;

//thetam_pu = tmp - ceil(tmp) + 1;
thetam_pu = tmp - floor(tmp);   // gives [0,1)

thetam = thetam_pu * _2pi;

//-------------------------------------------------------------------
va = (2*vab +   vbc)/3;
vb = (  -vab +   vbc)/3;
vc = (  -vab - 2*vbc)/3;
ic = -(ia + ib);

tmpcos = cos(thetam);
tmpsin = sin(thetam);

//valpha = sqrt(6)/2*va;
//vbeta = (sqrt(2)/2)*va + sqrt(2)*vb;
//vd = tmpcos*valpha + tmpsin*vbeta;
//vq = -tmpsin*valpha +tmpcos*vbeta;

ialpha = sqrt(6)/2 * ia;
ibeta  = (sqrt(2)/2)*ia + sqrt(2)*ib;
id =  tmpcos*ialpha + tmpsin*ibeta;
iq = -tmpsin*ialpha + tmpcos*ibeta;

// ------------- PI GAIN CALCULATION --------------
{
    double Tstl = 0.005; // 5ms
    double Mp   = 1.0;   // 1%
    double eW   = 1.0;   // 1%
    
    double lnMp = log(Mp/100.0);
    double zeta = -lnMp / sqrt(lnMp*lnMp + pi*pi);
    double lnEW = log(eW/100.0);
    double wn   = -lnEW / (Tstl * zeta);

    ki = wn * wn * m_Ls_tst;
    kp = 2.0 * zeta * wn * m_Ls_tst - m_Rs_tst;
}
// ------------------ CURRENT CONTROL (PI + decoupling + limiter + AW) ---------------------
// Error signals
ed = id_ref - id;
eq = iq_ref - iq;

// PI output using previous integrator state
vPI_d = (kp * ed) + xd_i;
vPI_q = (kp * eq) + xq_i;

// Raw dq voltage request (PI + decoupling/feedforward)
vc_d_raw = vPI_d - (wm * m_Ls_tst * iq);
vc_q_raw = vPI_q + (wm * (id * m_Ls_tst + psi_f));

// Voltage vector limiter (SVPWM linear region)
Vmag = sqrt(vc_d_raw*vc_d_raw + vc_q_raw*vc_q_raw);
Vmax = 0.95 * vdc * 0.5773502691896258; // 1/sqrt(3)

if (Vmag > eps) {
    S = Vmax / (Vmag + eps);
    if (S > 1.0) S = 1.0;
    if (S < 0.0) S = 0.0;
} else {
    S = 1.0;
}

vc_d = S * vc_d_raw;
vc_q = S * vc_q_raw;

// Back-calculation anti-windup (inject limiter error)
dv_d = vc_d - vc_d_raw;
dv_q = vc_q - vc_q_raw;

xd_i += (ki * ed + kaw * dv_d) * delt;
xq_i += (ki * eq + kaw * dv_q) * delt;

// ------------------ DECOUPLING (Per PSIM Schematic) ----------
// D-axis: vc_d = vPI_d + (w_e * -m_Ls * iq_f)
//vc_d = vPI_d - (wm * m_Ls_tst * iq);

// Q-axis: vc_q = vPI_q + (w_e * (id_f * m_Ls + psi_f))
//vc_q = vPI_q + (wm * (id * m_Ls_tst + psi_f));
//-------------------------------------------------------------

Tem_tst1 = m_p * psi_f * iq;
pelec = va*ia + vb*ib + vc*ic;   // instantaneous 3-phase power

//-------------------------------------------------------------
// --- Command path replication (dq -> alpha/beta -> abc) ---
vc_alpha =  tmpcos*vc_d - tmpsin*vc_q;
vc_beta  =  tmpsin*vc_d + tmpcos*vc_q;

// Inverse of custom Clarke mapping
vc_a = (2.0/sqrt(6.0)) * vc_alpha;
vc_b = (1.0/sqrt(2.0)) * vc_beta - 0.5 * vc_a;
vc_c = -(vc_a + vc_b);

// Optional: linear utilization estimate (NOT full SVPWM)
if (vdc > 1e-6) {
  ma_tst = 2.0 * vc_a / vdc;
  mb_tst = 2.0 * vc_b / vdc;
  mc_tst = 2.0 * vc_c / vdc;
} else {
  ma_tst = mb_tst = mc_tst = 0.0;
}
//-------------------OUTPUTS----------------------------
y1 = wMc_acc_pu;
y2 = wm;
y3 = thetam_pu;
y4 = thetam;

y5 = va;
y6 = vb;
y7 = vc;
y8 = valpha;
y9 = vbeta;
y10 = vd;
y11 = vq;
y12 = ic;
y13 = ialpha;
y14 = ibeta;
y15 = id;
y16 = iq;
y17 = Tem_tst1;
y18 = pelec;

y19 = vc_alpha;
y20 = vc_beta;
y21 = vc_a;
y22 = vc_b;
y23 = vc_c;
y24 = ma_tst;
y25 = mb_tst;
y26 = mc_tst;
