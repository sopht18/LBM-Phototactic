// const static int LX = 128 +2;
// const static int LY = 128 +2;
const static int LX = 128 +2;
const static int LY = 64 +2;
const static int dat_num = 1000;
const static int duration = 400;
const static int T = duration*dat_num; // step number
const static double a = -0.01;
const static double b = -a;
const static double kappa = -1.0*a;
const static double Gamma = 1.0; // mobility
const static double dt = 1.0;
const static double tau = 0.65; // viscosity : eta = (tau -0.5)/3
const static double grav = 0.00004;
const static double dx = 1.0;
const static double dy = dx;
const static double csq = dx*dx/(dt*dt);

const static double PI = 3.141592653589;
const static double phi_in = 1.0;
const static double zeta = -0.5;
const static double lmd  = 8.0; // damping length
const static double v0ph = 0.003; // constant for phototactic term
