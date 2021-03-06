//Parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie,ip,igofr;
double bin_size,nbins;
double stima_prop[m_props];

//Averages
double acc,att;

//Configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

//Thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

//Boolean-like input variables (self-explaining)
int old_config;
int rescale_vel;

//Simulation
int nstep, iprint, seed, blocksize;
const int nblocks=100;
double delta;

//Vectors for blocking average
double sum_prop[m_props];
double ave[m_props][nblocks];
double sum_prog[m_props][nblocks];
double sum2_prog[m_props][nblocks];
double err[m_props][nblocks];

//Functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfOld(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);

void MeasureEachStep(void);
