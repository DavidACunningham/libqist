// mqist.h
// header file for fortran interface with matlab in mqist.f90
//
// David Cunningham
// david.cunningham@utexas.edu
// Last revision: see git log
void M_INIT_N(char* namefile_c);

void M_STATE(double* tau, double* state_o);

void M_STM(double* tau, double* stm_o);

void M_STT(double* tau, double* stt_o);

void M_STM_I(double* tau, double* stm_i_o);

void M_STT_I(double* tau, double* stt_i_o);

void M_PROP_ONCE(double* ta, double* tb, double* xa, int* order, double* xb);

void M_STTS_AB(double* taua, double* taub, double* stm_o, double* stt_o);
