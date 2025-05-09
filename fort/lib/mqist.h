// mqist.h
// header file for fortran interface with matlab in mqist.f90
//
// David Cunningham
// david.cunningham@utexas.edu
// Last revision: see git log

<<<<<<< HEAD
void m_init_n(char* namefile_c);

void m_state(double* tau, double* state_o);

void m_stm(double* tau, double* stm_o);

void m_stt(double* tau, double* stt_o);

void m_stm_i(double* tau, double* stm_i_o);

void m_stt_i(double* tau, double* stt_i_o);

void m_prop_once(double* ta, double* tb, double* xa, int* order, double* xb);

void m_stts_ab(double* taua, double* taub, double* stm_o, double* stt_o);
