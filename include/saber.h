#ifdef ADD_EXPORTS
  #define ADDAPI __declspec(dllexport)
#else
  #define ADDAPI __declspec(dllimport)
#endif

/* Define calling convention in one place, for convenience. */
#define ADDCALL __cdecl

/* Make sure functions are exported with C linkage under C++ compilers. */

#ifdef __cplusplus
extern "C"
{
#endif

#include <saberApi.h>

enum call_tp { INIT,STEP,ACCEPT,START,END };

/* Declare our Add function using the above definitions. */
void IND(double* inp,int* ninp,int* ifl,int* nifl,double* out,int* nout,int *ofl,int *nofl,double *aundef,int *ier);
/* ADDAPI void ADDCALL IND(double* inp,int* ninp,int* ifl,int* nifl,double* out,int* nout,int *ofl,int *nofl,double *aundef,int *ier); */
struct Solver_t * do_setup(double * in);

#ifdef __cplusplus
} // __cplusplus defined.
#endif
