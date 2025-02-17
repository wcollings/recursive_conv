/*********************************************************************
 * This file created by Synopsys, Inc. for exclusive use with the    *
 * Saber(R) and Saber(R) HDL simulators.                             *
 *              Copyright 1999-2005, Synopsys, Inc.                  *
 * This file and the associated documentation are confidential and   *
 * proprietary to Synopsys, Inc. Your use or disclosure of this file *
 * is subject to the terms and conditions of a written license       *
 * agreement between you, or your company, and Synopsys, Inc.        *
 *********************************************************************/

#ifndef SABERAPI_H
#define SABERAPI_H

#ifdef __cplusplus
#   define AI_EXTERN_C extern "C"
#else 
#   define AI_EXTERN_C
#endif

#if defined(__linux)
#   define C_NAME(x) x
#   define C_CALL(x) C_NAME(x)
#   define F_NAME(x) x##_
#   define F_CALL(x) F_NAME(x)
#   define AI_EXPORT AI_EXTERN_C
#   define AI_IMPORT AI_EXTERN_C
#elif defined(_MSC_VER)
#   define C_NAME(x) x
#   define C_CALL(x) __cdecl C_NAME(x)
#   define F_NAME(x) x
#   define F_CALL(x) F_NAME(x)
#   define AI_EXPORT AI_EXTERN_C __declspec( dllexport )
#   define AI_IMPORT AI_EXTERN_C __declspec( dllimport )
#else
#   error "Unsupported platform"
#endif

/*********************************************************************
 * Foreign routine interface for C routines
 * ========================================
 * The Saber simulator supports one API for foreign routines written in C.
 * The SaberHDL simulator supports three APIs for foreign routines written
 * in C. One corresponds to the API supported by the Saber simulator, the
 * other two are specific to the SaberHDL simulator.
 *
 * API for C routines Common to Saber and SaberHDL
 * -----------------------------------------------
 * The name of the C routine must be in lower case on all Unix systems and 
 * in upper case on the Windows platforms. This API supports foreign 
 * routines for the following purposes:
 * - Saber simulator:
 *   - all foreign routines called for any purpose
 * - SaberHDL simulator:
 *   - foreign routines called from a MAST model or function 
 *   - foreign routines called from a VHDL-AMS model or function whose
 *     corresponding AI_FOREIGN attribute specifies the SABER_F calling 
 *     convention
 * Usage:
 *   To define a foreign routine that in the MAST template or the 
 *   AI_FOREIGN attribute is named my_routine:
 *   On UNIX systems:
 *	#include "saberApi.h"
 *	SABER_FOREIGN_ROUTINE(my_routine)
 *	{ ... }
 *   On Windows systems:
 *	#include "saberApi.h"
 *	SABER_FOREIGN_ROUTINE(MY_ROUTINE)
 *	{ ... }
 */
#define SABER_FOREIGN_ROUTINE(x) \
AI_EXPORT void F_CALL(x) ( \
        double* inp, \
        int*   ninp, \
        int*   ifl, \
        int*   nifl, \
        double* out, \
        int*   nout, \
        int*   ofl, \
        int*   nofl, \
        double* aundef, \
        int*   ier )

/* SaberHDL specific API for C Foreign Routines
 * --------------------------------------------
 * The name of the C routine can be in lower, upper, or mixed case. This
 * API supports foreign routines for the following purpose:
 * - SaberHDL simulator:
 *   - foreign routines called from a VHDL-AMS model or function whose
 *     corresponding AI_FOREIGN attribute specifies the SABER calling 
 *     convention. The name specified in the AI_FOREIGN attribute must
 *     match the name of the C routine exactly.
 * Usage:
 *   To define a foreign routine that in the AI_FOREIGN attribute is named 
 *   my_routine:
 *      #include "saberApi.h"
 *      SABERHDL_C_FOREIGN_ROUTINE(my_routine)
 *      { ... }
 */
#define SABERHDL_C_FOREIGN_ROUTINE(x) \
AI_EXPORT void C_CALL(x) ( \
        double* inp, \
        int*   ninp, \
        int*   ifl, \
        int*   nifl, \
        double* out, \
        int*   nout, \
        int*   ofl, \
        int*   nofl, \
        double* aundef, \
        int*   ier )

/* SaberHDL specific API for Newton Step limiting Functions
 * --------------------------------------------------------
 * The name of the C routine can be in lower, upper, or mixed case. This
 * API supports foreign routines for the following purpose:
 * - SaberHDL simulator:
 *   - foreign routines used for Newton step limiting in a VHDL-AMS model.
 *     The corresponding AI_FOREIGN attribute specifies the NS_LIMIT calling 
 *     convention. The name specified in the AI_FOREIGN attribute must
 *     match the name of the C routine exactly.
 * Usage:
 *   To define a foreign routine that in the AI_FOREIGN attribute is named 
 *   my_routine:
 *      #include "saberApi.h"
 *      SABER_FOREIGN_ROUTINE(my_routine)
 *      { ... }
 */
#define SABERHDL_NSLIMIT_ROUTINE(x) \
AI_EXPORT void C_CALL(x) ( \
        int    numberOfIndeps, \
        int    numberOfDeps, \
        double* prevIndepVector, \
        double* stepIndepVector, \
        double* depVector, \
        double* finalIndepVector, \
        int    numberOfArgs, \
        double* argVector, \
        int*   modifiedIndepFlagsVector, \
        int*   status )

/*********************************************************************
 * Foreign routine interface for Fortran routines
 * ==============================================
 * Such functions must have the following interface:
 * Unix:
 *	subroutine my_routine(inp,ninp,ifl,nifl,out,nout,ofl,nofl,undef,ier)
 *	integer ninp, ifl(*), nifl, nout(2), ofl(*), nofl, ier
 *	real*8 inp(*), out(*), undef
 * Windows:
 *	subroutine my_routine(inp,ninp,ifl,nifl,out,nout,ofl,nofl,undef,ier)
 *	!MS$ATTRIBUTES DLLEXPORT :: my_routine
 *	integer ninp, ifl(*), nifl, nout(2), ofl(*), nofl, ier
 *	real*8 inp(*), out(*), undef
 */

/*********************************************************************
 * Compilation of C foreign routine my_routine in file my_routine.c
 * ----------------------------------------------------------------
 * Linux:
 *	gcc -fPIC -fwritable-strings -Wall -c my_routine.c
 * Windows:
 *	cl /c my_routine.c
 *
 * Creation of a Fortran foreign routine my_routine in file my_routine.f
 * ---------------------------------------------------------------------
 * Linux:
 *	g77 -fPIC -malign-double -fno-second-underscore -c my_routine.f
 * Windows:
 *	fl32 /c my_routine.for
 *   or
 *	df /c my_routine.for
 *
 * Creation of a sharable library
 * ------------------------------
 * Linux:
 *	ld -E -shared -o my_routine.so my_routine.o <libraries>
 * Windows:
 *	link /DLL /OUT:my_routine.dll my_routine.obj <libraries>
 *   or, in one step:
 *	cl /LD my_routine.c <libraries>
 *	fl32 /DLL my_routine.for <libraries>
 *	df /DLL my_routine.for <libraries>
 * where:
 *	<libraries> are any libraries you may need, including:
 *	- linking for the Saber simulator
 *	    on Sun, Linux, IBM
 *		install_home/Saber/bin/libai_foreign.so
 *		install_home/Saber/bin/libai_saber.so
 *		install_home/Saber/bin/libai_analogy.so
 *	    on HP
 *		install_home/Saber/bin/libai_foreign.sl
 *		install_home/Saber/bin/libai_saber.sl
 *		install_home/Saber/bin/libai_analogy.sl
 *	    on Windows
 *		install_home/Saber/lib/libai_foreign.lib
 *		install_home/Saber/lib/libai_saber.lib	
 *		install_home/Saber/lib/libai_analogy.lib
 *	- linking for the SaberHDL simulator
 *	    on Sun
 *		install_home/SaberHDL/bin/libai_TheHDLforeign.so
 *		install_home/SaberHDL/bin/libai_analogy.so
 *		Note: 
 *	    on Windows
 *		install_home/SaberHDL/lib/libai_TheHDLforeign.lib
 *		install_home/SaberHDL/lib/libai_saberhdl.lib	
 *		install_home/SaberHDL/lib/libai_analogy.lib
 *	A library needs to be specified only if the foreign routine calls 
 *	any of the routines below.
 */

/*********************************************************************
 * Service routines that can be called from foreign routines
 * =========================================================
 * To call a service defined using the C_CALL macro, use the C_NAME macro
 * followed by the actual arguments of the function call. Similarly, to
 * call a service defined using the F_CALL macro, use the F_NAME macro
 * followed by the actual arguments of the function call. For example,
 * to get the value of the limited exponential for an actual argument 1.0:
 *	double d = C_NAME(c_limexp) (1.0);
 *
 * String handling routines
 *	requires libai_saber	or libai_saberhdl
 */
AI_IMPORT char * C_CALL(cgetstr) (double);
AI_IMPORT double C_CALL(csetstr) (const char *);

/*
 * Routines supporting statistical analyses
 *	requires libai_foreign	or libai_TheHDLforeign
 */
#if defined(_MSC_VER)
    /* __stdcall routine names must be in upper case on Windows */
#   define cdf CDF
#   define normal NORMAL
#   define pdf PDF
#   define pwl PWL
#   define uniform UNIFORM
#endif
AI_IMPORT void F_CALL(cdf) (double *, int *, int *, int *, 
		double *, int *, int *, int *, double *, int *);
AI_IMPORT void F_CALL(normal) (double *, int *, int *, int *, 
		double *, int *, int *, int *, double *, int *);
AI_IMPORT void F_CALL(pdf) (double *, int *, int *, int *, 
		double *, int *, int *, int *, double *, int *);
AI_IMPORT void F_CALL(pwl) (double *, int *, int *, int *, 
		double *, int *, int *, int *, double *, int *);
AI_IMPORT void F_CALL(uniform) (double *, int *, int *, int *, 
		double *, int *, int *, int *, double *, int *);
/*
 * Routines supporting statistical analyses
 *	requires libai_saber	or libai_saberhdl
 */
AI_IMPORT double C_CALL(c_statsv) (void);
AI_IMPORT double C_CALL(c_wcsv) (void);
AI_IMPORT double C_CALL(c_envrnd) (void);

/*
 * Miscellaneous routines
 *	requires libai_saber	or libai_saberhdl
 */
AI_IMPORT double C_CALL(c_limexp) (double);
AI_IMPORT double C_CALL(c_dsgnnm) (void);

/*
 * Error message handling routines
 *	requires libai_saber	or libai_saberhdl
 */
AI_IMPORT void C_CALL(c_erinst) (int);

/*
 * Error message handling routines
 *	requires libai_analogy
 */
AI_IMPORT void C_CALL(c_erfmt) (const char *, const char *);
AI_IMPORT void C_CALL(c_er) (int);
AI_IMPORT void C_CALL(c_erc) (int, const char *);
AI_IMPORT void C_CALL(c_erca) (int, const char *, int);
AI_IMPORT void C_CALL(c_erd) (int, double);
AI_IMPORT void C_CALL(c_eri) (int, int);
AI_IMPORT void C_CALL(c_ers) (int, const char *);
AI_IMPORT void C_CALL(c_ersys) (int, const char *);
AI_IMPORT void C_CALL(c_erstor) (const char *);
AI_IMPORT void C_CALL(c_ermsg) (int);
AI_IMPORT void C_CALL(c_ermore) (void);
AI_IMPORT void C_CALL(c_erflsh) (void);

/*
 * General output routines
 *	requires libai_analogy
 */
AI_IMPORT int C_CALL(c_ooopen) (const char *);
AI_IMPORT void C_CALL(c_ooclos) (int);
AI_IMPORT int C_CALL(c_oopush) (int);
AI_IMPORT void C_CALL(c_oopop) (void);
AI_IMPORT void C_CALL(c_ooseto) (const char *, int);
AI_IMPORT int C_CALL(c_ooquer) (const char *);
AI_IMPORT int C_CALL(c_oogso) (const char *, int);
AI_IMPORT void C_CALL(c_oob) (int, int);
AI_IMPORT void C_CALL(c_ooca) (const char *, int);
AI_IMPORT void C_CALL(c_oocs) (const char *);
AI_IMPORT void C_CALL(c_oocv) (const char *);
AI_IMPORT void C_CALL(c_ood) (double);
AI_IMPORT void C_CALL(c_ooi) (int);
AI_IMPORT void C_CALL(c_ootab) (int);
AI_IMPORT void C_CALL(c_oonewl) (int);
AI_IMPORT void C_CALL(c_oowrit) (void);
AI_IMPORT void C_CALL(c_ooflsh) (void);
AI_IMPORT void C_CALL(c_oofmt) (const char *);
AI_IMPORT int C_CALL(c_oogetb) (char *, int);
AI_IMPORT void C_CALL(c_oorese) (void);

/*
 * Exit handler support
 *	requires libai_tcl
 */
#ifndef _TCL
typedef void* ClientData;
typedef void (Tcl_ExitProc) (ClientData clientData);
AI_IMPORT void Tcl_CreateExitHandler(Tcl_ExitProc *proc, ClientData clientData);
#endif

#endif /* SABERAPI_H */
