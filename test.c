#include <tcl8.6/tcl.h>

static int hello_cmd(ClientData cdata, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
	Tcl_SetObjResult(interp,Tcl_NewStringObj("Hello World",-1));
	return TCL_OK;
}

int DLLEXPORT hello_init(Tcl_Interp * interp) {
	if (
	#ifdef USE_TCL_STUBS
			Tcl_InitStubs(interp,TCL_VERSION,0)
	#else
		Tcl_PkgRequire(interp,"Tcl","8.6",0)
	#endif
		==NULL) {
		return TCL_ERROR;
	}
	/* changed this to check for an error - GPS */
     if (Tcl_PkgProvide(interp, "Hello", "1.0") == TCL_ERROR) {
         return TCL_ERROR;
     }
     Tcl_CreateObjCommand(interp, "hello", hello_cmd, NULL, NULL);
     return TCL_OK;
}
