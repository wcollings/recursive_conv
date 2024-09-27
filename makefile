IDIR=include
cc=gcc
CFLAGS=-I$(IDIR) -fno-stack-protector -g
dll_cflags=-DADD_EXPORTS -DGW32
dll_eflags=-s -shared -Wl,--subsystem,windows -static-libgcc -DGW32
LIBS=-lm
_DEPS=pade.h sara.h poly.h linear.h deriv.h interpolate.h central.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))
ODIR=obj
_OBJ=pade.o sara.o poly.o linear.o deriv.o interpolate.o saber.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))
wcc=x86_64-w64-mingw32-gcc


sara_test:main.c $(OBJ)
	$(cc) -o $@ $^ $(CFLAGS) $(LIBS)


$(ODIR)/%.o: src/%.c $(DEPS)
	$(cc) -c -o $@ $< $(CFLAGS)

sara_lib:$(OBJ)
	$(cc) -shared -o libSARA.so $^ $(CFLAGS) $(LIBS)

saber:$(OBJ) src/saber.c
	$(wcc) -Wall -c -o obj/saber.o src/saber.c $(LIBS) -I -Iinclude\ $(dllflags) 
	$(wcc) -o ind.dll obj/saber.o $(OBJ) $(dll_eflags) 

print_obj: print_solver.c
	$(cc) print_solver.c obj/sara.o obj/pade.o obj/poly.o obj/linear.o obj/interpolate.o -lm
.phony: clean

clean:
	rm -rf $(ODIR)/*.o 
