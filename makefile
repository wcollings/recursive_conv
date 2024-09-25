IDIR=include
cc=gcc
CFLAGS=-I$(IDIR) -fPIC -g
LIBS=-lm
_DEPS=pade.h sara.h poly.h linear.h deriv.h interpolate.h central.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))
ODIR=obj
_OBJ=pade.o sara.o poly.o linear.o deriv.o interpolate.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

sara_test:main.c $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)


$(ODIR)/%.o: src/%.c $(DEPS)
	$(cc) -c -o $@ $< $(CFLAGS)

sara_lib:$(OBJ)
	gcc -shared -o libSARA.so $^ $(CFLAGS) $(LIBS)

saber:$(OBJ) saber.c
	gcc -m32 -fPIC -Wall -c saber.c -o ind.dll -I/usr/include/saber -Iinclude/
	# x86_64-w64-mingw32-gcc -m32 -fPIC -Wall -c saber.c -o ind.dll -I/usr/include/saber -Iinclude/

tcl:
	gcc -shared -o libhello.so -DUSE_TCL_STUBS -L/usr/share/tcltk/tcl8.6/ -ltclstub8.6 -I/usr/include/tcl8 -fPIC hello.c

.phony: clean

clean:
	rm -rf $(ODIR)/*.o 
