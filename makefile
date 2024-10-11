IDIR=include
cc=gcc
CFLAGS=-I$(IDIR) -g
dllflags=-static-libgcc
LIBS=-lm
_DEPS=pade.h sara.h poly.h linear.h deriv.h interpolate.h central.h saber.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))
ODIR=obj
_OBJ=pade.o sara.o poly.o linear.o deriv.o interpolate.o saber.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))
cc=x86_64-w64-mingw32-gcc


sara_test:main.c $(OBJ)
	$(cc) -o $@ $^ $(CFLAGS) $(LIBS)


$(ODIR)/%.o: src/%.c $(DEPS)
	$(cc) -c -o $@ $< $(CFLAGS)

sara_lib:$(OBJ)
	$(cc) -shared -o libSARA.so $^ $(CFLAGS) $(LIBS)

saber:$(OBJ) src/saber.c
	$(cc) -Wall -c -o obj/saber.o src/saber.c $(LIBS) -I -Iinclude\ -D ADD_EXPORTS
	$(cc) -o ind.dll $(OBJ) -s -shared -Wl,--subsystem,windows -static-libgcc

.phony: clean

clean:
	rm -rf $(ODIR)/*.o 
