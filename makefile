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

.phony: clean

clean:
	rm -rf $(ODIR)/*.o 
