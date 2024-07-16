IDIR=include
cc=gcc
CFLAGS=-I$(IDIR) -fPIC -g
LIBS=-lm
_DEPS=pade.h sara.h poly.h linear.h deriv.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))
ODIR=obj
_OBJ=pade.o sara.o poly.o linear.o deriv.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: src/%.c $(DEPS)
	$(cc) -c -o $@ $< $(CFLAGS)

sara_test:$(OBJ)
	gcc -o $@ main.c $^ $(CFLAGS) $(LIBS)

sara_lib:$(OBJ)
	gcc -shared -o libSARA.so $^ $(CFLAGS) $(LIBS)

.phony: clean

clean:
	rm -rf $(ODIR)/*.o 
