IDIR=include
cc=gcc
CFLAGS=-I$(IDIR) -fpermissive -fPIC -g
LIBS=-lm
_DEPS=pade.h sara.h poly.h linear.h ddiff.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))
ODIR=obj
_OBJ=pade.o sara.o poly.o linear.o ddiff.o
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
