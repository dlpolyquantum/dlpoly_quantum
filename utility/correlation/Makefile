#.PHONY all clean

PROG = correlation
FC = mpif90
PATHFLAGS = -I/opt/homebrew/Cellar/fftw/3.3.10_1/include
LIBFLAGS = -L/opt/homebrew/Cellar/fftw/3.3.10_1/lib -lfftw3
OBJS = corr.o
$(PROG) : $(OBJS)
	$(FC) -o $(PROG) $(OBJS)  $(PATHFLAGS) $(LIBFLAGS)

corr.o : corr.f90
	$(FC) -g $(PATHFLAGS) $(LIBFLAGS) -c corr.f90

debug: $(OBJS)
	$(FC) -o $(PROG) -g -Wall -fcheck=all $(OBJS)  $(PATHFLAGS) $(LIBFLAGS)

clean:
	rm -f core $(PROG) $(OBJS)
