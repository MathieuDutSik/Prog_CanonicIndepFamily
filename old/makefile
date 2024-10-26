CFILE = $(wildcard *.c)
OFILE = $(wildcard *.o)
PROG = $(patsubst %.c,%.prog, $(CFILE))
OBJS = $(patsubst %.c,%.o, $(CFILE))
CFLAGS = -Wall -O2

GMP_PREFIX = /opt/gmp-4.1.4_V2

%.o: %.c ; gcc $(CFLAGS) -I$(GMP_PREFIX)/include $< -c -o $@

%.prog: %.o ; gcc $< -L$(GMP_PREFIX)/lib -lgmp -o $@


all: $(PROG) $(OBJS)

clean:
	rm -f $(PROG) $(OBJS)
