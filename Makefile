
EDITOR = /usr/bin/vim
INDENT = /usr/bin/indent

CC        = clang
LDFLAGS   = -O
CFLAGS    = -Weverything -Wextra -pedantic $(LDFLAGS)
LDLIBS    = -lm -lgsl -lgslcblas
#LDLIBS    = $(shell gsl-config --libs)

.SUFFIXES:
.SUFFIXES:  .c .o .h

.PHONY: edit clean veryclean

target    = metropolis_calculator

$(target): $(target).c

edit : $(target).c
	$(EDITOR) $<
	$(INDENT) $<

clean : 
	rm -f *.o
	rm -f *.*~
	rm  -f $(target)
