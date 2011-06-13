## A basic makefile

CC = gcc
CFLAGS = -Wall
LDFLAGS = -lgsl -lgslcblas

PROG = curvefit
HDRS = fit_functions.h
SRCS = main.c fit_functions.c

OBJS = $(SRCS:.c=.o)

$(PROG) : $(OBJS)
	$(CC) $(LDFLAGS) -o $(PROG) $(OBJS)
         
main.o : main.c fit_functions.h

fit_functions.o : fit_functions.c fit_functions.h

.PHONY : clean
clean :
	-rm $(PROG) $(OBJS)
