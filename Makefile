.SUFFIXES : .cpp .o

CC = gcc
LD = gcc
INCDIR = stk
CFLAGS = -std=c++11 -I$(INCDIR) -O3
LFLAGS = -lm
LIBS =  -lstdc++
PROGRAM  = dpcmc

IFILES = stk/FileRead.h \
		 stk/FileWrite.h \
		 stk/Stk.h \
				 
CFILES = main.cpp \
		 fft4g.c \
		 stk/FileRead.cpp \
		 stk/FileWrite.cpp \
		 stk/Stk.cpp \

OFILES = $(CFILES:.cpp=.o)

all: $(PROGRAM)

$(PROGRAM):		$(OFILES)
				$(LD) $(LFLAGS) $? -o bin/$(PROGRAM) $(LIBS)

$(OFILES):	$(CFILES) $(IFILES) Makefile

.cpp.o : $(CFILES) $(IFILES)
	$(CC) $(CFLAGS) -c -o $*.o $*.cpp

clean:
	rm -f $(OFILES) $(PROGRAM)
