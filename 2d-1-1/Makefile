# $Id: Makefile,v 1.1 2009/12/10 17:34:47 pozdneev Exp $

.SUFFIXES:
.SUFFIXES: .c .h .o .d

CC = gcc

SILO_DIR = /opt/silo/4.8-bsd

FLAGS = \
	-Wall \
	-Wextra \
	-O0 \

CFLAGS = \
	$(FLAGS) \
	-I$(SILO_DIR)/include \
	-D_FILE_OFFSET_BITS=64 \
	-D_LARGEFILE64_SOURCE \

LFLAGS = \
	$(FLAGS) \
	-L$(SILO_DIR)/lib \

LIBS = \
       -lm \
       -lsilo \
       -lefence \

BIN = \
       2d \

OBJS = \
       main.o

.PHONY: all
all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(LFLAGS) $^ $(LIBS) -o $@

-include $(OBJS:.o=.d)

%.o: %.c
	$(CC) $(CFLAGS) -c $<
	$(CC) $(CFLAGS) -c -MM $< >$*.d

tags: *.c
	gnuctags *.c

.PHONY: clean
clean:
	rm -f *.o $(BIN) *.d

.PHONY: distclean
distclean: clean
	rm -f *~ tags *.bin *.silo
