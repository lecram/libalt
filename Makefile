CFLAGS = -Wall -Wextra -Werror -std=c99
DEBUG ?= 0
ifeq ($(DEBUG), 1)
    CFLAGS += -O0 -g
else
    CFLAGS += -O2
endif

all: libalt.so

libalt.so: alt.o
	gcc -shared -o $@ $<

alt.o: alt.c
	gcc ${CFLAGS} -c -fpic -o $@ $< -lm

.PHONY: clean
clean:
	rm -f alt.o
