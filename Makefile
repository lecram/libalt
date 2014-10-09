CFLAGS+= -Wall -Wextra -Werror -std=c99 -O2

all: libalt.so

libalt.so: alt.o
	gcc -shared -o $@ $<

alt.o: alt.c
	gcc ${CFLAGS} -c -o $@ $<

.PHONY: clean
clean:
	rm alt.o
