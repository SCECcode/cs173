# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

TARGETS = libcs173.a libcs173.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libcs173.so ${prefix}/lib
	cp libcs173.a ${prefix}/lib
	cp cs173.h ${prefix}/include

clean:
	rm -rf $(TARGETS)
	rm -rf *.o

libcs173.a: cs173_static.o cs173_gtl_static.o
	$(AR) rcs $@ $^

libcs173.so: cs173.o cs173_gtl.o
	$(CC) -shared $(AM_FCFLAGS) -o libcs173.so $^ $(AM_LDFLAGS)

cs173.o: cs173.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)

cs173_gtl.o: cs173_gtl.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
cs173_static.o: cs173.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

cs173_gtl_static.o: cs173_gtl.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	
