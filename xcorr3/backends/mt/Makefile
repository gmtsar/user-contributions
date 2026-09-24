# Makefile for xcorr_mt: multi-process (fork) version of GMTSAR xcorr.
# Link against the same numerical libraries as the reference GMTSAR build.
# See README for the limits of byte-level compatibility.

GMTSAR_HOME ?= /usr/local/GMTSAR
PREFIX      ?= $(HOME)/.local

CC          = gcc
CFLAGS      = -O2 -Wall -fPIC -fno-strict-aliasing -std=c99
INCLUDES    = -I/usr/include/gmt -I$(GMTSAR_HOME)/gmtsar -I/usr/local/include
# -z muldefs: gmtsar.h defines verbose/debug in every object (tentative defs)
LDFLAGS     = -z muldefs -s -Wl,-rpath,/usr/lib/x86_64-linux-gnu
LIBS        = -L$(GMTSAR_HOME)/gmtsar -lgmtsar \
              -L/usr/lib/x86_64-linux-gnu -lgmt \
              -llapack -lblas -lm -L/usr/local/lib -ltiff -lm

xcorr_mt: xcorr_mt.o
	$(CC) $(LDFLAGS) xcorr_mt.o $(LIBS) -o $@

xcorr_mt.o: xcorr_mt.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f xcorr_mt.o xcorr_mt

install: xcorr_mt
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 755 xcorr_mt "$(DESTDIR)$(PREFIX)/bin/xcorr_mt"

.PHONY: clean install
