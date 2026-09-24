PREFIX ?= $(HOME)/.local
GMTSAR_HOME ?= /usr/local/GMTSAR

# Build only the backend you need; installation never replaces stock xcorr.
all:
	@echo 'Choose: make mt, make cc, or make cc-system'

mt:
	$(MAKE) -C backends/mt GMTSAR_HOME="$(GMTSAR_HOME)"

cc:
	$(MAKE) -C backends/cc

cc-system:
	cd backends/cc && bash build_local.sh

install:
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 755 xcorr3 "$(DESTDIR)$(PREFIX)/bin/xcorr3"

install-mt: mt install
	install -m 755 backends/mt/xcorr_mt "$(DESTDIR)$(PREFIX)/bin/xcorr_mt"

# Use after either cc build variant, so system-CUDA settings are not lost.
install-cc: install
	install -m 755 backends/cc/xcorr_cc "$(DESTDIR)$(PREFIX)/bin/xcorr_cc"

test:
	python3 test/test_xcorr3.py

.PHONY: all mt cc cc-system install install-mt install-cc test
