PROJECT = sbmf

BUILDDIR = build

PROJ_SRCS = src/sbmf/sbmf.c

PROJ_OBJS = $(patsubst %.c, $(BUILDDIR)/%.o, $(PROJ_SRCS))

PROJ_LIBS = \
	-lm -fopenmp \
	-l:third_party/lib/libarpack.a \
	-l:third_party/lib/libopenblas.a \
	-l:third_party/lib/libfftw3.a \
	-lgfortran \
	-lgsl

RELEASE_FLAGS = -O3
DEBUG_FLAGS = -g \
			  -fsanitize=address \
			  -fsanitize=leak \
			  #-fsanitize=thread \#
			  -fsanitize=undefined \
			  -fsanitize=bool \
			  -fsanitize=enum \
			  -fsanitize=float-cast-overflow \
			  -fsanitize=signed-integer-overflow

PROJ_FLAGS = -c -fpic -pedantic -Wall -Wextra -Isrc -Ithird_party/include -Iinclude

.PHONY: release
release: PROJ_FLAGS += $(RELEASE_FLAGS)
release: $(BUILDDIR)/$(PROJECT).a

.PHONY: debug
debug: PROJ_FLAGS += $(DEBUG_FLAGS)
debug: $(BUILDDIR)/$(PROJECT).a

.PHONY: install
install:
	mkdir -p ~/ .local/lib
	mkdir -p ~/ .local/include/
	cp $(BUILDDIR)/$(PROJECT).a ~/.local/lib/
	cp -r include/sbmf ~/.local/include/

$(BUILDDIR):
	mkdir -p $(shell find src -type d | sed -e "s/^/$(BUILDDIR)\//")

third_party/lib:
	cd third_party && sh build.sh

$(BUILDDIR)/%.o: %.c
	$(CC) -o $@ $^ -static $(PROJ_FLAGS) $(PROJ_LIBS)

$(BUILDDIR)/$(PROJECT).a: $(BUILDDIR) $(PROJ_OBJS)
	mkdir -p $(BUILDDIR)/tmp
	ar x third_party/lib/libarpack.a --output=$(BUILDDIR)/tmp
	ar x third_party/lib/libopenblas.a --output=$(BUILDDIR)/tmp
	ar x third_party/lib/libfftw3.a --output=$(BUILDDIR)/tmp
	ar rcs $@ $(PROJ_OBJS) $(BUILDDIR)/tmp/*.o
	rm -r $(BUILDDIR)/tmp

.PHONY: clean
clean:
	rm -rf $(BUILDDIR)
