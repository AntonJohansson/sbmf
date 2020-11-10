PROJECT = sbmf

BUILDDIR = build

PROJ_SRCS = \
	src/sbmf/sbmf.c \
	src/sbmf/methods/quadgk_vec.c \
	src/sbmf/methods/item.c \
	src/sbmf/methods/gp2c.c \
	src/sbmf/methods/gp2c_gsl.c \
	src/sbmf/methods/best_meanfield.c \
	src/sbmf/math/find_eigenpairs.c \
	src/sbmf/math/matrix.c \
	src/sbmf/debug/profile.c

PROJ_OBJS = $(patsubst %.c, $(BUILDDIR)/%.o, $(PROJ_SRCS))

TEST_CORRECTNESS_SRCS = src/test/correctness.c
TEST_PERFORMANCE_SRCS = src/test/performance.c

PROJ_LIBS = \
	-lm -fopenmp \
	-l:$(BUILDDIR)/$(PROJECT).a \
	-l:third_party/lib/libarpack.a \
	-l:third_party/lib/libopenblas.a \
	-l:third_party/lib/libfftw3.a \
	-lgfortran \
	-lgsl
TEST_LIBS = \
	-fopenmp \
	-L . \
	-l:$(BUILDDIR)/$(PROJECT).a \
	-l:third_party/lib/libarpack.a \
	-l:third_party/lib/libopenblas.a \
	-l:third_party/lib/libfftw3.a \
	-lgfortran \
	-lcimgui \
	-L/home/aj/.local/lib \
	-lplot \
	-lgsl \
	-lcblas \
	-lm \
	-lpthread

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
TEST_FLAGS = -pedantic -Wall -Wextra -Isrc -Ithird_party/include -I/home/aj/.local/include -Iinclude

all: tests

.PHONY: release
release: PROJ_FLAGS += $(RELEASE_FLAGS)
release: TEST_FLAGS += $(RELEASE_FLAGS)
release: all

.PHONY: debug
debug: PROJ_FLAGS += $(DEBUG_FLAGS)
debug: TEST_FLAGS += $(DEBUG_FLAGS)
debug: all

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

$(BUILDDIR)/test_correctness: $(TEST_CORRECTNESS_SRCS) $(BUILDDIR)/$(PROJECT).a
	$(CC) $(TEST_CORRECTNESS_SRCS) -o $@ $(TEST_FLAGS) $(TEST_LIBS)

$(BUILDDIR)/test_performance: $(TEST_PERFORMANCE_SRCS) $(BUILDDIR)/$(PROJECT).a
	$(CC) $(TEST_PERFORMANCE_SRCS) -o $@ $(TEST_FLAGS) $(TEST_LIBS)


.PHONY: lib
lib: $(BUILDDIR)/$(PROJECT).a

.PHONY: tests
tests: $(BUILDDIR)/test_performance $(BUILDDIR)/test_correctness

.PHONY: clean
clean:
	rm -rf $(BUILDDIR)
