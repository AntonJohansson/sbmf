CC = gcc
PROJECT = sbmf

BUILDDIR = build

PROJ_SRCS = \
	src/sbmf/sbmf.c \
	src/sbmf/methods/quadgk_vec.c \
	src/sbmf/methods/quadgk_vec_inl.c \
	src/sbmf/methods/item.c \
	src/sbmf/methods/scim.c \
	src/sbmf/methods/gp2c.c \
	src/sbmf/math/find_eigenpairs.c \
	src/sbmf/math/matrix.c \
	src/sbmf/debug/profile.c \
	src/sbmf/debug/log.c \
	src/sbmf/threading/thread_pool.c

PROJ_OBJS = $(patsubst %.c, $(BUILDDIR)/%.o, $(PROJ_SRCS))

TEST_CORRECTNESS_SRCS = src/test/correctness.c
TEST_PERFORMANCE_SRCS = src/test/performance.c

PROJ_LIBS = \
	-lm -fopenmp \
	-l:$(BUILDDIR)/$(PROJECT).a \
	-l:third_party/lib/libarpack.a \
	-l:third_party/lib/libopenblas.a \
	-l:third_party/lib/libfftw3.a \
	-lgfortran
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

MODE_FLAGS = $(RELEASE_FLAGS)
PROJ_FLAGS = $(MODE_FLAGS) -c -fpic -pedantic -Wall -Wextra -Isrc -Ithird_party/include
TEST_FLAGS = $(MODE_FLAGS)          -pedantic -Wall -Wextra -Isrc -Ithird_party/include -I/home/aj/.local/include

all: tests

$(BUILDDIR):
	mkdir -p $(shell find src -type d | sed -e "s/^/$(BUILDDIR)\//")

third_party/lib:
	cd third_party && sh build.sh

$(BUILDDIR)/%.o: %.c
	$(CC) -o $@ $^ -static $(PROJ_FLAGS) $(PROJ_LIBS)

$(BUILDDIR)/$(PROJECT).a: $(BUILDDIR) $(PROJ_OBJS)
	ar rcs $@ $(PROJ_OBJS)

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
