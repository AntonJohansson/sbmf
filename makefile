CC = gcc
PROJECT = sbmf

BUILDDIR = build

PROJ_SRCS = \
	src/sbmf/sbmf.c \
	src/sbmf/methods/quadgk_vec.c \
	src/sbmf/methods/quadgk_vec_inl.c \
	src/sbmf/methods/item.c \
	src/sbmf/methods/scim.c \
	src/sbmf/math/find_eigenpairs.c \
	src/sbmf/math/matrix.c \
	src/sbmf/debug/profile.c \
	src/sbmf/debug/log.c \
	src/sbmf/threading/thread_pool.c

PROJ_OBJS = $(patsubst %.c, $(BUILDDIR)/%.o, $(PROJ_SRCS))

TEST_CORRECTNESS_SRCS = src/test/correctness.c
TEST_PERFORMANCE_SRCS = src/test/performance.c

PROJ_LIBS = \
	-lm -fopenmp -lgsl
TEST_LIBS = \
	-fopenmp \
	-lgsl \
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

#PROJ_FLAGS = -fsanitize=address -fsanitize=leak -g -fpic -Wall -Werror -Isrc
#TEST_FLAGS = -fsanitize=address -fsanitize=leak -g -Wall -Werror -Isrc -Isrc/test/plotting/third_party

#PROJ_FLAGS = -c -pg -g -fpic -O0 -Wall -Isrc -Ithird_party/include
#TEST_FLAGS = -pg -g -Isrc -Ithird_party/include -I/home/aj/.local/include

RELEASE_FLAGS = -O3
DEBUG_FLAGS = -g -O0
MODE_FLAGS = $(RELEASE_FLAGS)

PROJ_FLAGS = $(MODE_FLAGS) -c -fpic -pedantic -Wall -Isrc -Ithird_party/include
TEST_FLAGS = $(MODE_FLAGS) -pedantic -Wall -Isrc -Ithird_party/include -I/home/aj/.local/include

all: tests

$(BUILDDIR):
	mkdir -p $(shell find src -type d | sed -e "s/^/$(BUILDDIR)\//")

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
