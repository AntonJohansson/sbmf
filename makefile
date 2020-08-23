CC = gcc
PROJECT = sbmf

BUILDDIR = build

PROJ_SRCS = \
	src/sbmf/sbmf.c \
	src/sbmf/methods/quadgk.c \
	src/sbmf/methods/item.c \
	src/sbmf/methods/hob.c \
	src/sbmf/math/find_eigenpairs.c \
	src/sbmf/math/matrix.c \
	src/sbmf/debug/profile.c \
	src/sbmf/debug/log.c

PROJ_OBJS = $(patsubst %.c, $(BUILDDIR)/%.o, $(PROJ_SRCS))

TEST_SRCS = src/test/main.c

PROJ_LIBS = \
	-lm
TEST_LIBS = \
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

PROJ_FLAGS = -g -c -fpic -O3 -Wall -Isrc -Ithird_party/include
TEST_FLAGS = -g -O0 -Isrc -Ithird_party/include -I/home/aj/.local/include

all: tests

$(BUILDDIR):
	mkdir -p $(shell find src -type d | sed -e "s/^/$(BUILDDIR)\//")

$(BUILDDIR)/%.o: %.c
	$(CC) -o $@ $^ -static $(PROJ_FLAGS) $(PROJ_LIBS)

$(BUILDDIR)/$(PROJECT).a: $(BUILDDIR) $(PROJ_OBJS)
	ar rcs $@ $(PROJ_OBJS)

$(BUILDDIR)/tests: $(TEST_SRCS) $(BUILDDIR)/$(PROJECT).a
	$(CC) $(TEST_SRCS) -o $(BUILDDIR)/tests $(TEST_FLAGS) $(TEST_LIBS)

.PHONY: lib
lib: $(BUILDDIR)/$(PROJECT).a

.PHONY: tests
tests: $(BUILDDIR)/tests

.PHONY: clean
clean:
