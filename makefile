CC = gcc
PROJECT = sbmf

PROJ_SRCS = \
	src/sbmf/groundstate_solver/groundstate_solver.c  \
	src/sbmf/groundstate_solver/item.c \
	src/sbmf/common/profile.c \
	src/sbmf/common/eigenproblem.c \
	src/sbmf/common/matrix.c \
	src/sbmf/common/log.c \
	src/sbmf/common/common.c \
	src/sbmf/quadgk.c
TEST_SRCS = $(shell find src/test -name "*.c")

PROJ_LIBS = -lfftw3 -llapacke -larpack -lm
TEST_LIBS = -l$(PROJECT) -lm -larpack -lcblas -L. -L/home/aj/.local/lib -lplot -I/home/aj/.local/include
#PROJ_FLAGS = -fsanitize=address -fsanitize=leak -g -fpic -Wall -Werror -Isrc
#TEST_FLAGS = -fsanitize=address -fsanitize=leak -g -Wall -Werror -Isrc -Isrc/test/plotting/third_party
PROJ_FLAGS = -g -fpic -Wall -Werror -Isrc
TEST_FLAGS = -g -Isrc

default: lib$(PROJECT).so test

# -Ithird_party/include -Lthird_party/lib  required on school comp
#  How do I check for that?
lib$(PROJECT).so: $(PROJ_SRCS)
	$(CC) -o lib$(PROJECT).so -shared $(PROJ_SRCS) $(PROJ_FLAGS) $(PROJ_LIBS)

test: $(TEST_SRCS) lib$(PROJECT).so
	$(CC) $(TEST_SRCS) -o test $(TEST_FLAGS) $(TEST_LIBS)
