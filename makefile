CC = gcc
PROJECT = sbmf
BUILD_DIR = build

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
PLOT_SRCS = $(shell find src/plot -name "*.c")

PROJ_LIBS = -lfftw3 -llapacke -larpack -lm
PLOT_LIBS = -lglfw -lpthread -lm -ldl
TEST_LIBS = -l$(PROJECT) -lplot -lm -larpack -lcblas

#PROJ_FLAGS = -fsanitize=address -fsanitize=leak -g -fpic -Wall -Werror -Isrc
#TEST_FLAGS = -fsanitize=address -fsanitize=leak -g -Wall -Werror -Isrc -Isrc/test/plotting/third_party
PROJ_FLAGS = -g -fpic -Wall -Werror -Isrc
PLOT_FLAGS = -g -fpic -Isrc/plot/third_party/ -O3
TEST_FLAGS = -g -Isrc -fsanitize=address -fsanitize=leak

default: $(BUILD_DIR)/lib$(PROJECT) $(BUILD_DIR)/test $(BUILD_DIR)/libplot

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# -Ithird_party/include -Lthird_party/lib  required on school comp
#  How do I check for that?
$(BUILD_DIR)/lib$(PROJECT): $(BUILD_DIR) $(PROJ_SRCS)
	$(CC) -o $(BUILD_DIR)/lib$(PROJECT).so -shared $(PROJ_SRCS) $(PROJ_FLAGS) $(PROJ_LIBS)

$(BUILD_DIR)/libplot: $(BUILD_DIR) $(PLOT_SRCS)
	$(CC) -o $(BUILD_DIR)/libplot.so -shared $(PLOT_SRCS) $(PLOT_FLAGS) $(PLOT_LIBS)

$(BUILD_DIR)/test: $(TEST_SRCS) $(BUILD_DIR)/lib$(PROJECT) $(BUILD_DIR)/libplot
	$(CC) $(TEST_SRCS) -o $(BUILD_DIR)/test -Wl,-rpath=$(BUILD_DIR) -L$(BUILD_DIR) $(TEST_FLAGS) $(TEST_LIBS)

