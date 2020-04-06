CC = gcc
PROJECT = sbmf
BUILD_DIR = build

PROJ_SRCS = \
	src/sbmf/groundstate_solver/groundstate_solver.c  \
	src/sbmf/groundstate_solver/item.c \
	src/sbmf/common/profile.c
TEST_SRCS = $(shell find src/test -name "*.c")

PROJ_LIBS = -lfftw3 -llapacke -lm
TEST_LIBS = -l$(PROJECT) -llapacke -lm

default: $(BUILD_DIR)/$(PROJECT) $(BUILD_DIR)/test

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# -Ithird_party/include -Lthird_party/lib  required on school comp
#  How do I check for that?
$(BUILD_DIR)/$(PROJECT): $(BUILD_DIR) $(PROJ_SRCS)
	$(CC) -Wall -Werror -o $(BUILD_DIR)/lib$(PROJECT).so -shared -fpic $(PROJ_SRCS) -Isrc -g $(PROJ_LIBS)

$(BUILD_DIR)/test: $(TEST_SRCS) $(BUILD_DIR)/$(PROJECT)
	$(CC) $(TEST_SRCS) -o $(BUILD_DIR)/test -Wl,-rpath=$(BUILD_DIR) -Isrc -L$(BUILD_DIR) $(TEST_LIBS)
