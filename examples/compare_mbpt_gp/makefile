PROJECT = gp

SRCS = src/main.c
CC = /usr/local/cuda/bin/nvcc

DEBUG_FLAGS = -g \
              -fsanitize=address \
              -fsanitize=leak \
              -fsanitize=undefined \
              -fsanitize=bool \
              -fsanitize=enum \
              -fsanitize=float-cast-overflow \
              -fsanitize=signed-integer-overflow

$(PROJECT): $(SRCS)
	$(CC) -o $(PROJECT) $(SRCS) -O3 -L../../build -I../../include -l:sbmf.a -lm -lpthread -lgfortran -Xcompiler -fopenmp
