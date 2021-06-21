PROJECT = gp

CC = /usr/local/cuda/bin/nvcc
SRCS = src/main.c

$(PROJECT): $(SRCS)
	$(CC) -o $(PROJECT) $(SRCS) -O3 -L../../build -I../../include -l:sbmf.a -lm -lpthread -lgfortran -Xcompiler -fopenmp
