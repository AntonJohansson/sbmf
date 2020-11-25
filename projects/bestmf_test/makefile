PROJECT = bestmf_test

SRCS = src/main.c

$(PROJECT): $(SRCS)
	$(CC) -o $(PROJECT) $(SRCS) -g -O0 -L/home/aj/.local/lib/ -I/home/aj/.local/include -lgsl -l:sbmf.a -lplot -lm -lpthread -fopenmp -lgfortran
