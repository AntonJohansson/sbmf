PROJECT = bestmf_test

SRCS = src/main.c

DEBUG_FLAGS = -g \
              -fsanitize=address \
              -fsanitize=leak \
              #-fsanitize=thread \#
              -fsanitize=undefined \
              -fsanitize=bool \
              -fsanitize=enum \
              -fsanitize=float-cast-overflow \


$(PROJECT): $(SRCS)
	$(CC) -o $(PROJECT) $(SRCS) $(DEBUG_FLAGS) -O0 -L/home/aj/.local/lib/ -I/home/aj/.local/include -lgsl -l:sbmf.a -lplot -lm -lpthread -fopenmp -lgfortran -lgsl
