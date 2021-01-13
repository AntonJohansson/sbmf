PROJECT = pt

SRCS = src/main.c

DEBUG_FLAGS = -g \
              -fsanitize=address \
              -fsanitize=leak \
              -fsanitize=undefined \
              -fsanitize=bool \
              -fsanitize=enum \
              -fsanitize=float-cast-overflow \
              -fsanitize=signed-integer-overflow

RELEASE_FLAGS = -O3

$(PROJECT): $(SRCS)
	$(CC) -o $(PROJECT) $(SRCS) $(RELEASE_FLAGS) -L/home/aj/.local/lib/ -I/home/aj/.local/include -lgsl -l:sbmf.a -lplot -lm -lpthread -fopenmp -lgfortran
