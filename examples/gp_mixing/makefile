PROJECT = gp

SRCS = src/main.c

DEBUG_FLAGS = -g \
              -fsanitize=address \
              -fsanitize=leak \
              -fsanitize=undefined \
              -fsanitize=bool \
              -fsanitize=enum \
              -fsanitize=float-cast-overflow \
              -fsanitize=signed-integer-overflow

$(PROJECT): $(SRCS)
	$(CC) -o $(PROJECT) $(SRCS) -L/home/aj/.local/lib/ -I/home/aj/.local/include -lgsl -l:sbmf.a -lm -lpthread -fopenmp -lgfortran $(DEBUG_FLAGS)
