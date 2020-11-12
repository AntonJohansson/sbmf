SRCS = src/main.c
LIBS = -fopenmp -lgsl -lgfortran -L/home/aj/.local/lib -l:sbmf.a -l:cimgui.so -lplot -lglfw -lm

explore_gp: ${SRCS}
	gcc -o explore_gp ${SRCS} -g -I/home/aj/.local/include ${LIBS}
