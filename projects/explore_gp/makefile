SRCS = src/main.c
LIBS = -fopenmp -L/home/aj/.local/lib -l:sbmf.a -lplot -lglfw -lm -lcimgui

explore_gp: ${SRCS}
	gcc -o explore_gp ${SRCS} -g -I/home/aj/.local/include ${LIBS}
