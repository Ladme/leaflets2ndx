leaflets2ndx: main.c
	gcc main.c -I$(groan) -L$(groan) -D_POSIX_C_SOURCE=200809L -o leaflets2ndx -lgroan -lm -std=c99 -pedantic -Wall -Wextra -O3 -march=native

install: leaflets2ndx
	cp leaflets2ndx ${HOME}/.local/bin
