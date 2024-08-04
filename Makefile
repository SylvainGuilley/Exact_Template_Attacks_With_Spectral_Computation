ALL: main.x main_Template_coalescence.x

all: $(ALL)

main.x: main.c Matrice.h
	gcc -g -O0 -Wall main.c -o $@ -lm

main_Template_coalescence.x: main_Template_coalescence.c Matrice.h
	gcc -g -O0 -Wall main_Template_coalescence.c -o $@ -lm

clean:
	rm -f *~
	rm -f $(ALL)
