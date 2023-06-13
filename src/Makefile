CC      =   g++ --std=c++17
CFLAGS  = 	-Werror -Wall -Wextra -Wpedantic -Wcast-align -Wcast-qual \
		 	-Wconversion -Wctor-dtor-privacy -Wenum-compare -Wfloat-equal \
		 	-Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wredundant-decls \
		 	-Wsign-conversion -Wsign-promo 

build: clean
	$(CC) $(CFLAGS) *.cc -o main
	./main

clean: 
	@rm -rf main