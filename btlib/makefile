COMP_ATTRS=-O2 -Wall

O_FILES=balltree.o

compile: $(O_FILES)
	g++ main.cpp -o test $(O_FILES) $(COMP_ATTRS)

balltree.o: balltree.h balltree.cpp
	g++ -c balltree.cpp $(COMP_ATTRS)

clean:
	rm *.o
