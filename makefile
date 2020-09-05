COMP_ATTRS = -O2 -Wall
#COMP_ATTRS = -O2 -Wall -I /usr/include/eigen3

O_FILES = main.o core.o matrix.o model.o visualization.o lik.o player.o playerexperim.o pergen.o periodic.o dynrec.o ftsolver.o odestate.o cpc.o effdata.o geom.o ghost.o

#LINK_ATTRS = -L /usr/local/lib -lode -ldrawstuff -lm -lstdc++ -lGL -lGLU -lglut -lX11 -pthread -Lballtreelib/ -lballtree
#LINK_ATTRS = -L /usr/local/lib -lode -ldrawstuff -lm -lstdc++ -lGL -lGLU -lglut -lX11 -pthread -Lballtreelib1/ -lballtree
LINK_ATTRS = -L /usr/local/lib -lode -ldrawstuff -lm -lstdc++ -lGL -lGLU -lglut -lX11 -pthread -Lbtlib/ -lballtree
# -I /usr/include/eigen3
#PRECISION = dSINGLE

EIGEN = -I /usr/include/eigen3

all: coordconv

coordconv: $(O_FILES)
	g++ $(O_FILES) -o coordconv $(LINK_ATTRS) $(BALLTREE)

main.o: main.cpp
	g++ -c main.cpp $(COMP_ATTRS) $(EIGEN)

core.o: core.h core.cpp
	g++ -c core.cpp $(COMP_ATTRS)

matrix.o: matrix.h matrix.cpp
	g++ -c matrix.cpp $(COMP_ATTRS)

model.o: model.h model.cpp
	g++ -c model.cpp $(COMP_ATTRS)

visualization.o: visualization.h visualization.cpp
	g++ -c visualization.cpp $(COMP_ATTRS)

lik.o: lik.h lik.cpp
	g++ -c lik.cpp $(COMP_ATTRS)

player.o: player.h player.cpp
	g++ -c player.cpp $(COMP_ATTRS) $(EIGEN)

playerexperim.o: player.h playerexperim.cpp
	g++ -c playerexperim.cpp $(COMP_ATTRS) $(EIGEN)

pergen.o: pergen.h pergen.cpp
	g++ -c pergen.cpp $(COMP_ATTRS)

periodic.o: periodic.h periodic.cpp
	g++ -c periodic.cpp $(COMP_ATTRS) $(EIGEN)

dynrec.o: dynrec.h dynrec.cpp
	g++ -c dynrec.cpp $(COMP_ATTRS) $(EIGEN)

ftsolver.o: ftsolver.h ftsolver.cpp
	g++ -c ftsolver.cpp $(COMP_ATTRS) $(EIGEN)

odestate.o: odestate.h odestate.cpp
	g++ -c odestate.cpp $(COMP_ATTRS)

cpc.o: cpc.h cpc.cpp
	g++ -c cpc.cpp $(COMP_ATTRS) $(EIGEN)

effdata.o: effdata.h effdata.cpp
	g++ -c effdata.cpp $(COMP_ATTRS)

geom.o: geom.h geom.cpp
	g++ -c geom.cpp $(COMP_ATTRS)

ghost.o: ghost.h ghost.cpp
	g++ -c ghost.cpp $(COMP_ATTRS)


clean:
	rm *.o