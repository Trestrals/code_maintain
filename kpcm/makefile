CC=g++
CPPFLAGS=
LDFLAGS=
DEPS = graph.h utility.h experiment.h
OBJ =  graph.o main.o utility.o experiment.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++11 $(LDFLAGS) -c -O3 -o $@ $< $(CPPFLAGS)

kpcore: $(OBJ)
	$(CC) -std=c++11 -pthread $(LDFLAGS) -O3 -o $@ $^ $(CPPFLAGS)