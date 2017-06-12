CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp -DNDEBUG
LDFLAGS =
BIN = ./bin
SRC = ./src

.PHONY: steiner

all: steiner

rfw_timer:
	$(CXX) $(CXXFLAGS) -o $(BIN)/rfw_timer.o -c $(SRC)/rfw_timer.cpp

rfw_random:
	$(CXX) $(CXXFLAGS) -o $(BIN)/rfw_random.o -c $(SRC)/rfw_random.cpp

steiner: rfw_timer rfw_random
	$(CXX) $(CXXFLAGS) -o $(BIN)/steiner.o -c $(SRC)/steiner.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/steiner $(BIN)/steiner.o $(BIN)/rfw_random.o $(BIN)/rfw_timer.o
	rm $(BIN)/rfw_random.o $(BIN)/rfw_timer.o $(BIN)/steiner.o

clean:
	rm -f $(BIN)/*
