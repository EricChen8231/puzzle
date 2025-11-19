CXX = g++
CXXFLAGS = -std=c++17 -Wall -g -O0

SOURCES = main.cpp sudoku.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = puzzle_solver

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(OBJECTS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
