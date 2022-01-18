#CC      = g++
CC      = mpic++
CPPFLAGS= -Wno-deprecated -O3 -c -m64 -march=native -std=c++17 -DGCC -DRelease
#LDFLAGS = -static -O3 -m64
SOURCES = main-eval.cpp Util.cpp Graph.cpp UnitHeap.cpp UnitHeap.h
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE=SubGorder-eval

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o

