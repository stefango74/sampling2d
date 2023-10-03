CC=g++
CFLAGS=-c -Wall -O3 -I/usr/include/eigen3
LDFLAGS=-lann -lgmp
SOURCES=sample2d.cpp PrecTimer.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=sample2d

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
 
 clean:
	rm -rf $(OBJECTS) $(EXECUTABLE)
