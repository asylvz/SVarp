CC=g++
CFLAGS = -O0 -funroll-loops -g -Wall
LDFLAGS = -g
SOURCES = psvpan.cpp common.cpp cmdline.cpp reference.cpp alignment.cpp sv.cpp interval_tree.cpp
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = psvpan

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) *.o *~

#clang
