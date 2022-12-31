CC=clang
CFLAGS = -O3 -funroll-loops -g
LDFLAGS = -lz -lm -lpthread -llzma -lbz2
SOURCES = psvpan.c cmdline.c common.c alignment.c reference.c
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

