CC=gcc
CFLAGS = -O0 -funroll-loops -g
LDFLAGS = -g
SOURCES = psvpan.c cmdline.c common.c alignment.c reference.c interval_tree.c free.c sv.c
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
