TARGET := substring_search
CC := g++

CFLAGS = -Wall -Wextra
LDFLAGS =

SRCPATH = src

# HEADERS = $(shell find $(SRCPATH) -name '*.h' | sort -k 1nr | cut -f2-)
SOURCES = $(shell find $(SRCPATH) -name '*.cpp' | sort -k 1nr | cut -f2-)
OBJECTS = $(SOURCES:%.cpp=%.o)

.PHONY: all

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(TARGET)
	rm -f $(OBJECTS)
