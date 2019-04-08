TARGET := align_blosum62
CC ?= gcc

CFLAGS = -Wall -Wextra -g  -std=c99
LDFLAGS = -lm

SRCPATH = src

HEADERS = $(shell find $(SRCPATH) -name '*.h' | sort -k 1nr | cut -f2-)
SOURCES = $(shell find $(SRCPATH) -name '*.c' | sort -k 1nr | cut -f2-)
OBJECTS = $(SOURCES:%.c=%.o)

.PHONY: all

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJECTS): $(SOURCES) $(HEADERS)
	$(CC) -c $(CFLAGS) $(LDFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f $(TARGET)
	rm -f $(OBJECTS)
