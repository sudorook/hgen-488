TARGET := align_blosum62
CC ?= gcc

# CFLAGS = -Wall -Wextra -g $(shell pkg-config --cflags gsl) -I/usr/include/gsl
# LDFLAGS = $(shell pkg-config --libs gsl)
CFLAGS = -Wall -Wextra -g
LDFLAGS = -lm

SRC_PATH = src

HEADERS = $(shell find $(SRC_PATH) -name '*.h' | sort -k 1nr | cut -f2-)
SOURCES = $(shell find $(SRC_PATH) -name '*.c' | sort -k 1nr | cut -f2-)
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
