# Python script and output
DATA_SCRIPT = generate_data.py
DATA_FILE = test_data.h

# Default sizes
N ?= 65536
M ?= 32768

# Source files
SRC = main.c implementations.c
OBJ = $(SRC:.c=.o)
EXEC = main

# Compiler
CC = gcc
CFLAGS = -g -Wall
PYTHON = python3

all: $(EXEC)

# Generate test_data.h
$(DATA_FILE): $(DATA_SCRIPT)
	$(PYTHON) $(DATA_SCRIPT) $(N) $(M)

# Compile .c files into .o
%.o: %.c $(DATA_FILE)
	$(CC) $(CFLAGS) -c $< -o $@

# Link executable
$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

# Clean
clean:
	rm -f $(OBJ) $(EXEC) $(DATA_FILE)
