CXX ?= g++
CXXFLAGS ?= -std=c++20 -O3 -Wall
LDFLAGS ?= -L$(PACS_ROOT)/lib -lpacs
CPPFLAGS=-I$(PACS_ROOT)/include
SRCS = main.cpp

# Header files
HEADERS = $(wildcard *.hpp)

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
EXEC = main

# Default target
all: $(EXEC)

# Rule to link the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Rule to compile each source file into an object file
%.o: %.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(EXEC) $(OBJS)

# Phony target to avoid conflicts with files named clean
.PHONY: clean
