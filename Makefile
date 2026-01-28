# Default compiler settings (Linux/Windows WSL)
CXX      := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -fopenmp
LDFLAGS  := -fopenmp

# Detect OS
# If on macOS (Darwin), switch to clang++ and remove OpenMP flags to avoid errors
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    CXX      := clang++
    CXXFLAGS := -std=c++17 -O3 -Wall -Wextra
    LDFLAGS  :=
endif

TARGET  := phase_field_model
SOURCES := main.cpp PhaseFieldModel.cpp
OBJECTS := $(SOURCES:.cpp=.o)
HEADERS := PhaseFieldModel.h

.PHONY: all clean run delete

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

# Compile source files to object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run the program
run: $(TARGET)
	./$(TARGET)

# Clean build files (objects and executable)
clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *.dat

# Clean build files AND remove the data directory
delete: clean
	rm -rf data/
	rm *.txt