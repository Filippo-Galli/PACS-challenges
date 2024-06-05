# Define the compiler
CXX = mpicxx

# Define compiler flags
CXXFLAGS = -DMPI -Wall -Wextra -std=c++20 -fopenmp -lmuparser -lomp
DEBUG_CXXFLAGS = -DDEBUG

# Define the number of threads to use for the build
MAKEFLAGS += -j2

# Define the paths to the source, include and objects directories
SRC_DIR = src
INCLUDE_DIR = include
OBJ_DIR = obj

# Define the target executable name
TARGET = main

# Define the source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Define the object files with the object directory path
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# Define the include path
INCLUDES = -I$(INCLUDE_DIR)

# Check if obj directory exists or create it
$(shell mkdir -p $(OBJ_DIR))

# Check if vtk_files directory exists or create it
$(shell mkdir -p vtk_files)

# Define the rule to build the target
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

# Define the rule to compile the source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Define the rule to clean up the build artifacts
clean:
	rm -f $(OBJ_DIR)/*.o $(TARGET)

# Define the rule to build the target in debug mode
debug: CXXFLAGS += -Og $(DEBUG_CXXFLAGS)
debug: $(TARGET)

# Define the rule to compile with mpi the target
mpi: CXXFLAGS += -O3 
mpi: $(TARGET)

# Define the rule to compile and run the tests
test_compile: CXXFLAGS += -O3 -DTEST=1
test_compile: $(TARGET)
