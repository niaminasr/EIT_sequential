EXECUTABLE = solveHeat

SPARSE_DIR = #/home/niami/Code_C++/project/EIT_pb_directe/SuiteSparse
EIGEN_DIR = ./Eigen/Eigen/
BUILD_DIR = ./build

CXX = g++
CXXFLAGS = -std=c++14 -I$(SPARSE_DIR)/include -I$(EIGEN_DIR) -O3 # -DNDEBUG -Wno-maybe-uninitialized -Wno-unused-variable -Wno-sign-compare
LDFLAGS = #-L$(SPARSE_DIR)/lib -Wl,-R$(SPARSE_DIR)/lib '-Wl,-R'
LIBRARIES = #-lgfortran -lumfpack -lcholmod -lamd -lcolamd -lcamd -lccolamd -lmetis $(LAPACK) $(BLAS)

SRCS = parameters.cpp stringTools.cpp  grid.cpp  assembler.cpp assembler_forme.cpp heatProblem.cpp inverse.cpp inverseforme.cpp writer.cpp testHeat.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
DEPS =  $(SRCS:%.cpp=%.d)


.PHONY: all clean

# all:$(BUILD_DIR)/$(EXECUTABLE)

#$(BUILD_DIR)/$(EXECUTABLE) : $(OBJS)
#	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBRARIES)
#	cp params.in $(BUILD_DIR)/params.in

all:$(EXECUTABLE)

$(EXECUTABLE) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBRARIES) 

# Include dependencies
-include $(DEPS)

$(BUILD_DIR)/%.o : %.cpp
	mkdir -p '$(@D)'
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

clean:
	rm -r $(BUILD_DIR)
	rm -rf $(EXECUTABLE)
