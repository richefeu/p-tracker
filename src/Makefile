TOOFUSPATH = ~/toofus

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)

	# You can install clang-omp by using homebrew: brew install clang-omp
	CXX = g++-13
	
	# -Wshadow
	CXXFLAGS = -fopenmp -O3 -W -Wall -Wextra -pedantic -std=c++11 
	CXXFLAGS += -I $(TOOFUSPATH)
	CXXFLAGS += `pkg-config --cflags libraw`
	CXXFLAGS += `pkg-config --cflags libtiff-4`
	
	OPENMP_LIB = -l gomp

else

	CXX = g++
	CC_OPT = -fopenmp -O3 -pedantic -W -Wall -std=c++11
	CXXFLAGS += -I $(TOOFUSPATH)
	CXXFLAGS += `pkg-config --cflags libraw`
	CXXFLAGS += `pkg-config --cflags libtiff-4`
	
endif

# All the files that are part of the source code
SOURCES = p-tracker.cpp distortion_correction.cpp \
          grain_type_2D.cpp NCC_optimizer_functor.cpp interpol.cpp 
OBJECTS = $(SOURCES:.cpp=.o)

%.o:%.cpp
	@echo "\033[0;32m-> COMPILING" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c $< -o $@

all: ptracker

ptracker: run.cpp $(OBJECTS)
	@echo "\033[0;32m-> BUILDING" $@ "\033[0m"
	$(CXX) $(CXXFLAGS) -c run.cpp -o run.o
	$(CXX) -o ptracker run.o $(OBJECTS) `pkg-config --libs libraw` `pkg-config --libs libtiff-4` $(OPENMP_LIB)

clean:
	rm -f *.o ptracker
