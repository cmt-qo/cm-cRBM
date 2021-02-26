.PHONY: all clean

.DEFAULT: all


CXX = c++

CXXFLAGS = -ansi -pedantic -std=c++11 #-fopenmp

TARGETS = run_MC

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

$(TARGETS) : %: main.cc
	$(CXX) $(CXXFLAGS) -o $@ $<
