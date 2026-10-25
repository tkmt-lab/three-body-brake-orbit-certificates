CXX := clang++
CXXFLAGS := -std=c++17 -O2 -DFAST_INTERVAL -Wall -Wextra -I/opt/homebrew/include
LDLIBS := -L/opt/homebrew/lib -lmpfr -lgmp

.PHONY: all
all: bin/certify bin/probe bin/check bin/check_multiply

bin/%: src/%.cpp src/interval.hpp src/dynamics.hpp
	mkdir -p bin
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@
