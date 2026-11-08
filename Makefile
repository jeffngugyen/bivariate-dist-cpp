CXX ?= c++
PYTHON ?= python3
CXXFLAGS ?= -O3 -Wall -Wextra -Wpedantic
CPPFLAGS ?=
LDFLAGS ?=
BUILD_DIR ?= build

# Apple Clang builds serially by default; Linux GCC enables OpenMP.
# Override OPENMP=0 for a compiler without OpenMP, or OPENMP_FLAGS as needed.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
OPENMP ?= 0
MODULE_FLAGS = -shared -undefined dynamic_lookup
else
OPENMP ?= 1
MODULE_FLAGS = -shared
endif
OPENMP_FLAGS ?= -fopenmp
ifeq ($(OPENMP),1)
PARALLEL_FLAGS = $(OPENMP_FLAGS)
endif
COMPILE_FLAGS = $(CPPFLAGS) $(CXXFLAGS) -std=c++20 $(PARALLEL_FLAGS)

# An initialized, pinned submodule avoids requiring a pip installation.
ifneq ($(wildcard extern/pybind11/include/pybind11/pybind11.h),)
PYBIND_INCLUDES = -Iextern/pybind11/include -I$(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_path("include"))')
else
PYBIND_INCLUDES = $(shell $(PYTHON) -m pybind11 --includes 2>/dev/null)
endif
EXT_SUFFIX = $(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')

.PHONY: all dist example test python dist_py python-test clean
all: $(BUILD_DIR)/dist

dist: $(BUILD_DIR)/dist

$(BUILD_DIR):
	mkdir -p $@

# Recompile on flag changes as well as source changes. This small project has
# only two translation units; the phony prerequisite avoids stale serial builds.
.PHONY: FORCE
FORCE:

$(BUILD_DIR)/dist: dist.cpp dist.h main.cpp Makefile FORCE | $(BUILD_DIR)
	$(CXX) $(COMPILE_FLAGS) dist.cpp main.cpp $(LDFLAGS) -o $@

example: all
	./$(BUILD_DIR)/dist examples/four_qubit/H_X.txt examples/four_qubit/H_Z.txt 1 2

test: all
	OMP_NUM_THREADS=1 DIST_BINARY=$(abspath $(BUILD_DIR)/dist) $(PYTHON) -m unittest discover -s tests -v
	OMP_NUM_THREADS=4 DIST_BINARY=$(abspath $(BUILD_DIR)/dist) $(PYTHON) -m unittest discover -s tests -v

python dist_py: all
	@test -n "$(PYBIND_INCLUDES)" || (echo 'Run git submodule update --init extern/pybind11 or install pybind11 for $(PYTHON).' >&2; exit 1)
	$(CXX) $(COMPILE_FLAGS) -fPIC $(MODULE_FLAGS) $(PYBIND_INCLUDES) dist.cpp dist_py.cpp $(LDFLAGS) -o dist_py$(EXT_SUFFIX)

python-test: python
	OMP_NUM_THREADS=1 DIST_BACKEND=python DIST_BINARY=$(abspath $(BUILD_DIR)/dist) $(PYTHON) -m unittest discover -s tests -v
	OMP_NUM_THREADS=4 DIST_BACKEND=python DIST_BINARY=$(abspath $(BUILD_DIR)/dist) $(PYTHON) -m unittest discover -s tests -v

clean:
	rm -f $(BUILD_DIR)/dist dist_py*.so dist_py*.pyd
