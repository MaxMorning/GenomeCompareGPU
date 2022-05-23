NVCC := nvcc
GCC := g++
FLAGS := -arch=sm_75 -O3 
TARGET := compare.exe preprocess.exe

all: $(TARGET)

compare.exe : compare.cu
	$(NVCC) $< -o $@ $(FLAGS)

preprocess.exe : preprocess.cpp
	$(GCC) $< -o $@

.PHONY: clean
clean:
	-del compare.exe del compare.exp compare.lib
	-del preprocess.exe preprocess.obj