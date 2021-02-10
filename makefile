PROJECT = sbmf
BUILDDIR = build
INSTALLDIR = ~/.local
CUDADIR = /usr/local/cuda

CC = gcc
NVCC = $(CUDADIR)/bin/nvcc

PROJ_SRCS = src/sbmf/sbmf.c
PROJ_OBJS = $(patsubst %.c, $(BUILDDIR)/%.o, $(PROJ_SRCS))

CUDA_SRCS = src/sbmf/cuda/sbmf_cuda.cu
CUDA_OBJS = $(patsubst %.cu, $(BUILDDIR)/%.o, $(CUDA_SRCS))

PROJ_LIBS = 							\
	-lm 								\
	-fopenmp 							\
	-l:third_party/lib/libarpack.a 		\
	-l:third_party/lib/libopenblas.a 	\
	-lgfortran

RELEASE_FLAGS = -O3
DEBUG_FLAGS = -g \
			  -fsanitize=address \
			  -fsanitize=leak \
			  -fsanitize=undefined \
			  -fsanitize=bool \
			  -fsanitize=enum \
			  -fsanitize=float-cast-overflow \
			  -fsanitize=signed-integer-overflow

PROJ_FLAGS = -c -fPIC -pedantic -Wall -Wextra -Ithird_party/include -Iinclude
CUDA_CC_FLAGS = -c -fPIC -Wall -Wextra -Iinclude
CUDA_FLAGS = -c -I $(CUDADIR)/include

.PHONY: release
release: PROJ_FLAGS += $(RELEASE_FLAGS)
release: PROJ_FLAGS += $(PROJ_LIBS)
release: CUDA_FLAGS := $(CUDA_FLAGS) -Xcompiler "$(CUDA_CC_FLAGS) $(RELEASE_FLAGS) $(PROJ_LIBS)"
release: $(BUILDDIR)/$(PROJECT).a

.PHONY: debug
debug: PROJ_FLAGS += $(DEBUG_FLAGS)
debug: PROJ_FLAGS += $(PROJ_LIBS)
debug: CUDA_FLAGS := $(CUDA_FLAGS) -Xcompiler "$(CUDA_CC_FLAGS) $(DEBUG_FLAGS) $(PROJ_LIBS)"
debug: $(BUILDDIR)/$(PROJECT).a

.PHONY: install
install:
	mkdir -p $(INSTALLDIR)/lib
	mkdir -p $(INSTALLDIR)/include/
	cp $(BUILDDIR)/$(PROJECT).a $(INSTALLDIR)/lib/
	cp -r include/sbmf $(INSTALLDIR)/include/

$(BUILDDIR):
	mkdir -p $(shell find src -type d | sed -e "s/^/$(BUILDDIR)\//")

third_party/lib:
	cd third_party && sh build.sh

$(BUILDDIR)/%.o: %.c
	$(CC) -o $@ $^ $(PROJ_FLAGS)

$(BUILDDIR)/%.o: %.cu
	$(NVCC) -o $@ $^ $(CUDA_FLAGS)

$(BUILDDIR)/$(PROJECT).a: $(BUILDDIR) $(PROJ_OBJS) $(CUDA_OBJS)
	mkdir -p $(BUILDDIR)/tmp
	cd $(BUILDDIR)/tmp && ar x ../../third_party/lib/libarpack.a
	cd $(BUILDDIR)/tmp && ar x ../../third_party/lib/libopenblas.a
	ar rcs $@ $(PROJ_OBJS) $(CUDA_OBJS) $(BUILDDIR)/tmp/*.o
	rm -r $(BUILDDIR)/tmp


.PHONY: clean
clean:
	rm -rf $(BUILDDIR)
