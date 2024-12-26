CC = g++
TARGET = main

RELEASE := 1
ifeq ($(RELEASE),1)
CXXFLAGS += -O3 
else
CXXFLAGS += -g
endif

#CXXFLAGS += -DENABLE_WRITE_FILE

# Add the Eigen include path
CXXFLAGS += -I/usr/include/eigen3

LIBS_PATH = -L/usr/local/lib/OpenMesh
LIBS_PATH += -L/usr/local/lib

LIBS = -lGLEW -lGL -lGLU -lglfw -lX11  -lXrandr -lpthread -lXi
LIBS += -lOpenMeshCore -lOpenMeshTools

SRC := $(shell find . -name "*.cpp")

all:
	@make $(TARGET)


%.o : %.cpp
	@echo "Compiling $< ..."
	@$(CC) $(LIBS_PATH) $(LIBS) $(CXXFLAGS) -o $@ -c $<

$(TARGET): $(SRC:.cpp=.o)
	@echo "Linking $@..."
	@$(CC) -o $@ $(SRC:.cpp=.o) $(LIBS_PATH) $(LIBS) 

.PHONY: clean
clean:
	\rm -f *.o
	\rm -f $(TARGET)
