#
# Copyright 2017 Lime Microsystems Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
MAKEFLAGS += --silent
CCFLAGS := -Ofast -march=native
CCFLAGS += -ftree-vectorize
CCFLAGS += -funsafe-math-optimizations
#CCFLAGS += -fopt-info-vec-optimized
#CCFLAGS += -Wall
CC := g++
LDFLAGS :=
LD := g++
EXEC := LimeScan
EXEC2 := $(patsubst %, bin/%, $(EXEC))
LIBS := -lm -lfftw3 -lLimeSuite
CCINCLDS := -Iinclude
SRCDIR := src/
OBJDIR := obj/
SRC := LimeScan.cpp 
OBJ := $(patsubst %.c, $(OBJDIR)%.o, $(filter %.c, $(SRC)))
OBJ += $(patsubst %.cpp, $(OBJDIR)%.o, $(filter %.cpp, $(SRC)))

all: $(EXEC)
	@echo Done!

$(EXEC): $(OBJ)
	$(LD) $(OBJ) $(LIBS) $(LDFLAGS) -o $@ 

$(OBJDIR)%.o: $(SRCDIR)%.c
	$(CC) -c $(SRCDIR)$*.c $(LIBS) $(CCINCLDS) $(CCFLAGS) -o $(OBJDIR)$*.o 

$(OBJDIR)%.o: $(SRCDIR)%.cpp
	$(CC) -c $(SRCDIR)$*.cpp $(LIBS) $(CCINCLDS) $(CCFLAGS) -o $(OBJDIR)$*.o

clean:
	-rm -f $(OBJ)
	-rm -f $(EXEC)

diag:
	@echo OBJ: $(OBJ)
	@echo SRC: $(SRC)
	@echo LIBS: $(LIBS)
	@echo CCINCLDS: $(CCINCLDS)
#	@echo EXEC: $(EXEC2)
	


