# This is the makefile for testing the mesa colors library

#MESA_DIR = $(MESA_DIR)

#################################################################

# STEP 1: get the standard compiler setup

include $(MESA_DIR)/utils/makefile_header

ifeq ($(COLORS_OFFLOAD),YES)
CC += -Doffload
FC += -Doffload
endif

#################################################################

# STEP 2: info about the libraries

LOCAL_INCLUDES = -I$(LOCAL_LIB_DIR)
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lcolor
LOAD_OTHER = #-L$(MESA_LIB_DIR) -linterp_2d -linterp_1d -lnum $(LOAD_MATRIX)

#################################################################

# STEP 3: build

TEST = tester
TEST_OBJS = test_color.o

all : $(TEST)

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $(TEST_DIR)/$(TEST) $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

#################################################################

%.o: $(TEST_SRC_DIR)/%.f90
	$(TEST_COMPILE) $<

clean:
	-@rm -f *.o *.mod $(TEST_DIR)/$(TEST)
