# This is the makefile for testing the mesa colors library

MESA_DIR = /data/mesa7324

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
LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lcolors
LOAD_OTHER = -L$(MESA_LIB_DIR) -linterp_2d -linterp_1d -lnum $(LOAD_MATRIX)

#################################################################

# STEP 3: build

TEST = tester
TEST_OBJS = \
   test_colors.o

all : $(TEST)

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $(TEST_DIR)/$(TEST) $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

#################################################################

%.o: $(TEST_SRC_DIR)/%.f90
	$(TEST_COMPILE) $<

%.mod: %.o
	@true

clean:
	-@rm -f *.o *.mod $(TEST_DIR)/$(TEST)