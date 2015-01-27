# This is the makefile for the mesa atm library

MESA_DIR = /data/mesa7324

#################################################################
#
# PREAMBLE

include $(MESA_DIR)/utils/makefile_header

ifeq ($(CLR_OFFLOAD),YES)
CC += -Doffload
FC += -Doffload
endif

#################################################################
#
# SOURCE FILES

SRCS = color_def.f90 color_lib.f90 load_bc_tables.f90 eval_bc_tables.f90

#################################################################
#
# TARGETS

LIB = libcolor.a
LIB_OBJS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
$(LIB) : $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

#################################################################
#
# COMPILATION RULES

#COMPILE = $(COMPILE_TO_DEPLOY)
COMPILE = $(COMPILE_TO_TEST)



COMPILE_CMD = $(COMPILE) $(FCfree)

%.o : %.f90
	$(COMPILE_CMD) $<

%.mod : %.o
	@true

clean:
	-@rm -f *.o *.mod *.modmic *genmod.f90 *.a .depend

#################################################################
#
# DEPENDENCIES

# Paths to search for source & include files

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

# Automatic dependency generation using makedepf90

nodeps : $(.DEFAULT_GOAL)

NODEPS = $(or $(filter nodeps,$(MAKECMDGOALS)),$(filter clean,$(MAKECMDGOALS)))

ifeq ($(NODEPS),)

  DEPEND_DIRS = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

  .depend :
	makedepf90 -m %m.mod -I$(DEPEND_DIRS) $(SRCS) > .depend
  -include .depend

endif
