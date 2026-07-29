
###############################################################################
#
#  EGSnrc egs++ makefile to build kerma scoring object
#  Copyright (C) 2026 National Research Council Canada
#
#  This file is part of EGSnrc.
#
#  EGSnrc is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Affero General Public License as published by the
#  Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  EGSnrc is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
#  more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with EGSnrc. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
#
#  Author:          Ernesto Mainegra-Hing, 2026
#
#  Contributors:
#
###############################################################################

include $(EGS_CONFIG)
include $(SPEC_DIR)egspp.spec
include $(SPEC_DIR)egspp_$(my_machine).conf

DEFS = $(DEF1) -DBUILD_KERMA_SCORING_DLL

library  = egs_kerma_scoring
lib_files = egs_kerma_scoring

# Add the fluence scoring header directory to the include path
INC2 = -I$(IEGS2) -I..$(DSEP).. -I..$(DSEP)egs_fluence_scoring

# Link against egspp and the fluence scoring shared library
link2_libs = egspp egs_fluence_scoring

my_deps  = $(common_ausgab_deps)
extra_dep = $(addprefix $(DSOLIBS), $(my_deps))

include $(SPEC_DIR)egspp_libs.spec

$(make_depend)

clean:
	$(REMOVE) $(ABS_DSO)*$(library)*
	$(REMOVE) *.o
