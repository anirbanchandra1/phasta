# CMake generated Testfile for 
# Source directory: /lore/chanda5/phasta/phSolver/common/test
# Build directory: /lore/chanda5/phasta/build/phSolver/common/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(common_readHeader "mpirun" "-np" "4" "/lore/chanda5/phasta/build/bin/phIOreadheader" "2")
set_tests_properties(common_readHeader PROPERTIES  LABELS "phsolver_common" WORKING_DIRECTORY "/path/to/phastaCases//incompressible")
add_test(common_readIlwork "mpirun" "-np" "1" "/lore/chanda5/phasta/build/bin/phIOreadIlwork" "geombc.dat.1")
set_tests_properties(common_readIlwork PROPERTIES  LABELS "phsolver_common" WORKING_DIRECTORY "/path/to/phastaCases//crossflow/4-1chef/4-procs_case")
add_test(common_readHeaderMultiTopo "mpirun" "-np" "4" "/lore/chanda5/phasta/build/bin/phIOposixMultiTopo")
set_tests_properties(common_readHeaderMultiTopo PROPERTIES  LABELS "phsolver_common" WORKING_DIRECTORY "/path/to/phastaCases//crossflow/4-1chef/4-procs_case")
add_test(common_readDatablock "mpirun" "-np" "4" "/lore/chanda5/phasta/build/bin/phIOreaddatablock" "2")
set_tests_properties(common_readDatablock PROPERTIES  LABELS "phsolver_common" WORKING_DIRECTORY "/path/to/phastaCases//incompressible")
add_test(common_write "mpirun" "-np" "4" "/lore/chanda5/phasta/build/bin/phIOwrite" "2")
set_tests_properties(common_write PROPERTIES  LABELS "phsolver_common" WORKING_DIRECTORY "/path/to/phastaCases/")
add_test(common_readFtn "mpirun" "-np" "4" "/lore/chanda5/phasta/build/bin/phIOreadFtn")
set_tests_properties(common_readFtn PROPERTIES  LABELS "phsolver_common" WORKING_DIRECTORY "/path/to/phastaCases//incompressible/")
add_test(common_writeFtn "mpirun" "-np" "4" "/lore/chanda5/phasta/build/bin/phIOwriteFtn")
set_tests_properties(common_writeFtn PROPERTIES  LABELS "phsolver_common" WORKING_DIRECTORY "/path/to/phastaCases/")
