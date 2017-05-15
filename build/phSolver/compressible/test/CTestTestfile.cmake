# CMake generated Testfile for 
# Source directory: /lore/chanda5/phasta/phSolver/compressible/test
# Build directory: /lore/chanda5/phasta/build/phSolver/compressible/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(compressible__inpCfg "cp" "/lore/chanda5/phasta/phSolver/common/input.config" "/path/to/phastaCases//compressible")
set_tests_properties(compressible__inpCfg PROPERTIES  LABELS "phsolver_compressible" WORKING_DIRECTORY "/path/to/phastaCases//compressible")
add_test(compressible_native_solverInp "ln" "-snf" "/path/to/phastaCases//compressible/solver.inp.native" "/path/to/phastaCases//compressible/solver.inp")
set_tests_properties(compressible_native_solverInp PROPERTIES  LABELS "phsolver_compressible" WORKING_DIRECTORY "/path/to/phastaCases//compressible")
add_test(compressible_native_inpCfg_staticBubble "cp" "/lore/chanda5/phasta/phSolver/common/input.config" "/path/to/phastaCases//staticBubble/run")
set_tests_properties(compressible_native_inpCfg_staticBubble PROPERTIES  LABELS "phsolver_compressible" WORKING_DIRECTORY "/path/to/phastaCases//staticBubble/run")
add_test(compressible_native_resCfg_staticBubble "cp" "-r" "/path/to/phastaCases//staticBubble/run/../chef/8-procs_case" "/path/to/phastaCases//staticBubble/run")
set_tests_properties(compressible_native_resCfg_staticBubble PROPERTIES  LABELS "phsolver_compressible" WORKING_DIRECTORY "/path/to/phastaCases//staticBubble/run")
add_test(compressible_native_staticBubble "mpirun" "-np" "8" "/lore/chanda5/phasta/build/bin/phastaC.exe")
set_tests_properties(compressible_native_staticBubble PROPERTIES  LABELS "phsolver_compressible" WORKING_DIRECTORY "/path/to/phastaCases//staticBubble/run")
add_test(compressible_native_compareRestart-staticBubble "mpirun" "-np" "8" "/lore/chanda5/phasta/build/bin/checkphasta" "/path/to/phastaCases//staticBubble/run/8-procs_case/" "/path/to/phastaCases//staticBubble/run/8-procs_case_ref/" "8" "1e-6")
set_tests_properties(compressible_native_compareRestart-staticBubble PROPERTIES  LABELS "phsolver_compressible" WORKING_DIRECTORY "/path/to/phastaCases//staticBubble/run")
