# The set of languages for which implicit dependencies are needed:
set(CMAKE_DEPENDS_LANGUAGES
  "Fortran"
  )
# The set of files for implicit dependencies of each language:
set(CMAKE_DEPENDS_CHECK_Fortran
  "/lore/chanda5/phasta/phSolver/common/test/chdir_mod.f" "/lore/chanda5/phasta/build/phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o"
  "/lore/chanda5/phasta/phSolver/common/test/phIOread.f" "/lore/chanda5/phasta/build/phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o"
  )
set(CMAKE_Fortran_COMPILER_ID "GNU")

# Preprocessor definitions for this target.
set(CMAKE_TARGET_DEFINITIONS
  "DEBUG"
  "LINUX"
  "MPI"
  "MPICH_SKIP_MPICXX"
  "OMPI_SKIP_MPICXX"
  "OMPI_SKIP_MPICXX=1"
  "PARALLEL"
  )

# Targets to which this target links.
set(CMAKE_TARGET_LINKED_INFO_FILES
  "/lore/chanda5/phasta/build/phSolver/common/CMakeFiles/common.dir/DependInfo.cmake"
  "/lore/chanda5/phasta/build/phastaIO/CMakeFiles/phastaIO.dir/DependInfo.cmake"
  )

# Fortran module output directory.
set(CMAKE_Fortran_TARGET_MODULE_DIR "/lore/chanda5/phasta/build/phSolver/modules")

# The include file search paths:
set(CMAKE_C_TARGET_INCLUDE_PATH
  "."
  "../phastaIO"
  "../shapeFunction/src"
  "../phSolver/common"
  "phSolver/modules"
  "/usr/local/mpich3/3.1.2-thread-multiple/include"
  "../phSolver/common/phstreamEmpty"
  )
set(CMAKE_CXX_TARGET_INCLUDE_PATH ${CMAKE_C_TARGET_INCLUDE_PATH})
set(CMAKE_Fortran_TARGET_INCLUDE_PATH ${CMAKE_C_TARGET_INCLUDE_PATH})
set(CMAKE_ASM_TARGET_INCLUDE_PATH ${CMAKE_C_TARGET_INCLUDE_PATH})
