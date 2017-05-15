# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/cmake/3.0.0/bin/cmake

# The command to remove a file.
RM = /usr/local/cmake/3.0.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lore/chanda5/phasta

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lore/chanda5/phasta/build

# Include any dependencies generated for this target.
include phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/depend.make

# Include the progress variables for this target.
include phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/progress.make

# Include the compile flags for this target's objects.
include phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/flags.make

phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o: phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/flags.make
phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o: ../phSolver/common/test/phIOwrite.f
	$(CMAKE_COMMAND) -E cmake_progress_report /lore/chanda5/phasta/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o"
	cd /lore/chanda5/phasta/build/phSolver/common/test && /usr/local/gcc/4.9.2/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /lore/chanda5/phasta/phSolver/common/test/phIOwrite.f -o CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o

phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o.requires:
.PHONY : phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o.requires

phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o.provides: phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o.requires
	$(MAKE) -f phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/build.make phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o.provides.build
.PHONY : phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o.provides

phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o.provides.build: phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o

# Object files for target phIOwriteFtn
phIOwriteFtn_OBJECTS = \
"CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o"

# External object files for target phIOwriteFtn
phIOwriteFtn_EXTERNAL_OBJECTS =

bin/phIOwriteFtn: phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o
bin/phIOwriteFtn: phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/build.make
bin/phIOwriteFtn: lib/libcommon.a
bin/phIOwriteFtn: lib/libphastaIO.a
bin/phIOwriteFtn: /usr/local/mpich3/3.1.2-thread-multiple/lib/libmpicxx.so
bin/phIOwriteFtn: /usr/local/mpich3/3.1.2-thread-multiple/lib/libmpifort.so
bin/phIOwriteFtn: /usr/local/mpich3/3.1.2-thread-multiple/lib/libmpi.so
bin/phIOwriteFtn: /usr/lib64/librt.so
bin/phIOwriteFtn: /usr/lib64/libpthread.so
bin/phIOwriteFtn: phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking Fortran executable ../../../bin/phIOwriteFtn"
	cd /lore/chanda5/phasta/build/phSolver/common/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/phIOwriteFtn.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/build: bin/phIOwriteFtn
.PHONY : phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/build

phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/requires: phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/phIOwrite.f.o.requires
.PHONY : phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/requires

phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/clean:
	cd /lore/chanda5/phasta/build/phSolver/common/test && $(CMAKE_COMMAND) -P CMakeFiles/phIOwriteFtn.dir/cmake_clean.cmake
.PHONY : phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/clean

phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/depend:
	cd /lore/chanda5/phasta/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lore/chanda5/phasta /lore/chanda5/phasta/phSolver/common/test /lore/chanda5/phasta/build /lore/chanda5/phasta/build/phSolver/common/test /lore/chanda5/phasta/build/phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : phSolver/common/test/CMakeFiles/phIOwriteFtn.dir/depend

