# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/local/cmake/3.8.1/bin/cmake

# The command to remove a file.
RM = /usr/local/cmake/3.8.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lore/chanda5/phasta

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lore/chanda5/phasta/build

# Include any dependencies generated for this target.
include phSolver/common/test/CMakeFiles/phIOreadFtn.dir/depend.make

# Include the progress variables for this target.
include phSolver/common/test/CMakeFiles/phIOreadFtn.dir/progress.make

# Include the compile flags for this target's objects.
include phSolver/common/test/CMakeFiles/phIOreadFtn.dir/flags.make

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/flags.make
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o: ../phSolver/common/test/chdir_mod.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lore/chanda5/phasta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o"
	cd /lore/chanda5/phasta/build/phSolver/common/test && /usr/local/gcc/4.9.2/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /lore/chanda5/phasta/phSolver/common/test/chdir_mod.f -o CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/phIOreadFtn.dir/chdir_mod.f.i"
	cd /lore/chanda5/phasta/build/phSolver/common/test && /usr/local/gcc/4.9.2/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /lore/chanda5/phasta/phSolver/common/test/chdir_mod.f > CMakeFiles/phIOreadFtn.dir/chdir_mod.f.i

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/phIOreadFtn.dir/chdir_mod.f.s"
	cd /lore/chanda5/phasta/build/phSolver/common/test && /usr/local/gcc/4.9.2/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /lore/chanda5/phasta/phSolver/common/test/chdir_mod.f -o CMakeFiles/phIOreadFtn.dir/chdir_mod.f.s

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.requires:

.PHONY : phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.requires

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.provides: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.requires
	$(MAKE) -f phSolver/common/test/CMakeFiles/phIOreadFtn.dir/build.make phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.provides.build
.PHONY : phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.provides

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.provides.build: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o


phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/flags.make
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o: ../phSolver/common/test/phIOread.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lore/chanda5/phasta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o"
	cd /lore/chanda5/phasta/build/phSolver/common/test && /usr/local/gcc/4.9.2/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /lore/chanda5/phasta/phSolver/common/test/phIOread.f -o CMakeFiles/phIOreadFtn.dir/phIOread.f.o

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/phIOreadFtn.dir/phIOread.f.i"
	cd /lore/chanda5/phasta/build/phSolver/common/test && /usr/local/gcc/4.9.2/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /lore/chanda5/phasta/phSolver/common/test/phIOread.f > CMakeFiles/phIOreadFtn.dir/phIOread.f.i

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/phIOreadFtn.dir/phIOread.f.s"
	cd /lore/chanda5/phasta/build/phSolver/common/test && /usr/local/gcc/4.9.2/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /lore/chanda5/phasta/phSolver/common/test/phIOread.f -o CMakeFiles/phIOreadFtn.dir/phIOread.f.s

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o.requires:

.PHONY : phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o.requires

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o.provides: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o.requires
	$(MAKE) -f phSolver/common/test/CMakeFiles/phIOreadFtn.dir/build.make phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o.provides.build
.PHONY : phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o.provides

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o.provides.build: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o


# Object files for target phIOreadFtn
phIOreadFtn_OBJECTS = \
"CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o" \
"CMakeFiles/phIOreadFtn.dir/phIOread.f.o"

# External object files for target phIOreadFtn
phIOreadFtn_EXTERNAL_OBJECTS =

bin/phIOreadFtn: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o
bin/phIOreadFtn: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o
bin/phIOreadFtn: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/build.make
bin/phIOreadFtn: lib/libcommon.a
bin/phIOreadFtn: lib/libphastaIO.a
bin/phIOreadFtn: /usr/local/mpich3/3.1.2-thread-multiple/lib/libmpicxx.so
bin/phIOreadFtn: /usr/local/mpich3/3.1.2-thread-multiple/lib/libmpifort.so
bin/phIOreadFtn: /usr/local/mpich3/3.1.2-thread-multiple/lib/libmpi.so
bin/phIOreadFtn: /usr/lib64/librt.so
bin/phIOreadFtn: /usr/lib64/libpthread.so
bin/phIOreadFtn: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/lore/chanda5/phasta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking Fortran executable ../../../bin/phIOreadFtn"
	cd /lore/chanda5/phasta/build/phSolver/common/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/phIOreadFtn.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/build: bin/phIOreadFtn

.PHONY : phSolver/common/test/CMakeFiles/phIOreadFtn.dir/build

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/requires: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.requires
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/requires: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o.requires

.PHONY : phSolver/common/test/CMakeFiles/phIOreadFtn.dir/requires

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/clean:
	cd /lore/chanda5/phasta/build/phSolver/common/test && $(CMAKE_COMMAND) -P CMakeFiles/phIOreadFtn.dir/cmake_clean.cmake
.PHONY : phSolver/common/test/CMakeFiles/phIOreadFtn.dir/clean

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/depend:
	cd /lore/chanda5/phasta/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lore/chanda5/phasta /lore/chanda5/phasta/phSolver/common/test /lore/chanda5/phasta/build /lore/chanda5/phasta/build/phSolver/common/test /lore/chanda5/phasta/build/phSolver/common/test/CMakeFiles/phIOreadFtn.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : phSolver/common/test/CMakeFiles/phIOreadFtn.dir/depend
