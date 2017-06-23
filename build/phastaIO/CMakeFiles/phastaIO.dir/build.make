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
include phastaIO/CMakeFiles/phastaIO.dir/depend.make

# Include the progress variables for this target.
include phastaIO/CMakeFiles/phastaIO.dir/progress.make

# Include the compile flags for this target's objects.
include phastaIO/CMakeFiles/phastaIO.dir/flags.make

phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o: phastaIO/CMakeFiles/phastaIO.dir/flags.make
phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o: ../phastaIO/phiotmrc.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lore/chanda5/phasta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o"
	cd /lore/chanda5/phasta/build/phastaIO && /usr/local/gcc/4.9.2/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/phastaIO.dir/phiotmrc.cc.o -c /lore/chanda5/phasta/phastaIO/phiotmrc.cc

phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/phastaIO.dir/phiotmrc.cc.i"
	cd /lore/chanda5/phasta/build/phastaIO && /usr/local/gcc/4.9.2/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /lore/chanda5/phasta/phastaIO/phiotmrc.cc > CMakeFiles/phastaIO.dir/phiotmrc.cc.i

phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/phastaIO.dir/phiotmrc.cc.s"
	cd /lore/chanda5/phasta/build/phastaIO && /usr/local/gcc/4.9.2/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /lore/chanda5/phasta/phastaIO/phiotmrc.cc -o CMakeFiles/phastaIO.dir/phiotmrc.cc.s

phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o.requires:

.PHONY : phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o.requires

phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o.provides: phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o.requires
	$(MAKE) -f phastaIO/CMakeFiles/phastaIO.dir/build.make phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o.provides.build
.PHONY : phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o.provides

phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o.provides.build: phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o


phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o: phastaIO/CMakeFiles/phastaIO.dir/flags.make
phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o: ../phastaIO/phastaIO.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lore/chanda5/phasta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o"
	cd /lore/chanda5/phasta/build/phastaIO && /usr/local/gcc/4.9.2/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/phastaIO.dir/phastaIO.cc.o -c /lore/chanda5/phasta/phastaIO/phastaIO.cc

phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/phastaIO.dir/phastaIO.cc.i"
	cd /lore/chanda5/phasta/build/phastaIO && /usr/local/gcc/4.9.2/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /lore/chanda5/phasta/phastaIO/phastaIO.cc > CMakeFiles/phastaIO.dir/phastaIO.cc.i

phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/phastaIO.dir/phastaIO.cc.s"
	cd /lore/chanda5/phasta/build/phastaIO && /usr/local/gcc/4.9.2/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /lore/chanda5/phasta/phastaIO/phastaIO.cc -o CMakeFiles/phastaIO.dir/phastaIO.cc.s

phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o.requires:

.PHONY : phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o.requires

phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o.provides: phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o.requires
	$(MAKE) -f phastaIO/CMakeFiles/phastaIO.dir/build.make phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o.provides.build
.PHONY : phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o.provides

phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o.provides.build: phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o


# Object files for target phastaIO
phastaIO_OBJECTS = \
"CMakeFiles/phastaIO.dir/phiotmrc.cc.o" \
"CMakeFiles/phastaIO.dir/phastaIO.cc.o"

# External object files for target phastaIO
phastaIO_EXTERNAL_OBJECTS =

lib/libphastaIO.a: phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o
lib/libphastaIO.a: phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o
lib/libphastaIO.a: phastaIO/CMakeFiles/phastaIO.dir/build.make
lib/libphastaIO.a: phastaIO/CMakeFiles/phastaIO.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/lore/chanda5/phasta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library ../lib/libphastaIO.a"
	cd /lore/chanda5/phasta/build/phastaIO && $(CMAKE_COMMAND) -P CMakeFiles/phastaIO.dir/cmake_clean_target.cmake
	cd /lore/chanda5/phasta/build/phastaIO && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/phastaIO.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
phastaIO/CMakeFiles/phastaIO.dir/build: lib/libphastaIO.a

.PHONY : phastaIO/CMakeFiles/phastaIO.dir/build

phastaIO/CMakeFiles/phastaIO.dir/requires: phastaIO/CMakeFiles/phastaIO.dir/phiotmrc.cc.o.requires
phastaIO/CMakeFiles/phastaIO.dir/requires: phastaIO/CMakeFiles/phastaIO.dir/phastaIO.cc.o.requires

.PHONY : phastaIO/CMakeFiles/phastaIO.dir/requires

phastaIO/CMakeFiles/phastaIO.dir/clean:
	cd /lore/chanda5/phasta/build/phastaIO && $(CMAKE_COMMAND) -P CMakeFiles/phastaIO.dir/cmake_clean.cmake
.PHONY : phastaIO/CMakeFiles/phastaIO.dir/clean

phastaIO/CMakeFiles/phastaIO.dir/depend:
	cd /lore/chanda5/phasta/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lore/chanda5/phasta /lore/chanda5/phasta/phastaIO /lore/chanda5/phasta/build /lore/chanda5/phasta/build/phastaIO /lore/chanda5/phasta/build/phastaIO/CMakeFiles/phastaIO.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : phastaIO/CMakeFiles/phastaIO.dir/depend

