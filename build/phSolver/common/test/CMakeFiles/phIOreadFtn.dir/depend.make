# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0


phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.mod.proxy: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.provides
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod phSolver/modules/chdir_mod phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.provides.build
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/build: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.f.o.provides.build
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o: /usr/local/mpich3/3.1.2-thread-multiple/include/mpif.h

phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o.requires: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.mod.proxy
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o: phSolver/common/test/CMakeFiles/phIOreadFtn.dir/chdir_mod.mod.stamp
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o: phSolver/common/CMakeFiles/common.dir/phio.mod.stamp
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o: phSolver/common/CMakeFiles/common.dir/posixio.mod.stamp
phSolver/common/test/CMakeFiles/phIOreadFtn.dir/phIOread.f.o: phSolver/common/CMakeFiles/common.dir/syncio.mod.stamp
