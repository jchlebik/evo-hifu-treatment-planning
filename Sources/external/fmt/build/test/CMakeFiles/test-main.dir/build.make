# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build

# Include any dependencies generated for this target.
include test/CMakeFiles/test-main.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/test-main.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/test-main.dir/flags.make

test/CMakeFiles/test-main.dir/test-main.cc.o: test/CMakeFiles/test-main.dir/flags.make
test/CMakeFiles/test-main.dir/test-main.cc.o: ../test/test-main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/test-main.dir/test-main.cc.o"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test-main.dir/test-main.cc.o -c /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test/test-main.cc

test/CMakeFiles/test-main.dir/test-main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test-main.dir/test-main.cc.i"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test/test-main.cc > CMakeFiles/test-main.dir/test-main.cc.i

test/CMakeFiles/test-main.dir/test-main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test-main.dir/test-main.cc.s"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test/test-main.cc -o CMakeFiles/test-main.dir/test-main.cc.s

test/CMakeFiles/test-main.dir/gtest-extra.cc.o: test/CMakeFiles/test-main.dir/flags.make
test/CMakeFiles/test-main.dir/gtest-extra.cc.o: ../test/gtest-extra.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object test/CMakeFiles/test-main.dir/gtest-extra.cc.o"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test-main.dir/gtest-extra.cc.o -c /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test/gtest-extra.cc

test/CMakeFiles/test-main.dir/gtest-extra.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test-main.dir/gtest-extra.cc.i"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test/gtest-extra.cc > CMakeFiles/test-main.dir/gtest-extra.cc.i

test/CMakeFiles/test-main.dir/gtest-extra.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test-main.dir/gtest-extra.cc.s"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test/gtest-extra.cc -o CMakeFiles/test-main.dir/gtest-extra.cc.s

test/CMakeFiles/test-main.dir/util.cc.o: test/CMakeFiles/test-main.dir/flags.make
test/CMakeFiles/test-main.dir/util.cc.o: ../test/util.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object test/CMakeFiles/test-main.dir/util.cc.o"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test-main.dir/util.cc.o -c /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test/util.cc

test/CMakeFiles/test-main.dir/util.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test-main.dir/util.cc.i"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test/util.cc > CMakeFiles/test-main.dir/util.cc.i

test/CMakeFiles/test-main.dir/util.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test-main.dir/util.cc.s"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test/util.cc -o CMakeFiles/test-main.dir/util.cc.s

# Object files for target test-main
test__main_OBJECTS = \
"CMakeFiles/test-main.dir/test-main.cc.o" \
"CMakeFiles/test-main.dir/gtest-extra.cc.o" \
"CMakeFiles/test-main.dir/util.cc.o"

# External object files for target test-main
test__main_EXTERNAL_OBJECTS =

test/libtest-main.a: test/CMakeFiles/test-main.dir/test-main.cc.o
test/libtest-main.a: test/CMakeFiles/test-main.dir/gtest-extra.cc.o
test/libtest-main.a: test/CMakeFiles/test-main.dir/util.cc.o
test/libtest-main.a: test/CMakeFiles/test-main.dir/build.make
test/libtest-main.a: test/CMakeFiles/test-main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libtest-main.a"
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && $(CMAKE_COMMAND) -P CMakeFiles/test-main.dir/cmake_clean_target.cmake
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test-main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/test-main.dir/build: test/libtest-main.a

.PHONY : test/CMakeFiles/test-main.dir/build

test/CMakeFiles/test-main.dir/clean:
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test && $(CMAKE_COMMAND) -P CMakeFiles/test-main.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/test-main.dir/clean

test/CMakeFiles/test-main.dir/depend:
	cd /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/test /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test/CMakeFiles/test-main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/test-main.dir/depend

