# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.6.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.6.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/build

# Include any dependencies generated for this target.
include CMakeFiles/KSM_CalculationZip.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/KSM_CalculationZip.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/KSM_CalculationZip.dir/flags.make

CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o: CMakeFiles/KSM_CalculationZip.dir/flags.make
CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o: ../KSM_CalculationZip.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o -c /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/KSM_CalculationZip.cpp

CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/KSM_CalculationZip.cpp > CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.i

CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/KSM_CalculationZip.cpp -o CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.s

CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o.requires:

.PHONY : CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o.requires

CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o.provides: CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o.requires
	$(MAKE) -f CMakeFiles/KSM_CalculationZip.dir/build.make CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o.provides.build
.PHONY : CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o.provides

CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o.provides.build: CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o


CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o: CMakeFiles/KSM_CalculationZip.dir/flags.make
CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o: ../cmdline.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o   -c /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/cmdline.c

CMakeFiles/KSM_CalculationZip.dir/cmdline.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/KSM_CalculationZip.dir/cmdline.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/cmdline.c > CMakeFiles/KSM_CalculationZip.dir/cmdline.c.i

CMakeFiles/KSM_CalculationZip.dir/cmdline.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/KSM_CalculationZip.dir/cmdline.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/cmdline.c -o CMakeFiles/KSM_CalculationZip.dir/cmdline.c.s

CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o.requires:

.PHONY : CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o.requires

CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o.provides: CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o.requires
	$(MAKE) -f CMakeFiles/KSM_CalculationZip.dir/build.make CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o.provides.build
.PHONY : CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o.provides

CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o.provides.build: CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o


CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o: CMakeFiles/KSM_CalculationZip.dir/flags.make
CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o: ../jsoncpp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o -c /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/jsoncpp.cpp

CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/jsoncpp.cpp > CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.i

CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/jsoncpp.cpp -o CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.s

CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o.requires:

.PHONY : CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o.requires

CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o.provides: CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o.requires
	$(MAKE) -f CMakeFiles/KSM_CalculationZip.dir/build.make CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o.provides.build
.PHONY : CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o.provides

CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o.provides.build: CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o


# Object files for target KSM_CalculationZip
KSM_CalculationZip_OBJECTS = \
"CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o" \
"CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o" \
"CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o"

# External object files for target KSM_CalculationZip
KSM_CalculationZip_EXTERNAL_OBJECTS =

KSM_CalculationZip: CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o
KSM_CalculationZip: CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o
KSM_CalculationZip: CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o
KSM_CalculationZip: CMakeFiles/KSM_CalculationZip.dir/build.make
KSM_CalculationZip: CMakeFiles/KSM_CalculationZip.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable KSM_CalculationZip"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/KSM_CalculationZip.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/KSM_CalculationZip.dir/build: KSM_CalculationZip

.PHONY : CMakeFiles/KSM_CalculationZip.dir/build

CMakeFiles/KSM_CalculationZip.dir/requires: CMakeFiles/KSM_CalculationZip.dir/KSM_CalculationZip.cpp.o.requires
CMakeFiles/KSM_CalculationZip.dir/requires: CMakeFiles/KSM_CalculationZip.dir/cmdline.c.o.requires
CMakeFiles/KSM_CalculationZip.dir/requires: CMakeFiles/KSM_CalculationZip.dir/jsoncpp.cpp.o.requires

.PHONY : CMakeFiles/KSM_CalculationZip.dir/requires

CMakeFiles/KSM_CalculationZip.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/KSM_CalculationZip.dir/cmake_clean.cmake
.PHONY : CMakeFiles/KSM_CalculationZip.dir/clean

CMakeFiles/KSM_CalculationZip.dir/depend:
	cd /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/build /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/build /Users/cvardema/Workspace/Projects/cemodelservice/api/HAKOU/HAKOU/v3/source/src/build/CMakeFiles/KSM_CalculationZip.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/KSM_CalculationZip.dir/depend

