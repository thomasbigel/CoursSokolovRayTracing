# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = D:/VisualStudioC/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe

# The command to remove a file.
RM = D:/VisualStudioC/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing/build"

# Include any dependencies generated for this target.
include CMakeFiles/tinyraytracer.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/tinyraytracer.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/tinyraytracer.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tinyraytracer.dir/flags.make

CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.obj: CMakeFiles/tinyraytracer.dir/flags.make
CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.obj: ../tinyraytracer.cpp
CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.obj: CMakeFiles/tinyraytracer.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="D:/M1/Representation visuel de donnees/CoursSokolovRayTracing/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.obj"
	C:/msys64/mingw64/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.obj -MF CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.obj.d -o CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.obj -c "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing/tinyraytracer.cpp"

CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.i"
	C:/msys64/mingw64/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing/tinyraytracer.cpp" > CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.i

CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.s"
	C:/msys64/mingw64/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing/tinyraytracer.cpp" -o CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.s

# Object files for target tinyraytracer
tinyraytracer_OBJECTS = \
"CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.obj"

# External object files for target tinyraytracer
tinyraytracer_EXTERNAL_OBJECTS =

tinyraytracer.exe: CMakeFiles/tinyraytracer.dir/tinyraytracer.cpp.obj
tinyraytracer.exe: CMakeFiles/tinyraytracer.dir/build.make
tinyraytracer.exe: CMakeFiles/tinyraytracer.dir/linklibs.rsp
tinyraytracer.exe: CMakeFiles/tinyraytracer.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="D:/M1/Representation visuel de donnees/CoursSokolovRayTracing/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable tinyraytracer.exe"
	D:/VisualStudioC/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe -E rm -f CMakeFiles/tinyraytracer.dir/objects.a
	C:/msys64/mingw64/bin/ar.exe qc CMakeFiles/tinyraytracer.dir/objects.a @CMakeFiles/tinyraytracer.dir/objects1.rsp
	C:/msys64/mingw64/bin/c++.exe  -Wall -Wextra -pedantic -std=c++11 -O3 -fopenmp -g -Wl,--whole-archive CMakeFiles/tinyraytracer.dir/objects.a -Wl,--no-whole-archive -o tinyraytracer.exe -Wl,--out-implib,libtinyraytracer.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/tinyraytracer.dir/linklibs.rsp

# Rule to build all files generated by this target.
CMakeFiles/tinyraytracer.dir/build: tinyraytracer.exe
.PHONY : CMakeFiles/tinyraytracer.dir/build

CMakeFiles/tinyraytracer.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tinyraytracer.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tinyraytracer.dir/clean

CMakeFiles/tinyraytracer.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing" "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing" "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing/build" "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing/build" "D:/M1/Representation visuel de donnees/CoursSokolovRayTracing/build/CMakeFiles/tinyraytracer.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/tinyraytracer.dir/depend

