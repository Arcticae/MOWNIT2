# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.12

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2018.2.4\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2018.2.4\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/lab5_c.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/lab5_c.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lab5_c.dir/flags.make

CMakeFiles/lab5_c.dir/main.cpp.obj: CMakeFiles/lab5_c.dir/flags.make
CMakeFiles/lab5_c.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/lab5_c.dir/main.cpp.obj"
	C:\PROGRA~1\HASKEL~1\826561~1.1\mingw\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\lab5_c.dir\main.cpp.obj -c "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c\main.cpp"

CMakeFiles/lab5_c.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab5_c.dir/main.cpp.i"
	C:\PROGRA~1\HASKEL~1\826561~1.1\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c\main.cpp" > CMakeFiles\lab5_c.dir\main.cpp.i

CMakeFiles/lab5_c.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab5_c.dir/main.cpp.s"
	C:\PROGRA~1\HASKEL~1\826561~1.1\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c\main.cpp" -o CMakeFiles\lab5_c.dir\main.cpp.s

# Object files for target lab5_c
lab5_c_OBJECTS = \
"CMakeFiles/lab5_c.dir/main.cpp.obj"

# External object files for target lab5_c
lab5_c_EXTERNAL_OBJECTS =

lab5_c.exe: CMakeFiles/lab5_c.dir/main.cpp.obj
lab5_c.exe: CMakeFiles/lab5_c.dir/build.make
lab5_c.exe: CMakeFiles/lab5_c.dir/linklibs.rsp
lab5_c.exe: CMakeFiles/lab5_c.dir/objects1.rsp
lab5_c.exe: CMakeFiles/lab5_c.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable lab5_c.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\lab5_c.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lab5_c.dir/build: lab5_c.exe

.PHONY : CMakeFiles/lab5_c.dir/build

CMakeFiles/lab5_c.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\lab5_c.dir\cmake_clean.cmake
.PHONY : CMakeFiles/lab5_c.dir/clean

CMakeFiles/lab5_c.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c" "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c" "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c\cmake-build-debug" "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c\cmake-build-debug" "C:\Users\Timelock\Desktop\Studia Vol5\MOWNIT2\lab5_c\cmake-build-debug\CMakeFiles\lab5_c.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/lab5_c.dir/depend

