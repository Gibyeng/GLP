# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /home/changye/Documents/GLP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/changye/Documents/GLP

# Include any dependencies generated for this target.
include CMakeFiles/GPULP.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/GPULP.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/GPULP.dir/flags.make

CMakeFiles/GPULP.dir/GPULP_generated_main.cu.o: CMakeFiles/GPULP.dir/GPULP_generated_main.cu.o.depend
CMakeFiles/GPULP.dir/GPULP_generated_main.cu.o: CMakeFiles/GPULP.dir/GPULP_generated_main.cu.o.cmake
CMakeFiles/GPULP.dir/GPULP_generated_main.cu.o: main.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/changye/Documents/GLP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building NVCC (Device) object CMakeFiles/GPULP.dir/GPULP_generated_main.cu.o"
	cd /home/changye/Documents/GLP/CMakeFiles/GPULP.dir && /usr/bin/cmake -E make_directory /home/changye/Documents/GLP/CMakeFiles/GPULP.dir//.
	cd /home/changye/Documents/GLP/CMakeFiles/GPULP.dir && /usr/bin/cmake -D verbose:BOOL=$(VERBOSE) -D build_configuration:STRING= -D generated_file:STRING=/home/changye/Documents/GLP/CMakeFiles/GPULP.dir//./GPULP_generated_main.cu.o -D generated_cubin_file:STRING=/home/changye/Documents/GLP/CMakeFiles/GPULP.dir//./GPULP_generated_main.cu.o.cubin.txt -P /home/changye/Documents/GLP/CMakeFiles/GPULP.dir//GPULP_generated_main.cu.o.cmake

# Object files for target GPULP
GPULP_OBJECTS =

# External object files for target GPULP
GPULP_EXTERNAL_OBJECTS = \
"/home/changye/Documents/GLP/CMakeFiles/GPULP.dir/GPULP_generated_main.cu.o"

GPULP: CMakeFiles/GPULP.dir/GPULP_generated_main.cu.o
GPULP: CMakeFiles/GPULP.dir/build.make
GPULP: /usr/local/cuda/lib64/libcudart_static.a
GPULP: /usr/lib/x86_64-linux-gnu/librt.so
GPULP: libkernel.a
GPULP: /usr/local/cuda/lib64/libcudadevrt.a
GPULP: /usr/local/cuda/lib64/libcudart_static.a
GPULP: /usr/lib/x86_64-linux-gnu/librt.so
GPULP: /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
GPULP: /usr/lib/x86_64-linux-gnu/libpthread.so
GPULP: CMakeFiles/GPULP.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/changye/Documents/GLP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable GPULP"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GPULP.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/GPULP.dir/build: GPULP

.PHONY : CMakeFiles/GPULP.dir/build

CMakeFiles/GPULP.dir/requires:

.PHONY : CMakeFiles/GPULP.dir/requires

CMakeFiles/GPULP.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/GPULP.dir/cmake_clean.cmake
.PHONY : CMakeFiles/GPULP.dir/clean

CMakeFiles/GPULP.dir/depend: CMakeFiles/GPULP.dir/GPULP_generated_main.cu.o
	cd /home/changye/Documents/GLP && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/changye/Documents/GLP /home/changye/Documents/GLP /home/changye/Documents/GLP /home/changye/Documents/GLP /home/changye/Documents/GLP/CMakeFiles/GPULP.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/GPULP.dir/depend

