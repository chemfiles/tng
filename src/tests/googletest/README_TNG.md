# Google Test

This directory contains the source code for the Google C++ Testing and Mocking Frameworks for embedding in the Trajectory New Generation library.

This directory is an almost exact copy of the googletest 1.10 source directory.
Unnesecary files have been removed including.

*.clang-format
*.gitignore
*.travis.yml
*appveyor.yml
*BUILD.bazel
*ci/
*CONTRIBUTING.md
*library.json
*platformio.ini
*WORKSPACE

This README_TNG.md file has been added.

Additionally the enable_testing() option has been commented out in CMakeLists.txt as we will not test googletest as part of the TNG build.

We have also disabled the installation of googletest in CMakeLists.txt, prefering to build only as we are embedding the code in our own project.
