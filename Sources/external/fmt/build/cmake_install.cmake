# Install script for directory: /mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt" TYPE FILE FILES
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/fmt-config.cmake"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/fmt-config-version.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt/fmt-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt/fmt-targets.cmake"
         "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/CMakeFiles/Export/lib/cmake/fmt/fmt-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt/fmt-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt/fmt-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt" TYPE FILE FILES "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/CMakeFiles/Export/lib/cmake/fmt/fmt-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt" TYPE FILE FILES "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/CMakeFiles/Export/lib/cmake/fmt/fmt-targets-release.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/libfmt.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE OPTIONAL FILES
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/$<TARGET_PDB_FILE:fmt"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/fmt-header-only>"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/fmt" TYPE FILE FILES
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/chrono.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/color.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/compile.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/core.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/format.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/format-inl.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/locale.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/os.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/ostream.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/posix.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/printf.h"
    "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/include/fmt/ranges.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/fmt.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/doc/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/test/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/mnt/c/dev/Evo-Treatment-Plan/Sources/external/fmt/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
