####################################################################
#
# CMake Build Script for VirtualLeaf
#
# Frank T. Bergmann
#

cmake_minimum_required(VERSION 3.10)
project(VirtualLeaf)


set(VLEAF_VERSION_MAJOR)
set(VLEAF_VERSION_MINOR)
set(VLEAF_VERSION_PATCH)
set(VLEAF_VERSION_RELEASE)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION.txt")

    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION.txt" VersionString NEWLINE_CONSUME)
    string(STRIP "${VersionString}" VersionString)
    string(REPLACE "." ";" VersionString "${VersionString}" )
    string(REPLACE "-" ";" VersionString "${VersionString}" )
    list(LENGTH VersionString versionLength)
    list(GET VersionString 0 VLEAF_VERSION_MAJOR )
    list(GET VersionString 1 VLEAF_VERSION_MINOR )
    list(GET VersionString 2 VLEAF_VERSION_PATCH )

    if(${versionLength} GREATER 3)
        list(GET VersionString 3 VLEAF_VERSION_RELEASE )
    endif()

endif()

set(CPACK_PACKAGE_VERSION_MAJOR "${VLEAF_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${VLEAF_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${VLEAF_VERSION_PATCH}")
set(CPACK_SOURCE_GENERATOR "ZIP")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES
  "/.vs;/.vscode;/build*;screenshots;.DS_STORE;testmodels/.DS_Store;/build/;/.bzr/;~$;${CPACK_SOURCE_IGNORE_FILES}")
set(CPACK_NSIS_CREATE_ICONS "CreateShortCut '\$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\VirtualLeaf.lnk' '\$INSTDIR\\\\bin\\\\VirtualLeaf.exe'")
include(CPack)

find_package(LibXml2 REQUIRED)
find_package(ZLIB REQUIRED)

option(WITH_QWT "Use QWT." OFF)
if (WITH_QWT)
add_definitions(-DHAVE_QWT)
endif()

option(WITH_QT "Use QT." ON)
if (WITH_QT)
add_definitions(-DQTGRAPHICS)
endif()

option(WITH_LEAF_FLEMING "leaf_fleming" OFF)
if (WITH_LEAF_FLEMING)
add_definitions(-DFLEMING)
endif()

option(WITH_XFIG "xfig" OFF)
if (WITH_XFIG)
add_definitions(-DXFIGGRAPHICS)
endif()

option(WITH_X11 "x11 graphics" OFF)
if (WITH_X11)
add_definitions(-DX11GRAPHICS)
endif()

###############################################################################
#
# Need some variables set up, such as the name for the libSBML
# library and the Path and file separator characters. The
# MISC_PREFIX variable will cause libsbml bindings, examples and
# documentation to be installed in PREFIX/${MISC_PREFIX}.
#

if(UNIX OR CYGWIN)
    set(PATH_SEP "/")
    set(FILE_SEP ":")
else()
    set(PATH_SEP "\\")
    set(FILE_SEP ";")
endif()

set(EXTRA_LIBS "" CACHE STRING "List of additional libraries to link against." )
set(EXTRA_FLAGS "" CACHE STRING "List of extra flags." )
set(EXTRA_INCLUDE "" CACHE STRING "List of include diretories." )

if (EXTRA_FLAGS)
foreach(${flag} ${EXTRA_FLAGS})
	add_definitions(-D${flag})
endforeach()
endif()


if (APPLE)
add_definitions(-DDARWIN)
endif()

if (UNIX)
add_definitions(-DLINUX)
endif()



if(WIN32 AND NOT CYGWIN)
 add_definitions(-DLIBSBML_EXPORTS -DWIN32)
 if(MSVC)
        add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
        set(BUILD_DEFINITIONS "${BUILD_DEFINITIONS} -D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS")
        option(WITH_STATIC_RUNTIME "Compile using the static MSVC Runtime." OFF)
        if(WITH_STATIC_RUNTIME)
            foreach(flag_var
                CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
                CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO
                CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
                CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)

                if(${flag_var} MATCHES "/MD")
                    string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
                endif(${flag_var} MATCHES "/MD")
            endforeach(flag_var)
        add_definitions( -D_MT)
        endif(WITH_STATIC_RUNTIME)
  endif()
endif()

add_subdirectory(src)

message(STATUS         
"
 Project Virtual Leaf ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}
 ==========================
 
 Install Prefix : ${CMAKE_INSTALL_PREFIX}
                
 libxml2 include dir : ${LIBXML2_INCLUDE_DIR}
 libxml2 libs        : ${LIBXML2_LIBRARIES}
 libxml2 version     : ${LIBXML2_VERSION_STRING}
 
 zlib include dir : ${ZLIB_INCLUDE_DIRS}
 zlib libs        : ${ZLIB_LIBRARIES}
 
 C   FLAGS      : ${CMAKE_C_FLAGS}
 CXX FLAGS      : ${CMAKE_CXX_FLAGS}
 Extra Libs     : ${EXTRA_LIBS} 
 
"
)
