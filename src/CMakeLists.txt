####################################################################
#
# CMake Build Script for the  VirtualLeaf
#
# Frank T. Bergmann (fbergman@caltech.edu)
#


if (PREFIX_PATH)
set (CMAKE_PREFIX_PATH "${PREFIX_PATH}") 
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/ChooseQtVersion.cmake)
QT_FIND_MODULES(Core Gui OpenGL 3Support) 

file(GLOB CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp )
file(GLOB   H_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h )
file(GLOB  UI_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.ui )
file(GLOB RCC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.qrc )

list(REMOVE_ITEM CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/apoplastitem.cpp)
list(REMOVE_ITEM CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/rseed.cpp)

if (NOT WITH_X11)
list(REMOVE_ITEM CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/x11graph.cpp)
endif()

if (NOT WITH_QWT)
list(REMOVE_ITEM CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/data_plot.cpp)
list(REMOVE_ITEM H_FILES ${CMAKE_CURRENT_SOURCE_DIR}/data_plot.h)
endif()

if (Qt5_FOUND)
qt5_add_resources(RESOURCES_RCC ${RCC_FILES})
qt5_wrap_cpp(HEADERS_MOC ${H_FILES})
qt5_wrap_ui(FORMS_HEADERS ${UI_FILES})
elseif(Qt4_FOUND OR QT4_FOUND)
qt4_add_resources(RESOURCES_RCC ${RCC_FILES})
qt4_wrap_cpp(HEADERS_MOC ${H_FILES})
qt4_wrap_ui(FORMS_HEADERS ${UI_FILES})
endif()

#include(${QT_USE_FILE})
#add_definitions(${QT_DEFINITIONS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
if (EXTRA_INCLUDE)
include_directories(BEFORE ${EXTRA_INCLUDE})
endif()
include_directories(${LIBXML2_INCLUDE_DIR})
include_directories(${ZLIB_INCLUDE_DIRS})

add_definitions(-DQTGRAPHICS)

set(EXECUTABLE_FLAGS)

if(WIN32 AND NOT CYGWIN)
  option(ENABLE_UI_CONSOLE "Show the console window (to see debug messages) when running Spatial UI" OFF)
  if (NOT ENABLE_UI_CONSOLE)
    set(EXECUTABLE_FLAGS WIN32)
    if(MSVC)
      set(CPP_FILES ${CPP_FILES} VirtualLeaf.rc)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:mainCRTStartup ")
    endif(MSVC)
  endif(NOT ENABLE_UI_CONSOLE)
endif(WIN32 AND NOT CYGWIN)

if (APPLE)
 
set(EXECUTABLE_FLAGS MACOSX_BUNDLE)
set(MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/SpatialUI.plist )
set(MACOSX_BUNDLE_ICON_FILE leaficon.icns )
set(MACOSX_BUNDLE_BUNDLE_VERSION  "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(MACOSX_BUNDLE_INFO_STRING "VirtualLeaf ${MACOSX_BUNDLE_BUNDLE_VERSION} is distributed under the GPL." )
set(MACOSX_BUNDLE_BUNDLE_NAME "VirtualLeaf")
set(MACOSX_BUNDLE_LONG_VERSION_STRING "${MACOSX_BUNDLE_BUNDLE_VERSION}")
set(MACOSX_BUNDLE_SHORT_VERSION_STRING  "${MACOSX_BUNDLE_BUNDLE_VERSION}")
set(MACOSX_BUNDLE_COPYRIGHT ${MACOSX_BUNDLE_INFO_STRING})

endif(APPLE)


add_executable(VirtualLeaf ${EXECUTABLE_FLAGS} ${CPP_FILES}
  ${H_FILES}
  ${HEADERS_MOC}
  ${FORMS_HEADERS}
  ${RESOURCES_RCC}
  ${VLEAF_CPP_FILES}
)

QT_BIND_TO_TARGET(VirtualLeaf)

target_link_libraries(VirtualLeaf ${EXTRA_LIBS} ${LIBXML2_LIBRARIES} ${ZLIB_LIBRARIES})


  if (WIN32)
    target_link_libraries(VirtualLeaf Ws2_32.lib)
  endif()
  
set (VLEAF_BIN_DIR $<TARGET_FILE_DIR:VirtualLeaf>)


set(VLEAF_API
cellbase.cpp
 matrix.cpp
 output.cpp
 parameter.cpp
 parse.cpp
 random.cpp
 simplugin.cpp
 UniqueMessage.cpp
 vector.cpp
 wallbase.cpp
 warning.cpp
)

set(API_H_FILES
 cellbase.h
 matrix.h
 output.h
 parameter.h
 parse.h
 random.h
 simplugin.h
 UniqueMessage.h
 vector.h
 wallbase.h
 warning.h
)

if (Qt5_FOUND)
qt5_wrap_cpp(API_HEADERS_MOC ${API_H_FILES})
elseif(Qt4_FOUND OR QT4_FOUND)
qt4_wrap_cpp(API_HEADERS_MOC ${API_H_FILES})
endif()

add_library(vleaf ${VLEAF_API} ${H_FILES} ${API_HEADERS_MOC})
target_link_libraries(vleaf ${LIBXML2_LIBRARIES} ${ZLIB_LIBRARIES})
QT_BIND_TO_TARGET(vleaf)


SET(DIRS)

if (EXTRA_LIBS)
foreach(${lib} ${EXTRA_LIBS})
  if (EXISTS ${lib})
	  get_filename_component(dir ${lib} PATH)	
	  SET(DIRS  ${DIRS} ${dir})
  endif()
endforeach()
endif()

if (NOT DIRS)
   if (EXISTS "${EXTRA_LIBS}")
	   get_filename_component(DIRS ${EXTRA_LIBS} PATH)
   endif()
endif()

if (QT_BINARY_DIR)
  SET(DIRS ${QT_BINARY_DIR} ${DIRS})
endif()


if(NOT APPLE)

  SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/VirtualLeaf${CMAKE_EXECUTABLE_SUFFIX}" CACHE INTERNAL "")
  install(TARGETS VirtualLeaf DESTINATION bin)
  # ensure that next to the binaries also dependencies are copied
  install(CODE "
      include(BundleUtilities)
      fixup_bundle(\"${APPS}\" \"\" \"${DIRS}\")
      " COMPONENT Runtime)

else()

  install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/leaficon.icns DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/VirtualLeaf.app/Contents/Resources)
  # on OS X everything is different, here we need to first fix the app bundle,
  # and only then install it
  SET(APPS "${CMAKE_CURRENT_BINARY_DIR}/VirtualLeaf.app" CACHE INTERNAL "")
  # ensure that next to the binaries also dependencies are copied
  install(CODE "
      include(BundleUtilities)
      fixup_bundle(\"${APPS}\" \"\" \"${DIRS}\")
      " COMPONENT Runtime)
  install(TARGETS VirtualLeaf DESTINATION .)

endif()

include(InstallRequiredSystemLibraries)


add_subdirectory(TutorialCode/Wortel)