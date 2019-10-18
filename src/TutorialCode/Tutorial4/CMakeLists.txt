project(tutorial4)

if (Qt5_FOUND)
qt5_wrap_cpp(${PROJECT_NAME}_HEADERS_MOC ${PROJECT_NAME}.h)
elseif(Qt4_FOUND OR QT4_FOUND)

qt4_wrap_cpp(${PROJECT_NAME}_HEADERS_MOC ${PROJECT_NAME}.h)
endif()


add_library(${PROJECT_NAME} SHARED ${PROJECT_NAME}.cpp ${PROJECT_NAME}.h ${${PROJECT_NAME}_HEADERS_MOC})
target_link_libraries(${PROJECT_NAME} vleaf)
set_target_properties(${PROJECT_NAME}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "$<TARGET_FILE_DIR:VirtualLeaf>/../models"
    LIBRARY_OUTPUT_DIRECTORY "$<TARGET_FILE_DIR:VirtualLeaf>/../models"
    RUNTIME_OUTPUT_DIRECTORY "$<TARGET_FILE_DIR:VirtualLeaf>/../models"
)

QT_BIND_TO_TARGET(${PROJECT_NAME})