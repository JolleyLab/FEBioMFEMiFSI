# Tests organized into subdirectories

# Subdirectories
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

foreach (subdir ${SUBDIRS})
   add_subdirectory(${subdir})
endforeach()
