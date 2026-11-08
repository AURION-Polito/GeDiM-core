list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/FastMarchingMethod.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/BrentSegmentMinimization.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/BrentTriangleMinimization.cpp)

list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/FastMarchingMethod.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MinimizationResult.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/BrentSegmentMinimization.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/BrentTriangleMinimization.hpp)

list(APPEND gedim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gedim_src_sources ${gedim_src_sources} PARENT_SCOPE)
set(gedim_src_headers ${gedim_src_headers} PARENT_SCOPE)
set(gedim_src_includes ${gedim_src_includes} PARENT_SCOPE)
