add_executable(ArborX_HDBSCAN.exe
  parallel_boruvka.cpp
  edge_hierarchy.cpp
  computeDelta.cpp 
  incidence_matrix.cpp 
  main.cpp
)

set_target_properties(ArborX_HDBSCAN.exe PROPERTIES
  CXX_STANDARD 14
  CXX_EXTENSIONS OFF
)

add_executable(ArborX_mst2cluster.exe
  parallel_boruvka.cpp
  edge_hierarchy.cpp
  computeDelta.cpp 
  incidence_matrix.cpp 
  mst2cluster.cpp
)


set_target_properties(ArborX_mst2cluster.exe PROPERTIES
  CXX_STANDARD 14
  CXX_EXTENSIONS OFF
)
add_test(NAME ArborX_HDBSCAN_Example COMMAND ./ArborX_HDBSCAN.exe ${CMAKE_CURRENT_SOURCE_DIR}/notebooks/twoD-20.pts)