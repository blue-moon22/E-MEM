cmake_minimum_required (VERSION 3.4)
project (e-mem CXX)
 
  
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.71.0 COMPONENTS filesystem regex) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
add_executable(e-mem e-mem.cpp) 
target_link_libraries (e-mem seqan3::seqan3)
target_link_libraries(e-mem ${Boost_LIBRARIES})
