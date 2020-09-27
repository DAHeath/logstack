cmake_minimum_required (VERSION 3.10)
project (LogSGC)
set(NAME "LogSGC")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


add_definitions(-Wfatal-errors)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-pthread -Wall -march=native -O3")

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

find_package(Catch2 REQUIRED)
find_package(OpenSSL REQUIRED)
if(OPENSSL_FOUND)
  if(APPLE)
    include_directories(/usr/local/Cellar/openssl@1.1/1.1.1d/include)
    list(APPEND LIB_LIST /usr/local/Cellar/openssl@1.1/1.1.1d/lib/libssl.dylib)
    list(APPEND LIB_LIST /usr/local/Cellar/openssl@1.1/1.1.1d/lib/libcrypto.dylib)
  else()
    include_directories(${OPENSSL_INCLUDE_DIR})
    list(APPEND LIB_LIST ${OPENSSL_LIBRARIES})
  endif()
endif()

file(GLOB sources "src/*.cc")

add_executable(${NAME} ${sources})
target_include_directories(${NAME} PUBLIC src ${OPENSSL_INCLUDE_DIR})
target_link_libraries(${NAME} ${OPENSSL_LIBRARIES})

# add_executable(Test src/test.cc src/circuit.cc src/prf.cc)
# target_include_directories(Test PUBLIC src ${OPENSSL_INCLUDE_DIR})
# target_link_libraries(Test ${OPENSSL_LIBRARIES})
