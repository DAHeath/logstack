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
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if (APPLE)
  set(OPENSSL_ROOT_DIR "/usr/local/opt/openssl")
endif()
find_package(OpenSSL REQUIRED)
find_package(EMP-TOOL REQUIRED)

file(GLOB sources "src/*.cc")

add_executable(${NAME} ${sources})
target_include_directories(${NAME} PUBLIC src ${OPENSSL_INCLUDE_DIR} /usr/local/include)
target_link_libraries(${NAME} ${OPENSSL_LIBRARIES})
