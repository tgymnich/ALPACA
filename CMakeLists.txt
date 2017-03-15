cmake_minimum_required( VERSION 3.10 FATAL_ERROR )

set( ALPACA_ENV "$ENV{ALPACA_ENVIRONMENT}" )
# Set the compilers to be used to build the targets.
# Note, this has to be done before defining the project.
INCLUDE("./cmake/compiler.cmake")

PROJECT( ALPACA )
# To compile ALPACA, the c++17 standard is required.
set( CMAKE_CXX_STANDARD 17 )

INCLUDE( FindMPI REQUIRED )

FIND_PACKAGE( HDF5 REQUIRED )
include_directories( ${HDF5_INCLUDE_DIRS} )

# Define the present git hash as compiler macro.
INCLUDE("./cmake/git_hash.cmake")
# Define warning flags used for compilation.
INCLUDE("./cmake/warning_flags.cmake")
# Define performance flags used for compilation (machine dependent).
INCLUDE("./cmake/performance_flags.cmake")
# Define an option to chosse the dimension of the build.
INCLUDE("./cmake/dimension.cmake")
# Define a target to create the doxygen documentation.
INCLUDE("./cmake/documentation.cmake")

# Include directories to know all necessary ALPACA headers.
INCLUDE_DIRECTORIES(src)
# Include directories of third party libraries
INCLUDE_DIRECTORIES(3rdParty/expression_toolkit)
INCLUDE_DIRECTORIES(3rdParty/tiny_xml)
INCLUDE_DIRECTORIES(3rdParty/FakeIt/single_header/catch)
INCLUDE_DIRECTORIES(3rdParty/Catch2/single_include/catch2)

file(GLOB_RECURSE SOURCE_FILES "src/*.cpp")
list(APPEND SOURCE_FILES "./3rdParty/tiny_xml/tinyxml2.cpp")

# Exclude user_expression.cpp to compile it as a separate library
list(FILTER SOURCE_FILES EXCLUDE REGEX ".*user_expression.cpp$")

include(CheckIPOSupported)
check_ipo_supported( RESULT IPOPOSSIBLE )

# Define the ALPACA executable.
add_executable(ALPACA ${SOURCE_FILES})
if( IPOPOSSIBLE )
   set_property(TARGET ALPACA PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif( IPOPOSSIBLE )
set_target_properties(ALPACA PROPERTIES COMPILE_FLAGS "${ALPACA_CXX_FLAGS}")
set( ALLUSEDFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_COMPILE_FLAGS} ${CMAKE_CXX_LINK_FLAGS} ${ALPACA_CXX_FLAGS}" )
if( IPOPOSSIBLE )
   set( ALLUSEDFLAGS "${ALLUSEDFLAGS} ${CMAKE_CXX_COMPILE_OPTIONS_IPO}" )
endif( IPOPOSSIBLE )
MESSAGE( STATUS "Alpaca compiles and links using the following flags (without warranty, check VERBOSE build)\n-- ${ALLUSEDFLAGS}" )
target_compile_definitions(ALPACA PUBLIC TEST_VIRTUAL=)

target_link_libraries( ALPACA ${HDF5_LIBRARIES} )

# Define the UserExpressions executable, encapsulating the expression toolkit from the rest of the program
find_library(UserExpr NAMES UserExpressions PATHS ${CMAKE_CURRENT_BINARY_DIR})
if(NOT UserExpr)
   MESSAGE( STATUS "UserExpressions not found - will be build" )
   add_library(UserExpressions STATIC src/user_expression.cpp )
   if( IPOPOSSIBLE )
      set_property(TARGET UserExpressions PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
   endif( IPOPOSSIBLE )
   target_link_libraries(ALPACA UserExpressions)
else(NOT UserExpr)
   MESSAGE( STATUS "UserExpressions found at ${UserExpr}" )
   target_link_libraries(ALPACA ${UserExpr})
endif(NOT UserExpr)

INCLUDE_DIRECTORIES(test)

# Define a target for unit tests.
file(GLOB_RECURSE TEST_FILES "test/*.cpp")
list(FILTER SOURCE_FILES EXCLUDE REGEX "src/main.cpp$")
list(APPEND TEST_FILES "${SOURCE_FILES}")

# Use additional flags for test coverage extraction
set(PACO_CXX_FLAGS "-O0 -fprofile-arcs -ftest-coverage")
set(GCOV_LIBRARY "-lgcov")

add_executable(Paco EXCLUDE_FROM_ALL ${TEST_FILES})
set_target_properties(Paco PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(Paco PROPERTIES COMPILE_FLAGS "${PACO_CXX_FLAGS}")
target_compile_definitions(Paco PUBLIC TEST_VIRTUAL=virtual)

target_link_libraries( Paco ${HDF5_LIBRARIES} )
target_link_libraries(Paco ${GCOV_LIBRARY})
target_link_libraries(Paco UserExpressions)
