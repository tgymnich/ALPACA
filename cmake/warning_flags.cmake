MESSAGE( STATUS "Set warning flags" )

set(STANDARD_WARNINGS "-Wall -pedantic -W -Wformat -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wreturn-type -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STANDARD_WARNINGS}")
