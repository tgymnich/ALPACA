if( ALPACA_ENV STREQUAL "LRZ" )
    MESSAGE( STATUS "Compile on LRZ" )
    set( CMAKE_CXX_COMPILER "mpiCC" )
elseif( ALPACA_ENV STREQUAL "AER" )
    MESSAGE( STATUS "Compile on AER machine" )
    set(CMAKE_CXX_COMPILER "/global/mpich-3.1/bin/mpic++" )
else( ALPACA_ENV STREQUAL "LRZ" )
    MESSAGE( STATUS "Compile with user default compilers" )
    set( CMAKE_CXX_COMPILER "mpicxx" )
endif( ALPACA_ENV STREQUAL "LRZ" )

# Use gold linker, if it is installed and not disabled
option( USEGOLD "Gold Linker" ON )
if(UNIX AND NOT APPLE AND USEGOLD)
   execute_process(COMMAND ${CMAKE_CXX_COMPILER} -fuse-ld=gold -Wl,--version ERROR_QUIET OUTPUT_VARIABLE ld_version)
   if("${ld_version}" MATCHES "GNU gold")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold -Wl,--disable-new-dtags ")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=gold -Wl,--disable-new-dtags ")
      MESSAGE( STATUS "Using gold Linker" )
   endif("${ld_version}" MATCHES "GNU gold")
endif(UNIX AND NOT APPLE AND USEGOLD )