cmake_minimum_required(VERSION 3.2)

add_library(ini ini.c INIReader.cpp)
add_library(DB DB.c QV.c)
add_library(LA align.c)

add_library(PAF paf.c)
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( PAF ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

add_library(LAInterface LAInterface.cpp)
target_link_libraries(LAInterface LA DB PAF)


add_library(kmer_lookup kmer_lookup.c)
add_library(DW_banded DW_banded.c)
add_library(falcon falcon.c)
target_link_libraries(falcon kmer_lookup DW_banded)


