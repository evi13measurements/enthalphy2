# LibsWESTA.cmake
# hjaaj Sep. 2017

add_library(
    westa
    ${WESTA_FIXED_FORTRAN_SOURCES}
    )

add_executable(
    westa.x
    ${CMAKE_SOURCE_DIR}/WESTA/westa/westa.F
    )

set_property(TARGET westa.x PROPERTY LINKER_LANGUAGE Fortran)

get_target_property(_incdirs dalton INCLUDE_DIRECTORIES)
# DALTON/sirius/ needed for mxpdim.h
set(_incdirs ${_incdirs} ${CMAKE_SOURCE_DIR}/DALTON/sirius)
set_target_properties(westa PROPERTIES INCLUDE_DIRECTORIES "${_incdirs}")

target_link_libraries(
    westa.x
    westa
    dalton
    ${DALTON_LIBS}
    ${EXTERNAL_LIBS}
    )

