add_executable(create-htgfile create-htgfile.cpp)
add_executable(poisson poisson.cpp)
add_executable(tgv tgv.cpp)

target_link_libraries(create-htgfile PRIVATE oktal)
target_link_libraries(poisson PRIVATE oktal)
target_link_libraries(tgv PRIVATE oktal)