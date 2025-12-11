add_executable(create-htgfile create-htgfile.cpp)
add_executable(poisson poisson.cpp)

target_link_libraries(create-htgfile PRIVATE oktal)
target_link_libraries(poisson PRIVATE oktal)