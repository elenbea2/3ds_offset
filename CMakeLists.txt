cmake_minimum_required(VERSION 3.10)
project(Offset)

# Encontrar e configurar CGAL
find_package(CGAL REQUIRED)

# Adicionar o executável
add_executable(offset offset_mesh.cpp)

# Linkar com CGAL
target_link_libraries(offset CGAL::CGAL)

# Ativar C++17 (ou superior, se necessário)
set_target_properties(offset PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Incluir diretórios necessários
include_directories(${CGAL_INCLUDE_DIRS})

