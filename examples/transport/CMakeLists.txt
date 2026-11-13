cmake_minimum_required(VERSION 3.20)
project(CCGeant4Transport LANGUAGES CXX)
find_package(Geant4 11.4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(cc_transport main.cc
  src/DetectorConstruction.cc
  src/PhysicsList.cc
  src/PrimaryGeneratorAction.cc
  src/ActionInitialization.cc
  src/RunAction.cc
  src/EventAction.cc
  src/SteppingAction.cc)
target_include_directories(cc_transport PRIVATE include)
target_compile_features(cc_transport PRIVATE cxx_std_17)
target_link_libraries(cc_transport PRIVATE ${Geant4_LIBRARIES})
file(COPY macros DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
