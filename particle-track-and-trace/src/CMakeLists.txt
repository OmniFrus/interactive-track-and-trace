cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# Always use vcpkg toolchain
set(CMAKE_TOOLCHAIN_FILE "C:/Users/wesle/Documents/Universiteit/2024_2025/Thesis/Opdracht/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")

project(ParticleTrackTrace)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_definitions(_USE_MATH_DEFINES)

# --- VTK via find_package
find_package(VTK REQUIRED COMPONENTS
  GeovisCore
  CommonColor
  CommonCore
  CommonDataModel
  FiltersGeneral
  FiltersGeometry
  FiltersTexture
  FiltersProgrammable
  FiltersSources
  ImagingSources
  ImagingGeneral
  ImagingCore
  InteractionStyle
  IOImage
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

# --- SFML via find_package
find_package(SFML COMPONENTS Network Graphics Window Audio System CONFIG REQUIRED)
find_package(netCDF REQUIRED)

# --- Manual NetCDF setup (no find_package!)
set(NETCDF_CXX_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../../../vcpkg/installed/x64-windows/include")
set(NETCDF_CXX_LIB_DIR "${CMAKE_SOURCE_DIR}/../../../vcpkg/installed/x64-windows/lib")

# --- Executable
add_executable(ParticleTrackTrace MACOSX_BUNDLE
    main.cpp
    layers/BackgroundImage.cpp
    layers/BackgroundImage.h
    layers/EulerGlyphs.cpp
    layers/EulerGlyphs.h
    layers/LagrangeGlyphs.cpp
    layers/LagrangeGlyphs.h
    layers/Layer.cpp
    layers/Layer.h
    Program.cpp
    Program.h
    commands/TimerCallbackCommand.h
    commands/TimerCallbackCommand.cpp
    commands/SpawnPointCallback.h
    commands/SpawnPointCallback.cpp
    commands/CameraMoveCallback.cpp
    commands/CameraMoveCallback.h
    CartographicTransformation.cpp
    advection/kernel/AdvectionKernel.h
    advection/kernel/EulerAdvectionKernel.cpp
    advection/kernel/EulerAdvectionKernel.h
    advection/interpolate.cpp
    advection/interpolate.h
    advection/readdata.cpp
    advection/readdata.h
    advection/kernel/RK4AdvectionKernel.cpp
    advection/kernel/RK4AdvectionKernel.h
    advection/UVGrid.cpp
    advection/UVGrid.h
    advection/Vel.cpp
    advection/Vel.h
    advection/kernel/SnapBoundaryConditionKernel.h
    advection/kernel/SnapBoundaryConditionKernel.cpp
    layers/Character.h
    layers/Character.cpp
    commands/CharacterMoveCallback.h
    commands/CharacterMoveCallback.cpp
    layers/ParticleCollision.h
    layers/ParticleCollision.cpp
    collisions/ParticleCollisionCallback.h
    collisions/DummyCollisionHandler.h
    collisions/DummyCollisionHandler.cpp
    layers/Vessels.h
    layers/Vessels.cpp
    vesselroutes/VesselRoute.h
    vesselroutes/QuadraticBezierRoute.cpp
    vesselroutes/QuadraticBezierRoute.h
    vesselroutes/VesselRouteFactory.cpp
    vesselroutes/VesselRouteFactory.h
    layers/FoodSpawner.cpp
    layers/FoodSpawner.h
    collisions/ParticleRemover.cpp
    collisions/ParticleRemover.h
    vesselroutes/StraightLineRoute.cpp
    vesselroutes/StraightLineRoute.h
    layers/Health.cpp
    layers/Health.h
    collisions/FoodPickup.cpp
    collisions/FoodPickup.h
    collisions/DebrisPickup.cpp
    collisions/DebrisPickup.h
    layers/Camera.cpp
    layers/Camera.h
    layers/GameOverScreen.cpp
    layers/GameOverScreen.h
    gameovers/GameoverCallback.h
    layers/Timer.cpp
    layers/Timer.h
    commands/DismissGameoverScreen.cpp
    commands/DismissGameoverScreen.h
    layers/DayCounter.cpp
    layers/DayCounter.h
    layers/Badges.cpp
    layers/Badges.h
    badges/Badge.h
    badges/Badge.cpp
    badges/StaticBadge.cpp
    badges/StaticBadge.h
    statistics/Statistic.h
    collisions/TrackedCollision.cpp
    collisions/TrackedCollision.h
    badges/Achievement.h
    badges/SimpleAchievement.cpp
    badges/SimpleAchievement.h
    badges/FileBadgeFactory.cpp
    badges/FileBadgeFactory.h
    layers/StatisticsManager.cpp
    layers/StatisticsManager.h
    statistics/Statistic.cpp
    layers/BadgeCounter.cpp
    layers/BadgeCounter.h
    statistics/BadgesAcquired.cpp
    statistics/BadgesAcquired.h
    badges/TrackedAchievement.cpp
    badges/TrackedAchievement.h
    commands/CharacterMouseMoveCallback.cpp
    commands/CharacterMouseMoveCallback.h
    layers/CharacterNoStop.h
    layers/CharacterNoStop.cpp
    collisions/SoundEffect.cpp
    collisions/SoundEffect.h
    layers/SplashScreen.cpp
    layers/SplashScreen.h
    commands/DisplaySplashScreenCommand.cpp
    commands/DisplaySplashScreenCommand.h
)

# --- Include + link NetCDF manually
target_include_directories(ParticleTrackTrace PRIVATE ${NETCDF_CXX_INCLUDE_DIR})
target_link_directories(ParticleTrackTrace PRIVATE ${NETCDF_CXX_LIB_DIR})

target_link_libraries(ParticleTrackTrace PRIVATE
    netcdf
    netcdf-cxx4
    ${VTK_LIBRARIES}
    sfml-audio
)

# --- Required by VTK
vtk_module_autoinit(
  TARGETS ParticleTrackTrace
  MODULES ${VTK_LIBRARIES}
)

if(MSVC)
  add_compile_options(/wd4267 /wd4244)  # Optional: suppress noisy warnings for now
endif()