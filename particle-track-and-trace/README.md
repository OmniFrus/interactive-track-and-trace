## Interactive track and trace
This folder contains the Vtk program which actually displays the simulated data. The code is driven by the `Program` class, which contains the upper level of the vtk pipeline: the class has attributes for a vtkRenderWindow and vtkRenderWindowInteractor. vtkRenderers are managed through an abstract `Layer` class, which the program keeps track of trough a vector attribute.

Each layer implementation contains and manages one vtkRenderer, this includes managing which layer of the vtkrenderwindow ths layer renders to. Currently implemented are three such layers:
 * the `BackgroundImage` class reads in image data and displays this to the screen on the 0th layer - the background.
 * the `EulerGlyphs` class renders a visualization of the Eulerian flow-velocities as a grid of arrow-glyphs (in which the direction and length of the glyph represents the direction and strength of the velocity at that point). Right now it spoofs the data for these glyphs, but this class will interface with the code for reading h5 data to accurately display the velocities at a given timestamp.
 * the `LagrangeGlyphs` class renders a given set of particles as circular glyphs. These particles are advected according to an advection function, which in this implementation is spoofed. Like the EglyphLayer class, this layer will interact with the code for advecting particles according to the actual dataset, to accurately simulate its particles.

The `LagrangeGlyphs` deserves some more explanation, as it depends on the `SpawnpointCallback` class to place particles in its dataset. The `SpawnpointCallback` makes use of the vtkCallbackCommand class and the vtk observer pattern to create new particles on mouseclick. It does so through a shared reference to the LagrangeGlyphs's `data` and `points` attributes, which the SpawnpointCallback then edits directlr. 

The program also adds a second observer to the vtk pattern through the `TimerCallbackCommand`. This class subscribes to a vtkTimerEvent to manage the simulation of the program. To this end the TimerCallbackCommand has attributes for a timestep (dt) and current time (time). On every controller, the current time is updated according to the dt attribute, and this change is propagated to the layers containing the data by use of the program and layer's `updateData()` functions.

## Location of data
The data path is hardcoded such that the following tree structure is assumed:
```
data/
  grid.h5
  hydrodynamic_U.h5
  hydrodynamic_V.h5
  map_2071-2067.png
  spawn_locations.csv
  
interactive-track-and-trace/
  particle-track-and-trace/
     ...
```
The csv file you find in the data folder are the spawn locations of all the litter particles for this simulation. They will always spawn in the same location to avoid different test results based on spawn location.
You can also spawn particles by left mouse clicking inside the simulation.
Particles are colored Red if they are considered beached and Blue if they are not.

## Compiling
You will need [vtk](https://gitlab.kitware.com/vtk/vtk), [netcdf-cxx](https://github.com/Unidata/netcdf-cxx4/),
[SFML](https://github.com/SFML/SFML).

I downloaded everything using vcpkg.
./vcpkg install vtk
./vcpkg install netcdf
./vcpkg install sfml
./vcpkg install hdf5:x64-windows
install PROJ.db -> https://download.osgeo.org/osgeo4w/

Let the current directory be the `src` directory. Run:
```shell
cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE="C:/Users/wesle/Documents/Universiteit/2024_2025/Thesis/Opdracht/vcpkg/scripts/buildsystems/vcpkg.cmake"
cmake --build build --config Release
```
Then run C:\Users\wesle\Documents\Universiteit\2024_2025\Thesis\Opdracht\interactive-track-and-trace\particle-track-and-trace\src\build\Release\ParticleTrackTrace.exe from the file explorer.


We can make changes to our simulation in the main.cpp file, where we can choose which boundary condition we want to run and which beaching condition we want to use.
Furthermore, we can run a python script which gives us all the statistics of what was ran before. This could be for a particle in particular or all at once.

