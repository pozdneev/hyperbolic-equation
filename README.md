# README

The simple code to solve hyperbolic equation in 2D and 3D for Dirichlet boundary conditions (designated as `1`) and periodic boundary conditions (designated as `p`).

Type `make` in each subdirectory to obtain the binaries.

Modify `Makefile` if you do not want to include Electric Fence memory debugger.

You need to have `libsilo` installed in your system. Specify the path to it by setting the `SILO_DIR` variable in a makefile.

The binaries produce a sequence of silo-files. You can visualize them with the [VisIt](http://visit.llnl.gov) software. Modify the paths in `visit.session*` files accordingly. Some video files obtained from these calculations are available on [YouTube](https://www.youtube.com/playlist?list=PLDudW2rTDAfoB61KTFHjlcTeHj7-3R0gR).

The code was prepared in 2011 for the students studying HPC at the [Faculty of Computational Mathematics and Cybernetics](http://cmc.msu.ru) at [Lomonosov Moscow State University](http://msu.ru).
