# Multi-layer solver for XPBD
Matlab code for the multi-layer solver project. This code base provides an implementation of the paper by building onto the adaptive rigidification codebase of Siggraph 2022. I note, that this is simply a quick snapshot of the code to provide a prototype. The code itself requires heavy modifications for production use. Likewise, I turned some functions in the xpbdLayer classes to imperative so I could use the matlab gpu code generation. The generated gpu code was not great, and the changes might have broken some of parts of the code. I do not plan to maintain this codebase as I plan to move to a faster open source language for future work. However, this should give a good idea of how to get an implementation going in your own language.

Matlab is intuitive which makes the project easy to setup and results quickly reproducable. To test the various scenes featured in the paper, simply follow the instructions below then open an example file and click play in the editor. Alternatively the name of the example can be called in the command prompt which will launch the simulation. The scenes feature heterogenous materials, multiple objects, scripted animations, surface mapping and more!

## How to use

### Dependencies
required: opengl, visual studio (or any C++ compilers), cmake, matlab, gptoolbox(for 3D)
GPtoolbox might require: 
	embree https://github.com/embree/embree
	libigl https://github.com/libigl/libigl
	boost 1.48 https://www.boost.org/users/history/version_1_48_0.html
	to be detected automatically on windows use the path C:\local\boost_1_48_0\boost\*
Tip: install visual studio before cmake so its compiler is automatically added to cmake

GpToolBox is already included in the repo as a subrepo. Simply call `git submodule update --init --recursive`.
cmake is needed to build gptoolbox's mex. You will need the boost library. 
I recommend not to build gptoolbox with eltopo if you are on windows.

We depend on GPToolbox for contacts(signed distances) as well as loading some mesh files.

Some matlab addons are needed or recommended such as 
Aerospace toolbox, 
Computer Vision Toolbox, 
Deep Learning Toolbox, 
Image processing Toolbox
Matlab Compiler, 
Parallel Computing Toolbox,
and Signal Processing Toolbox. Those can be found and downloaded from matlab's add-on menu.

### Installation
#### run the following commands in a command prompt
```
git submodule update --init --recursive ;
#if the recursive update fails then manually go to lib/gptoolbox/mex and build
cd lib/gptoolbox/mex ;
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DWITH_ELTOPO=false ;
cmake --build build --target install --config Release ;
xcopy /F /Y Release ./ ;
```

#### run the following commands in the matlab console
Once the project is set up, add the files to path by calling `addAdaptiveRigidificationToPath();` in the command prompt.
Followed by `mexAdaptiveRigidification();` to build the project's mex code.
```
addAdaptiveRigidificationToPath();
savepath();
mexAdaptiveRigidification();
```
By now you should be all set.
Try running examples in 2d/examples or 3d/examples.

### FAQ
If you get a signed-distance error when running the 3D code, then you probably did not set up the compiled mex from gptoolbox properly.
