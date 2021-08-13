# VItA

# Linux Installation

As a first step, we must clone the github repository in our local drive using the following commands:

    cd <my_vita_folder>
    mkdir vita_source
    mkdir vita_build
    git clone https://github.com/GonzaloMaso/VItA.git vita_source
    cd vita_build
    ccmake ../vita_source

In the ccmake interface, press key "c" to configure the cmake project, complete the CMAKE_INSTALL_PREFIX = <my_vita_folder>, press key "c" again, and then press "g" to generate make files. Once back in the terminal, we will build and install the library with the following commands:

    make
    make install

For parallel compilation with N threads, switch the command "make" by "make -jN", e.g., to parallelise with 8 threads execute "make -j8". Note that the "make" command may take more than 30 minutes as is also installing VTK 8.1 (a dependency of this library).

# Using VItA library

Once you create your project, cpp files using VItA should be compiled and linked as follows

    g++ <cpp_filename>.cpp -Wall -std=c++11 -O3 -I<vita_folder>/vita_build/include/vtk-8.1 -I<vita_folder>/include/vita_source -L<vita_folder>/vita_build/lib -L<vita_folder>/lib -o <executable_filename> -lVItA -lvtkCommonCore-8.1 -lvtkCommonDataModel-8.1 -lvtkCommonExecutionModel-8.1 -lvtkFiltersModeling-8.1 -lvtkIOCore-8.1 -lvtkIOLegacy-8.1 -lvtkIOXML-8.1 -lvtkIOGeometry-8.1 -lvtkInfovisCore-8.1 -lvtkFiltersGeneral-8.1 -lvtkFiltersCore-8.1 -lvtkCommonTransforms-8.1 -lvtkIOXMLParser-8.1

An example can be found in the released version 0.2 (see tag VItA v0.2 - https://github.com/GonzaloMaso/VItA/releases/tag/v0.2).

# About

The Virtual ITerative Angiogenesis library allows the generation of synthetic vasculatures mimicking the angiogenesis process. The models support the usage of in-vivo or experimental prior constraints over the geometrical description.

# Publications

Maso Talou, G. D., et al. "Adaptive constrained constructive optimisation for complex vascularisation processes." Scientific Reports 11.1 (2021): 1-22. - https://www.nature.com/articles/s41598-021-85434-9
