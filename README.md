# VItA

# Linux Installation

As a first step, we must clone the github repository in our local drive using the following commands:

    cd <my_vita_folder>
    mkdir vita_source
    mkdir vita_build
    git clone https://github.com/GonzaloMaso/VItA.git vita_source
    cd vita_build
    ccmake ../vita_source

Complete the CMAKE_INSTALL_PREFIX = <my_vita_folder>, press key "c" to configure and after "g" to generate make files. Once back in the terminal, we will build and install the library with the following commands:

    make
    make install

For parallel compilation with N threads, switch the command "make" by "make -jN", e.g., to parallelise with 8 threads execute "make -j8". Note that the "make" command will take more than 30 minutes as is also installing VTK 8.1 (a dependency of this library).



# Using VItA library



# About

The Virtual ITerative Angiogenesis library allows the generation of synthetic vasculatures mimicking the angiogenesis process. The models support the usage of in-vivo or experimental prior constraints over the geometrical description.
