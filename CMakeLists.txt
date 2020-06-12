cmake_minimum_required(VERSION 3.10)

project(VItA)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_CXX_EXTENSIONS false)

add_compile_options(-Wall)

# Setting components of the library

set(CORE_CPP 
    core/FixedPerfusionRadiusTreeGenerator.cpp 
    core/FRRSTreeGenerator.cpp 
    core/GeneratorData.cpp 
    core/GeneratorDataMonitor.cpp 
    core/StagedFixedPerfusionRadiusTreeGenerator.cpp
    core/StagedFRROTreeGenerator.cpp 
    core/StagedFRRSTreeGenerator.cpp)

set(CORE_HPP
    core/FixedPerfusionRadiusTreeGenerator.h
    core/FRRSTreeGenerator.h
    core/GeneratorData.h
    core/GeneratorDataMonitor.h
    core/StagedFixedPerfusionRadiusTreeGenerator.h
    core/StagedFRROTreeGenerator.h
    core/StagedFRRSTreeGenerator.h
    )

set(CREATORS_CPP
    creators/AbstractCreator.cpp
    creators/CylinderCreator.cpp
    creators/ParallelepipedCreator.cpp
    creators/SphereCreator.cpp)

set(CREATORS_HPP
    creators/AbstractCreator.h
    creators/CylinderCreator.h
    creators/ParallelepipedCreator.h
    creators/SphereCreator.h)

set(FILTERS_CPP
    filters/AbstractVesselFilter.cpp
    filters/VesselFilterByBranchingMode.cpp
    filters/VesselFilterByStage.cpp
    filters/VesselFilterComposite.cpp)

set(FILTERS_HPP
    filters/AbstractVesselFilter.h
    filters/VesselFilterByBranchingMode.h
    filters/VesselFilterByStage.h
    filters/VesselFilterComposite.h)

set(IO_CPP
    io/CSVWriter.cpp
    io/EXFTreeWriter.cpp
    io/StagedFRROTreeGeneratorLogger.cpp
    io/VTKConverter.cpp
    io/VTKObjectTreeElementalWriter.cpp
    io/VTKObjectTreeNodalWriter.cpp
    io/VTKObjectTreeSplinesNodalWriter.cpp
    io/VTKObjectTreeWriter.cpp)

set(IO_HPP
    io/CSVWriter.h
    io/EXFTreeWriter.h
    io/StagedFRROTreeGeneratorLogger.h
    io/VTKConverter.h
    io/VTKObjectTreeElementalWriter.h
    io/VTKObjectTreeNodalWriter.h
    io/VTKObjectTreeSplinesNodalWriter.h
    io/VTKObjectTreeWriter.h)

set(IO_TASK_CPP
    io/task/AbstractSavingTask.cpp
    io/task/CheckpointSavingTask.cpp
    io/task/VisualizationSavingTask.cpp
    )

set(IO_TASK_HPP
    io/task/AbstractSavingTask.h
    io/task/CheckpointSavingTask.h
    io/task/VisualizationSavingTask.h)

set(CONSTRAINTS_HPP
    constrains/AbstractConstraintFunction.h
    constrains/bloodViscosity.h
    constrains/ConstantConstraintFunction.h
    constrains/ConstantPiecewiseConstraintFunction.h
    constrains/potentialLaw.h
    constrains/symmetryLaw.h
    )

set(STATS_CPP
    stats/ObjectTreeIndexCreator.cpp
    stats/ObjectTreeStatsManager.cpp
    stats/StructuredTreeIndexCreator.cpp
    stats/StructuredTreeStatsManager.cpp
    stats/VesselObjectHandler.cpp
    stats/VesselStructHandler.cpp)

set(STATS_HPP
    stats/ObjectTreeIndexCreator.h
    stats/ObjectTreeStatsManager.h
    stats/StructuredTreeIndexCreator.h
    stats/StructuredTreeStatsManager.h
    stats/VesselObjectHandler.h
    stats/VesselStructHandler.h)

set(STATS_MANIPULATORS_CPP
    stats/manipulators/Abstract0DStatManipulator.cpp
    stats/manipulators/Abstract1DStatManipulator.cpp
    stats/manipulators/AbstractStructuredStatManipulator.cpp
    stats/manipulators/MeanStatManipulator.cpp
    stats/manipulators/MeanStructStatManipulator.cpp
    stats/manipulators/PercentileStatManipulator.cpp
    stats/manipulators/PercentileStructStatManipulator.cpp
    stats/manipulators/SeriesStatManipulator.cpp
    stats/manipulators/StdStatManipulator.cpp
    stats/manipulators/StdStructStatManipulator.cpp)

set(STATS_MANIPULATORS_HPP
    stats/manipulators/Abstract0DStatManipulator.h
    stats/manipulators/Abstract1DStatManipulator.h
    stats/manipulators/AbstractStructuredStatManipulator.h
    stats/manipulators/MeanStatManipulator.h
    stats/manipulators/MeanStructStatManipulator.h
    stats/manipulators/PercentileStatManipulator.h
    stats/manipulators/PercentileStructStatManipulator.h
    stats/manipulators/SeriesStatManipulator.h
    stats/manipulators/StdStatManipulator.h
    stats/manipulators/StdStructStatManipulator.h)

set(STRUCTURES_HPP
    structures/CCOCommonStructures.h)

set(STRUCTURES_DOMAIN_CPP
    structures/domain/AbstractDomain.cpp
    structures/domain/CompositeDistributionGenerator.cpp
    structures/domain/DistributionGenerator.cpp
    structures/domain/DomainNVR.cpp
    structures/domain/DummyDomain.cpp
    structures/domain/IDomainObservable.cpp
    structures/domain/IDomainObserver.cpp
    structures/domain/IntersectionVascularizedDomain.cpp
    structures/domain/NormalDistributionGenerator.cpp
    structures/domain/PartiallyVascularizedDomain.cpp
    structures/domain/SimpleDomain2D.cpp
    structures/domain/SimpleDomain.cpp
    structures/domain/StagedDomain.cpp
    structures/domain/TreeProjector.cpp
    structures/domain/UniformDistributionGenerator.cpp)

set(STRUCTURES_DOMAIN_HPP
    structures/domain/AbstractDomain.h
    structures/domain/CompositeDistributionGenerator.h
    structures/domain/DistributionGenerator.h
    structures/domain/DomainNVR.h
    structures/domain/DummyDomain.h
    structures/domain/IDomainObservable.h
    structures/domain/IDomainObserver.h
    structures/domain/IntersectionVascularizedDomain.h
    structures/domain/NormalDistributionGenerator.h
    structures/domain/PartiallyVascularizedDomain.h
    structures/domain/SimpleDomain2D.h
    structures/domain/SimpleDomain.h
    structures/domain/StagedDomain.h
    structures/domain/TreeProjector.h
    structures/domain/UniformDistributionGenerator.h)

set(STRUCTURES_TREE_CPP
    structures/tree/AbstractCostEstimator.cpp
    structures/tree/AbstractObjectCCOTree.cpp
    structures/tree/AbstractStructuredCCOTree.cpp
    structures/tree/AdimSproutingVolumetricCostEstimator.cpp
    structures/tree/FRRCCOSTree.cpp
    structures/tree/FRRVariableViscosityCCOSTree.cpp
    structures/tree/FRRVaViOptCCOSTree.cpp
    structures/tree/SingleVesselCCOOTree.cpp
    structures/tree/SproutingVolumetricCostEstimator.cpp
    structures/tree/VolumetricCostEstimator.cpp)

set(STRUCTURES_TREE_HPP
    structures/tree/AbstractCostEstimator.h
    structures/tree/AbstractObjectCCOTree.h
    structures/tree/AbstractStructuredCCOTree.h
    structures/tree/AdimSproutingVolumetricCostEstimator.h
    structures/tree/FRRCCOSTree.h
    structures/tree/FRRVariableViscosityCCOSTree.h
    structures/tree/FRRVaViOptCCOSTree.h
    structures/tree/SingleVesselCCOOTree.h
    structures/tree/SproutingVolumetricCostEstimator.h
    structures/tree/VolumetricCostEstimator.h)

set(STRUCTURES_VASCULAR_ELEMENTS_CPP
    structures/vascularElements/AbstractVascularElement.cpp
    structures/vascularElements/MultiSegmentVessel.cpp
    structures/vascularElements/SingleVessel.cpp)

set(STRUCTURES_VASCULAR_ELEMENTS_HPP
    structures/vascularElements/AbstractVascularElement.h
    structures/vascularElements/MultiSegmentVessel.h
    structures/vascularElements/SingleVessel.h)

set(UTILS_CPP
    utils/MemoryMonitor.cpp)

set(UTILS_HPP
    utils/MemoryMonitor.h)

set(SOURCE_CPP ${CORE_CPP} ${CREATORS_CPP} ${FILTERS_CPP} ${IO_CPP}  ${IO_TASK_CPP} ${STATS_CPP}
    ${STATS_MANIPULATORS_CPP} ${STRUCTURES_DOMAIN_CPP} ${STRUCTURES_TREE_CPP}
    ${STRUCTURES_VASCULAR_ELEMENTS_CPP} ${UTILS_CPP})

set(SOURCE_HPP ${CORE_HPP} ${CREATORS_HPP} ${FILTERS_HPP} ${IO_HPP}  ${IO_TASK_HPP} ${STATS_HPP}
    ${STATS_MANIPULATORS_HPP} ${STRUCTURES_DOMAIN_HPP} ${STRUCTURES_TREE_HPP}
    ${STRUCTURES_VASCULAR_ELEMENTS_HPP} ${UTILS_HPP} ${CONSTRAINTS_HPP} ${STRUCTURES_HPP})

option(BUILD_SHARED_LIBS "Build VItA as a shared library." ON)

option(OPEN_MP "Uses OpenMP Parallelization." ON)
if(OPEN_MP)
    add_compile_options(-fopenmp)
endif()

option(DEBUG "Allows debugging." ON)
if(DEBUG)
    add_compile_options(-g)
endif()

option(OPTM "Allow optimization." ON)
if(OPTM)
    add_compile_options(-O3)
endif()

add_library(vita ${SOURCE_CPP})

#[[
set(VTK_LIBS 
    vtkCommonCore-9.0
    vtkCommonDataModel-9.0
    vtkCommonExecutionModel-9.0
    vtkFiltersFlowPaths-9.0
    vtkFiltersModeling-9.0
    vtkIOCore-9.0
    vtkIOLegacy-9.0
    vtkIOXML-9.0
    vtkIOGeometry-9.0)

set(VTK_INCLUDE /usr/local/include/vtk-9.0)
]]

set(VTK_LIBS 
    vtkCommonCore-8.1
    vtkCommonDataModel-8.1
    vtkCommonExecutionModel-8.1
    vtkFiltersFlowPaths-8.1
    vtkFiltersModeling-8.1
    vtkIOCore-8.1
    vtkIOLegacy-8.1
    vtkIOXML-8.1
    vtkIOGeometry-8.1)

set(VTK_INCLUDE /usr/local/include/vtk-8.1)

target_include_directories(vita PUBLIC ${VTK_INCLUDE})
target_link_libraries(vita PUBLIC ${VTK_LIBS} gomp)

# set_target_properties(vita PROPERTIES PUBLIC_HEADER "${SOURCE_HPP}")

install(TARGETS vita
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install(DIRECTORY "${CMAKE_SOURCE_DIR}/"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include/vita"
        FILES_MATCHING PATTERN "*.h"
        REGEX "doc/*" EXCLUDE
        REGEX "build_/*" EXCLUDE
        REGEX ".git/*" EXCLUDE)
