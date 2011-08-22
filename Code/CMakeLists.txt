INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

PROJECT(QBallReconstruction)
ENABLE_TESTING()
INCLUDE(Dart)

SET(Required_Packages
ITK
)

FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

SET(Libraries
ITKCommon
ITKIO
)

SET(CurrentExe "ImageCompare")
ADD_EXECUTABLE(${CurrentExe} ImageCompare.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
SET(CurrentExe "QBallReconstructionTest")
ADD_EXECUTABLE(${CurrentExe} itkDiffusionQballReconstructionImageFilterTest.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
SET(CurrentExe "QBallReconstructionExample")
ADD_EXECUTABLE(${CurrentExe} DiffusionQballReconstructionImageFilter.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
                              
ADD_TEST(ReconstructImage QBallReconstructionTest ${CMAKE_SOURCE_DIR}/dwi.nrrd)
ADD_TEST(CompareImage1 ImageCompare ${CMAKE_SOURCE_DIR}/odfs.nrrd ${CMAKE_SOURCE_DIR}/dwi.nrrd.odfs.nrrd)
ADD_TEST(CompareImage2 ImageCompare ${CMAKE_SOURCE_DIR}/b0.nrrd ${CMAKE_SOURCE_DIR}/dwi.nrrd.b0.nrrd)

