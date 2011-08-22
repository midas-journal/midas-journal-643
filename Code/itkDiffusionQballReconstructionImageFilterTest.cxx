/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffusionTensor3DReconstructionImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2006-03-09 03:37:16 $
  Version:   $Revision: 1.10 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkDiffusionQballReconstructionImageFilter.h"

// itk includes
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaDataObject.h"


int main(int argc, char **argv)
{

  if( argc < 2 )
  {
    std::cout << "[FAILED]" << std::endl;
    std::cerr << "You must supply a filename of a diffusion weighted image" 
      << std::endl;
    return EXIT_FAILURE;
  }

  typedef short int          ReferencePixelType;
  typedef short int          GradientPixelType;
  typedef float              OdfPrecisionType;

  // The angular resolution of the resulting ODFs depends on the number of
  // directions distributed over the sphere
  const int nOdfDirections = 252;

  // Some typedefs
  typedef itk::DiffusionQballReconstructionImageFilter< 
      ReferencePixelType, GradientPixelType, OdfPrecisionType, nOdfDirections > 
                                          QBallReconstructionImageFilterType;
  typedef QBallReconstructionImageFilterType::GradientImagesType 
                                          GradientImageType;
  typedef QBallReconstructionImageFilterType::OutputImageType 
                                          OutputImageType;
  typedef QBallReconstructionImageFilterType::BZeroImageType 
                                          BZeroImageType;

  typedef itk::ImageFileReader< GradientImageType > 
                                          ReaderType;
  typedef itk::ImageFileWriter< BZeroImageType >    
                                          BZeroWriterType;
  typedef itk::VectorImage<OdfPrecisionType, 3>     
                                          VarVecImgType;
  typedef itk::ImageFileWriter< VarVecImgType >     
                                          OdfWriterType;

  typedef vnl_vector_fixed< double, 3 >   GradientDirectionType;
  typedef itk::VectorContainer< unsigned int, 
    GradientDirectionType >               GradientDirectionContainerType;

  // Read the input file
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);

  try
  {
    reader->UpdateLargestPossibleRegion();
  }
  catch (itk::ExceptionObject& e)
  {
    std::cout << "[FAILED]" << std::endl;
    std::cerr << "Exception detected while reading " << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
  }

  GradientImageType::Pointer gradImg = reader->GetOutput();

  // read gradient information from meta data dictionary
  itk::MetaDataDictionary imgMetaDictionary = gradImg->GetMetaDataDictionary();    
  std::vector<std::string> imgMetaKeys = imgMetaDictionary.GetKeys();
  std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
  std::string metaString;

  GradientDirectionType vect3d;
  GradientDirectionContainerType::Pointer gradientDirections = 
    GradientDirectionContainerType::New();
  double bValue;

  int numberOfImages = 0;
  int numberOfGradientImages = 0;
  bool readb0 = false;

  for (; itKey != imgMetaKeys.end(); itKey ++)
  {
    double x,y,z;

    itk::ExposeMetaData<std::string> (imgMetaDictionary, *itKey, metaString);
    if (itKey->find("DWMRI_gradient") != std::string::npos)
    { 
      std::cout << *itKey << " ---> " << metaString << std::endl;      
      sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
      vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;
      gradientDirections->InsertElement( numberOfImages, vect3d );
      ++numberOfImages;
      // If the direction is 0.0, this is a reference image
      if (vect3d[0] == 0.0 &&
        vect3d[1] == 0.0 &&
        vect3d[2] == 0.0)
      {
        continue;
      }
      ++numberOfGradientImages;;
    }
    else if (itKey->find("DWMRI_b-value") != std::string::npos)
    {
      std::cout << *itKey << " ---> " << metaString << std::endl;      
      readb0 = true;
      bValue = atof(metaString.c_str());
    }
  }

  std::cout << "Number of gradient images: "
    << numberOfGradientImages
    << " and Number of reference images: "
    << numberOfImages - numberOfGradientImages
    << std::endl;

  if(!readb0)
  {
    std::cerr << "BValue not specified in header file" << std::endl;
  }
  else
  {
    std::cout << "B-Value: " << bValue << std::endl;
  }

  // init and run the reconstruction
  QBallReconstructionImageFilterType::Pointer qBallReconstructionFilter = 
    QBallReconstructionImageFilterType::New();
  qBallReconstructionFilter->SetThreshold(70.0);
  qBallReconstructionFilter->SetNormalizationMethod(
    QBallReconstructionImageFilterType::QBR_STANDARD);
  qBallReconstructionFilter->SetBValue(bValue);
  qBallReconstructionFilter->SetGradientImage( 
    gradientDirections, reader->GetOutput() );
  std::cout << std::endl << "This filter is using " << 
    qBallReconstructionFilter->GetNumberOfThreads() << " threads " << std::endl;
  qBallReconstructionFilter->Update();

  // convert result to vector image that can be written to file
  // (simple image with pixeltype vector not supported by writer)
  VarVecImgType::Pointer vecImg = VarVecImgType::New();
  vecImg->SetSpacing( qBallReconstructionFilter->GetOutput()->GetSpacing() );
  vecImg->SetOrigin( qBallReconstructionFilter->GetOutput()->GetOrigin() );
  vecImg->SetDirection( qBallReconstructionFilter->GetOutput()->GetDirection());
  vecImg->SetLargestPossibleRegion( 
    qBallReconstructionFilter->GetOutput()->GetLargestPossibleRegion());
  vecImg->SetBufferedRegion( 
    qBallReconstructionFilter->GetOutput()->GetLargestPossibleRegion() );
  vecImg->SetVectorLength(nOdfDirections);
  vecImg->Allocate();
  
  itk::ImageRegionIterator<VarVecImgType> ot (
    vecImg, vecImg->GetLargestPossibleRegion() );
  ot = ot.Begin();
  
  itk::ImageRegionIterator<OutputImageType> it (
    qBallReconstructionFilter->GetOutput(), 
    qBallReconstructionFilter->GetOutput()->GetLargestPossibleRegion() );
  it = it.Begin();
  
  for (it = it.Begin(); !it.IsAtEnd(); ++it)
  {
    itk::Vector<OdfPrecisionType,nOdfDirections> vec = it.Get();
    VarVecImgType::PixelType varvec(vec.GetDataPointer(), nOdfDirections);
    ot.Set(varvec);
    ++ot;
  }

  // init file writer and write output to disk
  OdfWriterType::Pointer odfWriter = OdfWriterType::New();
  odfWriter->SetInput(vecImg);
  ::itk::OStringStream baseName;
  baseName << argv[1] << ".odfs.nrrd";

  try
  {
    odfWriter->SetFileName( baseName.str().c_str() );
    odfWriter->Update();
  }
  catch (itk::ExceptionObject& e)
  {
    std::cout << "[FAILED]" << std::endl;
    std::cerr << "Error during write of " << baseName.str();
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
  }

  // init second file writer and write extracted b-zero image to disk
  BZeroWriterType::Pointer b0Writer = BZeroWriterType::New();
  b0Writer->SetInput(qBallReconstructionFilter->GetBZeroImage());
  ::itk::OStringStream baseNameB0;
  baseNameB0 << argv[1] << ".b0.nrrd";

  try
  {
    b0Writer->SetFileName( baseNameB0.str().c_str() );
    b0Writer->Update();
  }
  catch (itk::ExceptionObject& e)
  {
    std::cout << "[FAILED]" << std::endl;
    std::cerr << "Error during write of " << baseNameB0.str();
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
  }

  std::cout << "[PASSED]" << std::endl;    
  return EXIT_SUCCESS;
}

