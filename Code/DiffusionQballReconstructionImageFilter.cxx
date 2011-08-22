/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: DiffusionTensor3DReconstructionImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2007-12-08 18:26:26 $
  Version:   $Revision: 1.0 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
// 
// This example shows how to use the DiffusionQballReconstructionImageFilter
// to reconstruct an image of Orientation Distribution Functions (ODFs)
// from Diffusion weighted images. See the documentation of 
// DiffusionQballReconstructionImageFilter first.
//
// The example takes diffusion weighted images in the Nrrd format, writes out
// the reconstructed ODF image and the mean of all reference images.
// 
// Acquiring sample datasets:
//  1. Get the DWI datasets from 
//        ftp://public.kitware.com/pub/namic/DTI/Data/dwi.nhdr
//        ftp://public.kitware.com/pub/namic/DTI/Data/dwi.img.gz (gunzip this)
//     These datasets contain a reference T1 image and 30 diffusion weighted
//     images. See the nrrd header for details such as B value etc..
//
//  2. Run the example with the following args
//       dwi.nhdr 80 Tensors.mhd FractionalAnisotropy.mhd RelativeAnisotropy.mhd 1
//  
//  3. You should find 30 gradient images, 1 reference image, the FA and RA images
//     in your working directory, which you can fire up in your favourite volume
//     browser.
//
// This work is part of the National Alliance for Medical Image
// Computing (NAMIC), funded by the National Institutes of Health
// through the NIH Roadmap for Medical Research, Grant U54 EB005149.
//
// Additional documentation: For details on the Nrrd format for DTI, see
// http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
//

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// remove the following line to disable timing
#define __QBR_TIMING__

#include "itkDiffusionQballReconstructionImageFilter.h"

// itk includes
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaDataObject.h"

int main(int argc, char **argv)
{

  if(argc < 3)
  {
    std::cerr << "Reconstruct ODFs and reference b-zero image from NRRD file "
      << "and write them as images" << std::endl << std::endl;
      
    std::cerr << "Usage: " << argv[0] << " NrrdFileName(.nhdr/.nrrd) threshold(on B0)"
      << " [ODF-ImageFileName] [b-zero-ImageFileName] [#RepetitionsForTimingAnalysis] [#Threads]" << std::endl << std::endl;
      
    std::cerr << "\tExample args: dwi1.nrrd 70 odfs.nrrd b0.nrrd 10 2" << std::endl;
    std::cerr << "\tExample args: dwi2.nrrd 80" << std::endl;
    return EXIT_FAILURE;
  }

  typedef short int          ReferencePixelType;
  typedef short int          GradientPixelType;
  typedef float              OdfPrecisionType;

  // The angular resolution of the resulting ODFs depends on the number of
  // directions distributed over the sphere
  const int nOdfDirections = 252;

  // Here we instantiate the DiffusionQballReconstructionImageFilter class.
  // The class is templated over the pixel types of the reference, gradient
  // and the to be created ODF pixel's precision (We use float here). Further
  // template parameters define the number of directions of the resulting ODFs
  // and the number of basis function kernels used for interpolation.
  typedef itk::DiffusionQballReconstructionImageFilter< 
    ReferencePixelType, GradientPixelType, OdfPrecisionType, 
    nOdfDirections, nOdfDirections > 
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
    std::cerr << "Exception detected while reading " << argv[1];
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
  }

  GradientImageType::Pointer gradImg = reader->GetOutput();

  // -------------------------------------------------------------------------
  // Parse the Nrrd headers to get the B value and the gradient directions used
  // for diffusion weighting. 
  // 
  // The Nrrd headers should look like :
  // The tags specify the B value and the gradient directions. If gradient 
  // directions are (0,0,0), it indicates that it is a reference image. 
  //
  // DWMRI_b-value:=3500
  // DWMRI_gradient_0000:= 0 0 0
  // DWMRI_gradient_0001:=-1.000000       0.000000        0.000000
  // DWMRI_gradient_0002:=-0.166000       0.986000        0.000000
  // DWMRI_gradient_0003:=0.110000        0.664000        0.740000
  // ...
  // 
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

  OutputImageType::Pointer odfImage;
  BZeroImageType::Pointer  b0Image;

  int rep = argc > 5 ? atof(argv[5]) : 1;
  for(int i=0; i<rep;i++)
  {

    // init and run the reconstruction
    QBallReconstructionImageFilterType::Pointer qBallReconstructionFilter = 
      QBallReconstructionImageFilterType::New();

    qBallReconstructionFilter->SetThreshold(static_cast< 
      QBallReconstructionImageFilterType::ReferencePixelType >( 
      atof(argv[2])));

    qBallReconstructionFilter->SetNormalizationMethod(
      QBallReconstructionImageFilterType::QBR_STANDARD);

    qBallReconstructionFilter->SetBValue(bValue);

    // An alternate way to provide the inputs, when you have the reference and
    // gradient images in seperate itk::Image< type, 3 > is  :
    // 
    //   qBallReconstructionFilter->SetReferenceImage( image0 );
    //   qBallReconstructionFilter->AddGradientImage( direction1, image1 );
    //   qBallReconstructionFilter->AddGradientImage( direction2, image2 );
    qBallReconstructionFilter->SetGradientImage( 
      gradientDirections, reader->GetOutput() );
    
    if(argc > 6)
      qBallReconstructionFilter->SetNumberOfThreads(atof(argv[6]));

    std::cout << std::endl << "This filter is using " << 
      qBallReconstructionFilter->GetNumberOfThreads() << " threads " << std::endl;

#ifdef __QBR_TIMING__
    QBallReconstructionImageFilterType::m_Collector.Start( "TotalTime" );
#endif //__QBR_TIMING__

    qBallReconstructionFilter->Update();

#ifdef __QBR_TIMING__
    QBallReconstructionImageFilterType::m_Collector.Stop( "Reconstruction" );
    QBallReconstructionImageFilterType::m_Collector.Stop( "TotalTime" );
#endif //__QBR_TIMING__

    odfImage = qBallReconstructionFilter->GetOutput();
    b0Image  = qBallReconstructionFilter->GetBZeroImage();

  }

#ifdef __QBR_TIMING__
  QBallReconstructionImageFilterType::m_Collector.Report();
#endif //__QBR_TIMING__

  // convert result to vector image that can be written to file
  // (simple image with pixeltype vector not supported by writer)
  VarVecImgType::Pointer vecImg = VarVecImgType::New();
  vecImg->SetSpacing( odfImage->GetSpacing() );
  vecImg->SetOrigin( odfImage->GetOrigin() );
  vecImg->SetDirection( odfImage->GetDirection());
  vecImg->SetLargestPossibleRegion( 
    odfImage->GetLargestPossibleRegion());
  vecImg->SetBufferedRegion( 
    odfImage->GetLargestPossibleRegion() );
  vecImg->SetVectorLength(nOdfDirections);
  vecImg->Allocate();

  itk::ImageRegionIterator<VarVecImgType> ot (
    vecImg, vecImg->GetLargestPossibleRegion() );
  ot = ot.Begin();

  itk::ImageRegionIterator<OutputImageType> it (
    odfImage, 
    odfImage->GetLargestPossibleRegion() );
  it = it.Begin();

  for (it = it.Begin(); !it.IsAtEnd(); ++it)
  {
    itk::Vector<OdfPrecisionType,nOdfDirections> vec = it.Get();
    VarVecImgType::PixelType varvec(vec.GetDataPointer(), nOdfDirections);
    ot.Set(varvec);
    ++ot;
  }

  // -------------------------------------------------------------------------
  // Write out the image of ODFs.
  OdfWriterType::Pointer odfWriter = OdfWriterType::New();
  odfWriter->SetInput(vecImg);
  ::itk::OStringStream baseName;
  if(argc < 4)
  {
    baseName << argv[1] << ".odfs.nrrd";
  }
  else
  {
    baseName << argv[3];
  }

  try
  {
    odfWriter->SetFileName( baseName.str().c_str() );
    odfWriter->Update();
  }
  catch (itk::ExceptionObject& e)
  {
    std::cerr << "Error during write of " << baseName.str();
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
  }

  // init second file writer and write extracted b-zero image to disk
  BZeroWriterType::Pointer b0Writer = BZeroWriterType::New();
  b0Writer->SetInput(b0Image);
  ::itk::OStringStream baseNameB0;
  if(argc < 5)
  {
    baseNameB0 << argv[1] << ".b0.nrrd";
  }
  else
  {
    baseNameB0 << argv[4];
  }

  try
  {
    b0Writer->SetFileName( baseNameB0.str().c_str() );
    b0Writer->Update();
  }
  catch (itk::ExceptionObject& e)
  {
    std::cerr << "Error during write of " << baseNameB0.str();
    std::cerr << " : "  << e.GetDescription();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

  
  