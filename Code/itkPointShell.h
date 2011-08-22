#ifndef __itkPointShell_h__
#define __itkPointShell_h__

namespace itk
{

#define DIST_POINTSHELL_PI 3.14159265358979323846

  /** \class PointShell
  * \brief Point distribution on a sphere by generating n-fold subdivisions 
  * of an icosahedron in case of 12, 42, 92, 162, 252, 362, 492, 642, 812 or
  * 1002 points. All other cases: computation following 
  *
  * \par References:
  * \li E. A. Rakhmanov, E. B. Saff, and Y. M. Zhou
  * MINIMAL DISCRETE ENERGY ON THE SPHERE
  * Mathematical Research Letters 1, 647–662 (1994)
  *
  * \par Template parameters
  * The class is templated over
  * \li the number of points to distribute
  * \li the output matrix type
  *
  * \par Example usage for distribution of 92 points:
  * vnl_matrix_fixed<double, 3, 92>* U =
  *   itk::PointShell<92, vnl_matrix_fixed<double, 3, 92> >
  *   ::DistributePointShell();
  *
  */
  template<int NPoints, class TMatrixType >
  class PointShell
  {
  public:
    static TMatrixType *DistributePointShell();
  };

}

#include "itkPointShell.txx"

#endif //__itkPointShell_h__
