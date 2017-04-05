#ifndef ALICEO2_TPC_PAINTER_H_
#define ALICEO2_TPC_PAINTER_H_

///
/// \file   Painter.h
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de
///

class TH1;
class TH2;

namespace o2
{
namespace TPC 
{

//template <class T>
//class CalDet;

//template <class T>
//class CalArray;

/// \namespace Painter
/// \brief Drawing helper functions
///
/// In this namespace drawing function for calibration objects are implemented
///
/// origin: TPC
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

namespace Painter
{
  using T=float;
  /// Drawing of a CalDet object
  /// \param CalDet object to draw
  //template <class T>
  void draw(const CalDet<T>& calDet);

  /// Drawing of a CalDet object
  /// \param CalArray object to draw
  //template <class T>
  void draw(const CalArray<T>& calArray);

  /// get 2D histogram for CalDet object
  /// \param CalDet object with data
  /// \param side side which to get the histogram for
  /// \return 2D histogram with data
  //template <class T>
  TH2* getHistogram2D(const CalDet<T>& calDet, Side side);

  /// get 2D histogram for CalArray object
  /// \param CalDet object with data
  /// \param side side which to get the histogram for
  /// \return 2D histogram with data
  //template <class T>
  TH2* getHistogram2D(const CalArray<T>& calArray);

} // namespace Painter

} // namespace TPC

} // namespace AliceO2

#endif // ALICEO2_TPC_PAINTER_H_
