/// \file HwFixedPoint.h
/// \brief Class to Fixed Point calculations as it would be done in Hardware
#ifndef ALICEO2_TPC_HWFIXEDPOINT_H
#define ALICEO2_TPC_HWFIXEDPOINT_H

#include "Rtypes.h"                             // for Double_t, ULong_t, etc
#include "TObject.h"
#include "TMath.h"

#include <iostream>

namespace AliceO2{
  namespace TPC{

    class HwFixedPoint : public TObject {
    public:

      typedef int T;                // base type to store value
      typedef long long int T2;     // extended type for * and / operation

      // Constructors
      HwFixedPoint(UShort_t totalWidth = 22, UShort_t decPrec = 12);
      HwFixedPoint(Int_t val, UShort_t totalWidth = 22, UShort_t decPrec = 12);
      HwFixedPoint(Float_t val, UShort_t totalWidth = 22, UShort_t decPrec = 12);
      HwFixedPoint(Double_t val, UShort_t totalWidth = 22, UShort_t decPrec = 12);
      HwFixedPoint(const HwFixedPoint& val);
      HwFixedPoint(const HwFixedPoint& val, UShort_t totalWidth, UShort_t decPrec);

      // Destructors
      ~HwFixedPoint();

      //
      // Cast operators
      //
      operator bool()       const { return (bool) mValue; };
      operator int()        const { return (int) (mValue>>mDecPrecision); };
      operator unsigned()   const { return (unsigned) (mValue>>mDecPrecision); };
      operator float()      const { return ((float) mValue) / (1<<mDecPrecision); };
      operator double()     const { return ((double) mValue) / (1<<mDecPrecision); };

      //
      // Assignment operators
      //
    template<typename TT>
      HwFixedPoint& operator =  (const TT& rhs)                   { setValue((T) (rhs*(1<<mDecPrecision))); return (*this); };
      HwFixedPoint& operator =  (const HwFixedPoint& rhs)         { mValue = rhs.mValue; mMask = rhs.mMask; mTotalWidth = rhs.mTotalWidth; 
                                                                    mDecPrecision = rhs.mDecPrecision; return (*this); };

      //
      // Addition operators
      //
      HwFixedPoint& operator += (const HwFixedPoint& rhs)         { if(rhs.mDecPrecision > mDecPrecision) { mValue += (rhs.mValue>>(rhs.mDecPrecision-mDecPrecision)); }
                                                                                                     else { mValue += (rhs.mValue<<(mDecPrecision-rhs.mDecPrecision)); };
                                                                    return (*this); };
      HwFixedPoint  operator +  (const HwFixedPoint& rhs)   const { HwFixedPoint tmp(*this); return tmp += rhs; };
    template<typename TT>
      HwFixedPoint& operator += (const TT& rhs)                   { (*this) += HwFixedPoint(rhs,mTotalWidth,mDecPrecision); return (*this); };
    template<typename TT>
      HwFixedPoint  operator +  (const TT& rhs)             const { HwFixedPoint tmp(*this); return tmp += HwFixedPoint(rhs,mTotalWidth,mDecPrecision); };
      HwFixedPoint& operator ++ ()                                { (*this) += 1; return (*this); };
      HwFixedPoint  operator ++ (int)                             { HwFixedPoint tmp(*this); ++(*this); return tmp; };

      //
      // Subtraction operators
      //
      HwFixedPoint& operator -= (const HwFixedPoint& rhs)         { if(rhs.mDecPrecision > mDecPrecision) { mValue -= (rhs.mValue>>(rhs.mDecPrecision-mDecPrecision)); }
                                                                                                     else { mValue -= (rhs.mValue<<(mDecPrecision-rhs.mDecPrecision)); };
                                                                    return (*this); };
      HwFixedPoint  operator -  (const HwFixedPoint& rhs)   const { HwFixedPoint tmp(*this); return tmp -= rhs; };
    template<typename TT>
      HwFixedPoint& operator -= (const TT& rhs)                   { (*this) -= HwFixedPoint(rhs,mTotalWidth,mDecPrecision); return (*this); };
    template<typename TT>
      HwFixedPoint  operator -  (const TT& rhs)             const { HwFixedPoint tmp(*this); return tmp -= HwFixedPoint(rhs,mTotalWidth,mDecPrecision); };
      HwFixedPoint& operator -- ()                                { (*this) -= 1; return (*this); };
      HwFixedPoint  operator -- (int)                             { HwFixedPoint tmp(*this); --(*this); return tmp; };

      //
      // Muliplication operators
      //
      HwFixedPoint& operator *= (const HwFixedPoint& rhs)         { mValue = (T) (((T2) mValue * (T2) rhs.mValue) >> rhs.mDecPrecision); 
                                                                    return (*this); };
      HwFixedPoint  operator *  (const HwFixedPoint& rhs)   const { HwFixedPoint tmp(*this); return tmp *= rhs; };
    template<typename TT>
      HwFixedPoint& operator *= (const TT& rhs)                   { (*this) *= HwFixedPoint(rhs); return (*this); };
    template<typename TT>
      HwFixedPoint  operator *  (const TT& rhs)             const {HwFixedPoint tmp(*this); return tmp *= HwFixedPoint(rhs); };

      //
      // Division operators
      //
      HwFixedPoint& operator /= (const HwFixedPoint& rhs)         { mValue = (T) ( ((T2)mValue << rhs.mDecPrecision) / ((T2) rhs.mValue) );
                                                                    return (*this); };
      HwFixedPoint  operator /  (const HwFixedPoint& rhs)   const { HwFixedPoint tmp(*this); return tmp /= rhs; };
    template<typename TT>
      HwFixedPoint& operator /= (const TT& rhs)                   { (*this) /= HwFixedPoint(rhs); return (*this); };
    template<typename TT>
      HwFixedPoint  operator /  (const TT& rhs)             const {HwFixedPoint tmp(*this); return tmp /= HwFixedPoint(rhs); };

      //
      // Comparator operators
      //
      bool          operator <  (const HwFixedPoint& rhs)   const { return mValue < rhs.mValue; };
      bool          operator >  (const HwFixedPoint& rhs)   const { return rhs < (*this); };
      bool          operator >= (const HwFixedPoint& rhs)   const { return ! ((*this) < rhs); };
      bool          operator <= (const HwFixedPoint& rhs)   const { return ! ((*this) > rhs); };
      bool          operator == (const HwFixedPoint& rhs)   const { return mValue == rhs.mValue; };
      bool          operator != (const HwFixedPoint& rhs)   const { return ! ((*this) == rhs); };
    template<typename TT>
      bool          operator <  (const TT& rhs)     const { return (*this) < HwFixedPoint(rhs); };
    template<typename TT>
      bool          operator >  (const TT& rhs)     const { return (*this) > HwFixedPoint(rhs); };
    template<typename TT>
      bool          operator >= (const TT& rhs)     const { return (*this) >= HwFixedPoint(rhs); };
    template<typename TT>
      bool          operator <= (const TT& rhs)     const { return (*this) <= HwFixedPoint(rhs); };
    template<typename TT>
      bool          operator == (const TT& rhs)     const { return (*this) == HwFixedPoint(rhs); };
    template<typename TT>
      bool          operator != (const TT& rhs)     const { return (*this) != HwFixedPoint(rhs); };
      
      //
      // Getter methods
      //
      UShort_t GetDecPrecision()    const { return mDecPrecision; };
      UShort_t GetTotalWidth()      const { return mTotalWidth; };
      T        GetValue()           const { return mValue; };
      T        GetMask()            const { return mMask; };

      //
      // Print operator
      //
    template<typename OUT>
      friend OUT& operator<< (OUT& out, const HwFixedPoint &fp) { return fp.print(out); };

    private:

      //
      // Print funciton
      //
      template<typename OUT>
        OUT& print(OUT& out) const { return ( out << ((double)mValue)/(1<<mDecPrecision)
                                                  << " <" << mTotalWidth
                                                  << "," << mDecPrecision << ">"// ("
                                                  //<< mValue << ")"
                                                  //<< "(" << (sizeof(T)*8)-mTotalWidth-mDecPrecision
                                                  //<< " bits not used)"
                                                  ); };

      void setValue(T val);

      T mValue;
      T mMask;
      UShort_t mDecPrecision;
      UShort_t mTotalWidth;

      ClassDef(HwFixedPoint, 1);
    };
  }
}
#endif

