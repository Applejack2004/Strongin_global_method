/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      evolvent.h                                                  //
//                                                                         //
//  Purpose:   Header file for evolvent classes                            //
//                                                                         //
//  Author(s): Barkalov K., Sysoyev A.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __EVOLVENT_H__
#define __EVOLVENT_H__

#include "common.h"
#include "extended.h"
#include <vector>

// ------------------------------------------------------------------------------------------------
class TEvolvent
{
protected:
  int      m;             // accuracy of decomposition of hypercube
  int      N;             // dimension
  double   A[MaxDim];     // left and
  double   B[MaxDim];     // right bounds of search area

  double*  y;             // y point from hypercube [-1/2, 1/2]^N

  const Extended extNull; // = 0.0; //Extended(0.0);
  const Extended extOne;  // = 1.0; //Extended(1.0);
  const Extended extHalf; // = 0.5; //Extended(0.5);
  Extended nexpExtended;

  void CalculateNumbr(Extended *s, int *u, int *v, int *l);
  void CalculateNode(Extended is, int n, int *u, int *v, int *l);
  void transform_P_to_D(); // transformation from hypercube P to hyperinterval D
  void transform_D_to_P(); // transformation from hyperinterval D to hypercube P
  double* GetYOnX(const Extended& _x);
  Extended GetXOnY();
public:
  TEvolvent(int _N = 2, int _m = 10);
  TEvolvent(const TEvolvent& evolvent);
  virtual ~TEvolvent();
  void SetBounds(const double* _A, const double* _B);
  //x-->y
  virtual void GetImage(const Extended& x, double* _y, int EvolventNum = 0);
  //y-->x
  void GetInverseImage(double* _y, Extended& x);
  //y-->x
  virtual void GetPreimages(double* _y, Extended* x);
  TEvolvent& operator=(const TEvolvent& evolvent);

  /// Вычисляет функцию существования точки в развертки EvolventNum для y, <0 - существует
  virtual double ZeroConstraintCalc(const double* _y, int EvolventNum = 0);

};

class TShiftedEvolvent : public TEvolvent
{
protected:
  int    L; // number of evolvents - 1, 0 by default, means one evolvent is used
  double PowOf2[MaxDim * MaxL + 1]; // degrees of 2
  void transform_P_to_Pl(int EvolventNum); //
  void transform_Pl_to_P(int EvolventNum);
  double ZeroConstraint();
public:
  TShiftedEvolvent(int _N = 2, int _m = 10, int _L = 0);
  virtual ~TShiftedEvolvent();
  virtual void GetImage(const Extended& x, double* _y, int EvolventNum = 0);
  virtual void GetPreimages(double* _y, Extended* x);
  /// Вычисляет функцию существования точки в развертки EvolventNum для y, <0 - существует
  virtual double ZeroConstraintCalc(const double* _y, int EvolventNum = 0);
};

class TRotatedEvolvent : public TEvolvent
{
protected:
  int L;          // number of evolvents, 1 by default, means one evolvent is used
  int PlaneCount; // current number of planes
  int Planes[MaxDim * (MaxDim - 1) / 2][2]; // axes array
  double PowOfHalf[MaxM + 2]; // degrees of 1/2

  void GetAllPlanes();
public:
  TRotatedEvolvent(int _N = 2, int _m = 10, int _L = 0);
  virtual ~TRotatedEvolvent();
  virtual void GetImage(const Extended& x, double* _y, int EvolventNum = 0);
  virtual void GetPreimages(double* _y, Extended* x);
};

class TSmoothEvolvent : public TEvolvent
{
protected:
  double h; //smoothness
  int smoothPointCount;
  bool continuously;
  std::vector<double> tmp_y, tmp_y_;
  void GetYOnXSmooth(double x, std::vector<double>&y, std::vector<double>&y_);
public:
  TSmoothEvolvent(int _N = 2, int _m = 10, double _h=0.25);
  virtual ~TSmoothEvolvent();
  virtual void GetImage(const Extended& x, double* _y, int EvolventNum = 0);
  void GetInverseImage(double* _y, Extended& x);
};

class TLinearEvolvent : public TEvolvent
{
protected:
  Extended mneExtended;
public:
  TLinearEvolvent(int _N = 2, int _m = 10);
  virtual ~TLinearEvolvent();
  virtual void GetImage(const Extended& x, double* _y, int EvolventNum = 0);
  void GetInverseImage(double* _y, Extended& x);
};

class TNoninjectiveEvolvent : public TEvolvent
{
protected:
  int max_preimages;
  Extended mneExtended;
public:
  TNoninjectiveEvolvent(int _N = 2, int _m = 10, int _max_preimages=64);
  virtual ~TNoninjectiveEvolvent();
  virtual void GetImage(const Extended& x, double* _y, int EvolventNum = 0);
  virtual int GetNoninjectivePreimages(double* _y, Extended* x);
  int GetMaxPreimagesNumber() const { return max_preimages; };
};

#endif
// - end of file ----------------------------------------------------------------------------------
