/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      evolvent.cpp                                                //
//                                                                         //
//  Purpose:   Source file for evolvent classes                            //
//                                                                         //
//  Author(s): Barkalov K., Sysoyev A.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
#pragma warning(disable:4996) 

#include "evolvent.h"
#include "exception.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// ------------------------------------------------------------------------------------------------
TEvolvent::TEvolvent(int _N, int _m) : extNull(0.0), extOne(1.0), extHalf(0.5)
{
  int i;
  if ((_N < 1) || (_N > MaxDim))
  {
    throw EXCEPTION("N is out of range");
  }
  N = _N;
  y = new double[N];

  if ((_m < 2) || (_m > MaxM))
  {
    throw EXCEPTION("m is out of range");
  }
  m = _m;

  for (nexpExtended = extOne, i = 0; i < N; nexpExtended += nexpExtended, i++)
    ;
}

// ------------------------------------------------------------------------------------------------
TEvolvent::~TEvolvent()
{
  delete[] y;
}
// ------------------------------------------------------------------------------------------------
TEvolvent::TEvolvent(const TEvolvent& evolvent) : extNull(0.0), extOne(1.0), extHalf(0.5)
{
  //Считаем развертку evolvent корректной, проверка ее параметров не требуется
  N = evolvent.N;
  y = new double[N];
  m = evolvent.m;
  nexpExtended = evolvent.nexpExtended;

  for (int i = 0; i < N; i++)
  {
    A[i] = evolvent.A[i];
    B[i] = evolvent.B[i];
  }
}

// ------------------------------------------------------------------------------------------------
TEvolvent& TEvolvent::operator=(const TEvolvent& evolvent)
{
  //Считаем развертку evolvent корректной, проверка ее параметров не требуется
  if (N != evolvent.N)
  {
    N = evolvent.N;
    if (y)
    {
      delete[] y;
    }
    y = new double[N];
  }
  m = evolvent.m;
  nexpExtended = evolvent.nexpExtended;

  for (int i = 0; i < N; i++)
  {
    A[i] = evolvent.A[i];
    B[i] = evolvent.B[i];
  }

  return *this;
}


// ----------------------------------------------------------------------------
void TEvolvent::SetBounds(const double* _A, const double* _B)
{
  for (int i = 0; i < N; i++)
  {
    A[i] = _A[i];
    B[i] = _B[i];
  }
}

// ----------------------------------------------------------------------------
void TEvolvent::CalculateNumbr(Extended *s, int *u, int *v, int *l)
// calculate s(u)=is,l(u)=l,v(u)=iv by u=iu
{
  int i, k1, k2, l1;
  Extended is, iff;

  iff = nexpExtended;
  is = extNull;
  k1 = -1;
  k2 = 0;
  l1 = 0;
  for (i = 0; i < N; i++)
  {
    iff = iff / 2;
    k2 = -k1 * u[i];
    v[i] = u[i];
    k1 = k2;
    if (k2 < 0)
      l1 = i;
    else
    {
      is += iff;
      *l = i;
    }
  }
  if (is == extNull)
    *l = N - 1;
  else
  {
    v[N - 1] = -v[N - 1];
    if (is == (nexpExtended - extOne))
      *l = N - 1;
    else
    {
      if (l1 == (N - 1))
        v[*l] = -v[*l];
      else
        *l = l1;
    }
  }
  *s = is;
}

// ----------------------------------------------------------------------------
void TEvolvent::CalculateNode(Extended is, int n, int *u, int *v, int *l)
// вычисление вспомогательного центра u(s) и соответствующих ему v(s) и l(s)
// calculate u=u[s], v=v[s], l=l[s] by is=s
{
  int n1, i, j, k1, k2, iq;
  Extended iff;

  iq = 1;
  n1 = n - 1;
  *l = 0;
  if (is == 0)
  {
    *l = n1;
    for (i = 0; i < n; i++)
    {
      u[i] = -1;
      v[i] = -1;
    }
  }
  else if (is == (nexpExtended - extOne))
  {
    *l = n1;
    u[0] = 1;
    v[0] = 1;
    for (i = 1; i < n; i++)
    {
      u[i] = -1;
      v[i] = -1;
    }
    v[n1] = 1;
  }
  else
  {
    iff = nexpExtended;
    k1 = -1;
    for (i = 0; i < n; i++)
    {
      iff = iff / 2;
      if (is >= iff)
      {
        if ((is == iff) && (is != extOne))
        {
          *l = i;
          iq = -1;
        }
        is -= iff;
        k2 = 1;
      }
      else
      {
        k2 = -1;
        if ((is == (iff - extOne)) && (is != extNull))
        {
          *l = i;
          iq = 1;
        }
      }
      j = -k1 * k2;
      v[i] = j;
      u[i] = j;
      k1 = k2;
    }
    v[*l] = v[*l] * iq;
    v[n1] = -v[n1];
  }
}

// ------------------------------------------------------------------------------------------------
void TEvolvent::transform_P_to_D()
{
  //if (N == 1) return;
  // transformation from hypercube P to hyperinterval D
  for (int i = 0; i < N; i++)
    y[i] = y[i] * (B[i] - A[i]) + (A[i] + B[i]) / 2;
}

// ----------------------------------------------------------------------------
void TEvolvent::transform_D_to_P()
{
  //if (N == 1) return;
  // transformation from hyperinterval D to hypercube P
  for (int i = 0; i < N; i++)
    y[i] = (y[i] - (A[i] + B[i]) / 2) / (B[i] - A[i]);
}

// ----------------------------------------------------------------------------
double* TEvolvent::GetYOnX(const Extended& _x)
{
  if (N == 1)
  {
    y[0] = _x.toDouble() - 0.5;
    return y;
  }

  int iu[MaxDim];
  int iv[MaxDim];
  int l;
  Extended d;
  int mn;
  double r;
  int iw[MaxDim];
  int it, i, j;
  Extended is;

  d = _x;
  r = 0.5;
  it = 0;
  mn = m * N;
  for (i = 0; i < N; i++)
  {
    iw[i] = 1;
    y[i] = 0.0;
  }
  for (j = 0; j < m; j++)
  {
    if (_x == extOne)
    {
      is = nexpExtended - extOne;
      d = extNull; // d = 0.0;
    }
    else
    {
      //Код из старой версии - уточнить работоспособность при N > 32
      d *= nexpExtended;
      is = (int)d.toDouble();
      d -= is;
    }
    CalculateNode(is, N, iu, iv, &l);
    i = iu[0];
    iu[0] = iu[it];
    iu[it] = i;
    i = iv[0];
    iv[0] = iv[it];
    iv[it] = i;
    if (l == 0)
      l = it;
    else if (l == it)
      l = 0;
    r *= 0.5;
    it = l;
    for (i = 0; i < N; i++)
    {
      iu[i] *= iw[i];
      iw[i] *= -iv[i];
      y[i] += r * iu[i];
    }
  }
  return y;
}

//-----------------------------------------------------------------------------
Extended TEvolvent::GetXOnY()
{
  int u[MaxDim], v[MaxDim];
  Extended x, r1;
  if (N == 1)
  {
    x = y[0] + 0.5;
    return x;
  }

  double  r;
  int w[MaxDim];
  int i, j, it, l;
  Extended is;

  for (i = 0; i < N; i++)
    w[i] = 1;
  r = 0.5;
  r1 = extOne;
  x = extNull;
  it = 0;
  for (j = 0; j < m; j++)
  {
    r *= 0.5;
    for (i = 0; i < N; i++)
    {
      u[i] = (y[i] < 0) ? -1 : 1;
      y[i] -= r * u[i];
      u[i] *= w[i];
    }
    i = u[0];
    u[0] = u[it];
    u[it] = i;
    CalculateNumbr(&is, u, v, &l);
    i = v[0];
    v[0] = v[it];
    v[it] = i;
    for (i = 0; i < N; i++)
      w[i] *= -v[i];
    if (l == 0)
      l = it;
    else
      if (l == it)
        l = 0;
    it = l;
    r1 = r1 / nexpExtended;
    x += r1 * is;
  }
  return x;
}

//----------------------------------------------------------------------------
void TEvolvent::GetImage(const Extended& x, double* _y, int EvolventNum)
{

  // в одиночной развертке EvolventNum не используется
  //   введен, чтобы работал полиморфизм в множественных развертках


  if ((x.toDouble() < 0) || (x.toDouble() > 1))
  {
    throw EXCEPTION("x is out of range");
  }
  // x ---> y
  GetYOnX(x); // it saves return value to y, so no need to call operator= again

  transform_P_to_D();

  memcpy(_y, y, N * sizeof(double));
}

void TEvolvent::GetInverseImage(double* _y, Extended& x)
{
  // y ---> x
  memcpy(y, _y, N * sizeof(double));
  transform_D_to_P();
  x = GetXOnY();
}

//----------------------------------------------------------------------------
void TEvolvent::GetPreimages(double* _y, Extended* x)
{
  // y ---> x
  memcpy(y, _y, N * sizeof(double));
  transform_D_to_P();
  x[0] = GetXOnY();
}

// ------------------------------------------------------------------------------------------------
/// Вычисляет функцию существования точки в развертки EvolventNum для y, <0 - существует
double TEvolvent::ZeroConstraintCalc(const double* _y, int EvolventNum)
{
  return -1;
}

// ------------------------------------------------------------------------------------------------
TShiftedEvolvent::TShiftedEvolvent(int _N, int _m, int _L) :
  TEvolvent(_N, _m)
{
  if ((_L < 0) || (_L >= m))
  {
    throw EXCEPTION("L is out of range");
  }
  L = _L;
  //Инициализация массива степеней двойки
  PowOf2[0] = 1;
  for (int i = 1; i <= L * N; i++)
	  PowOf2[i] = PowOf2[i - 1] * 2;

}

// ------------------------------------------------------------------------------------------------
TShiftedEvolvent::~TShiftedEvolvent()
{
}

// ------------------------------------------------------------------------------------------------
void TShiftedEvolvent::GetImage(const Extended& x, double* _y, int EvolventNum)
{
  /*
  double* tmp = NULL;
  int i;
  //непонятно, зачем в GlobalExpert данная точка была вынесена в отдельную ветку

  if (x == extNull)
  {
    for (i = 0; i < N; i++)
    {
      _y[i] = 0.0;
    }
  }
  else
    */
  GetYOnX(x);

  transform_P_to_Pl(EvolventNum);
  transform_P_to_D();

  memcpy(_y, y, N * sizeof(double));
}

// ------------------------------------------------------------------------------------------------
void TShiftedEvolvent::transform_P_to_Pl(int EvolventNum)
{
  //  if (N == 1) return;
    // transformation from hypercube P to hypercube P[l]
  double temp;
  if (EvolventNum == 0)
  {
    temp = 0.0;
  }
  else
  {
    temp = 1.0 / PowOf2[EvolventNum]; // temp = 1 / 2^l (l = 1,...,L)
  }
  for (int i = 0; i < N; i++)
  {
    y[i] = y[i] * 2 + 0.5 - temp;
  }
}

// ------------------------------------------------------------------------------------------------
void TShiftedEvolvent::transform_Pl_to_P(int EvolventNum)
{
  //  if (N == 1) return;
    // transformation from hypercube P to hypercube P[l]
  double temp;
  if (EvolventNum == 0)
  {
    temp = 0;
  }
  else
  {
    temp = 1.0 / PowOf2[EvolventNum]; // temp = 1 / 2^l (l = 1,...,L)
  }
  for (int i = 0; i < N; i++)
  {
    y[i] = (y[i] - 0.5 + temp) / 2;
  }
}

// ------------------------------------------------------------------------------------------------
double TShiftedEvolvent::ZeroConstraint()
{
  double CurZ = -MaxDouble;
  for (int i = 0; i < N; i++)
  {
    if (fabs(y[i]) - 0.5 > CurZ)
    {
      CurZ = fabs(y[i]) - 0.5;
    }
  }
  return CurZ;
}


// ------------------------------------------------------------------------------------------------
void TShiftedEvolvent::GetPreimages(double* _y, Extended *x)
{
  for (int i = 0; i <= L; i++)
  {
    memcpy(y, _y, N * sizeof(double));
    transform_D_to_P();
    transform_Pl_to_P(i);
    x[i] = GetXOnY();
  }

}

// ------------------------------------------------------------------------------------------------
/// Вычисляет функцию существования точки в развертки EvolventNum для y, <0 - существует
double TShiftedEvolvent::ZeroConstraintCalc(const double* _y, int EvolventNum)
{
  // копируем y
  memcpy(y, _y, N * sizeof(double));
  // центрируем и нормируем область
  transform_D_to_P();
  // сдвигаем область в соответствии с разверткой
  transform_Pl_to_P(EvolventNum);
  // вычисляем функционал
  return ZeroConstraint();
}

// ------------------------------------------------------------------------------------------------
TRotatedEvolvent::TRotatedEvolvent(int _N, int _m, int _L) :
  TEvolvent(_N, _m)
{
  L = _L;
  // !!!!!!!!!!!!!
  if (N == 1)
    return;
  // !!!!!!!!!!!!!
  PlaneCount = N * (N - 1) / 2;
  if ((L < 1) || (L > 2 * PlaneCount + 1))
  {
    throw EXCEPTION("L is out of range");
  }
  GetAllPlanes();
  PowOfHalf[0] = 1;
  for (int i = 1; i < m + 2; i++)
    PowOfHalf[i] = PowOfHalf[i - 1] / 2;
}

// ------------------------------------------------------------------------------------------------
TRotatedEvolvent::~TRotatedEvolvent()
{
}

// ------------------------------------------------------------------------------------------------
void TRotatedEvolvent::GetAllPlanes()
{
  const int k = 2; // Подмножества из двух элементов
  int plane[k];    // Два номера под элементы

  for (int i = 0; i < k; i++)
    plane[i] = i;

  if (N <= k)
  {
    for (int i = 0; i < k; i++)
    {
      Planes[0][i] = plane[i];
    }
    return;
  }
  int p = k - 1;
  int counter = 0; //счетчик числа перестановок
  while (p >= 0)
  {
    for (int i = 0; i < k; i++)
    {
      Planes[counter][i] = plane[i];
    }
    counter++;

    if (plane[k - 1] == N - 1)
    {
      p--;
    }
    else
    {
      p = k - 1;
    }

    if (p >= 0)
    {
      for (int i = k - 1; i >= p; i--)
      {
        plane[i] = plane[p] + i - p + 1;
      }
    }
  }
}

// ------------------------------------------------------------------------------------------------
void TRotatedEvolvent::GetImage(const Extended& x, double* _y, int EvolventNum)
{
  if (L == 1 || EvolventNum == 0)
  {
    TEvolvent::GetImage(x, _y);
    return;
  }

  int PlaneIndex = EvolventNum - 1; // теперь PlaneIndex - номер перестановки
  PlaneIndex = PlaneIndex % PlaneCount;

  GetYOnX(x);

  // shift to center for convenient rotation
  //for (int i = 0; i < N; i++)
  //  y[i] += PowOfHalf[m + 1];

  // rotate
  double tmpCoord = y[Planes[PlaneIndex][1]];
  y[Planes[PlaneIndex][1]] = y[Planes[PlaneIndex][0]];
  y[Planes[PlaneIndex][0]] = -tmpCoord;

  //Меняем знак преобразования, если число разверток больше числа плоскостей
  if (EvolventNum > PlaneCount)
  {
    y[Planes[PlaneIndex][0]] = -y[Planes[PlaneIndex][0]];
    y[Planes[PlaneIndex][1]] = -y[Planes[PlaneIndex][1]];
  }

  // shift back to corner
  //for (int i = 0; i < N; i++)
  //  y[i] -= PowOfHalf[m + 1];

  transform_P_to_D();
  memcpy(_y, y, N * sizeof(double));
}

// ------------------------------------------------------------------------------------------------
void TRotatedEvolvent::GetPreimages(double* _y, Extended *x)
{
  memcpy(y, _y, N * sizeof(double));
  transform_D_to_P();
  // прообраз для первой развертки
  x[0] = GetXOnY();

  if (L == 1)
    return;

  for (int i = 1; i < L; i++)
  {
    memcpy(y, _y, N * sizeof(double));
    transform_D_to_P();
    // обратное преобразование координат
    int PlaneIndex = (i - 1) % PlaneCount;

    double tmpCoord = y[Planes[PlaneIndex][1]];
    y[Planes[PlaneIndex][1]] = -y[Planes[PlaneIndex][0]];
    y[Planes[PlaneIndex][0]] = tmpCoord;

    if (i > PlaneCount)//Меняем знак преобразования, если число разверток больше числа плоскостей
    {
      y[Planes[PlaneIndex][0]] = -y[Planes[PlaneIndex][0]];
      y[Planes[PlaneIndex][1]] = -y[Planes[PlaneIndex][1]];
    }

    // прообраз для i - 1 развертки
    x[i] = GetXOnY();
  }
}

namespace
{
// ------------------------------------------------------------------------------------------------
int _Pow_int(int x, int n)
{
    int val = 1;
    for(int i = 0; i < n; i++)
        val *= x;
    return val;
}

// ------------------------------------------------------------------------------------------------
double Hermit(double y0, double d0, double y1, double d1, double h, double x)
{
  return y0 + (x + h) * (d0 + (x + h) * (d0 - (y0 - y1) / (-2 * h) +
        (x - h) * (d0 - 2 * (y0 - y1) / (-2 * h) + d1)/(-2 * h)) / (-2 * h));
}

// ------------------------------------------------------------------------------------------------
double HermitDer(double y0, double d0, double y1, double d1, double h, double x)
{
  return pow(h, -3) * (d0 * h * (-0.25 * h * h - 0.5 * h * x + 0.75 * x * x)
                    + d1 * h * (-0.25 * h * h + 0.5 * h * x + 0.75 * x * x)
                    - 0.75 * h * h * y0 + 0.75 * h * h * y1
                    + 0.75 * x * x * y0
                    - 0.75 * x * x * y1);
}

// ------------------------------------------------------------------------------------------------
int node_smooth(int is, int n, int& iq, int nexp, int* iu, int* iv)
{
 /* calculate iu=u[s], iv=v[s], l=l[s] by is=s */
  iq = 0;
  int k1, k2, iff;
  static int l = 0;
  if (is == 0)
  {
    l = n - 1;
    for (int i = 0; i < n; i++)
    {
      iu[i] = -1;
      iv[i] = -1;
    }
  }
  else
    if (is == (nexp - 1))
    {
      l = n - 1;
      iu[0] = 1;
      iv[0] = 1;
      for (int i = 1; i < n; i++)
      {
        iu[i] = -1;
        iv[i] = -1;
      }
      iv[n - 1] = 1;
    }
    else
    {
      iff = nexp;
      k1 = -1;
      for (int i = 0; i < n; i++)
      {
        iff /= 2;
        if (is >= iff)
        {
          if ((is == iff) && (is != 1))  { l = i; iq = -1; }
          is = is - iff;
          k2 = 1;
        }
        else
        {
          k2 = -1;
          if ((is == (iff - 1)) && (is != 0))  { l = i; iq = 1; }
        }
        iu[i] = iv[i] = -k1 * k2;
        k1 = k2;
      }
      iv[l] *= iq;
      iv[n - 1] = -iv[n - 1];
    }
    return l;
}

// ------------------------------------------------------------------------------------------------
void SmoothEvolventDer(double x, int n, int m, std::vector<double>&y, std::vector<double>&y_, bool c)
{
  if (y.size() != n)
    y.assign(n, .0);
  if (y_.size() != n)
    y_.assign(n, .0);
  int l = 0, iq = 0;
  std::vector<int> iu(n, 0);
  std::vector<int> iv(n, 0);
  int nexp = _Pow_int(2, n); // nexp=2**n */
  double mnexp = _Pow_int(nexp, m); // mnexp=2**(nm)
  double d = 1.0 / mnexp;
  std::vector<int> iw(n, 0);
  double xd = x;
  int it = 0;
  double dr = nexp;
  for (int i = 0; i < n; i++)
  {
    iw[i] = 1;
    y[i] = 0;
  }
  int k = 0;
  double r = 0.5;
  int ic;
  for (int j = 0; j < m; j++)
  {
    if (x == 1.0 - d)
    {
      ic = nexp - 1;
      xd = 0.0;
    }
    else
    {
      xd = xd * nexp;
      ic = (int)xd;
      xd = xd - ic;
    }
    iq = 0;
    l = node_smooth(ic, n, iq, nexp, iu.data(), iv.data());
    int swp = iu[it];
    iu[it] = iu[0];
    iu[0] = swp;
    swp = iv[it];
    iv[it] = iv[0];
    iv[0] = swp;
    if (l == 0)
      l = it;
    else if (l == it)
      l = 0;
    if ((iq > 0) || ((iq == 0) && (ic == 0)))
      k = l;
    else if (iq < 0)
    {
      if (it == n - 1)
        k = 0;
      else
        k = n - 1;
    }
    r *= 0.5;
    it = l;
    for (int i = 0; i < n; i++)
    {
      iu[i] *= iw[i];
      iw[i] *= -iv[i];
      y[i] += r * iu[i];
    }
  }
  if (c)
  {
    if (ic == (nexp - 1))
    {
      y[k] += 2 * iu[k] *r * xd;
      y_[k] = ((iu[k] > 0) ? 1 : (iu[k] < 0) ? -1 : 0) * pow(2, m * (n - 1));
    }
    else
    {
      y[k] -= 2 * iu[k] * r * xd;
      y_[k] = - ((iu[k] > 0) ? 1 : (iu[k] < 0) ? -1 : 0) * pow(2, m * (n - 1));
    }
    if (x == 1.0 - d)
    {
      std::vector<double> y0(n, 0);
      y_[k] = 0;
      ::SmoothEvolventDer(x - d / 2, n, m, y0, y_, true);
    }
  }
}
}

// ------------------------------------------------------------------------------------------------
TSmoothEvolvent::TSmoothEvolvent(int _N, int _m, double _h) :
  TEvolvent(_N, _m)
{
  if (_N > 2)
  {
    print << "Warning: smooth evolvent is very slow when problem dimension > 2\n";
  }
  h = _h;
  if (h < 0 || h > 1)
  {
    throw EXCEPTION("h is out of range");
  }
  continuously = h != 1. ?  true : false;
  smoothPointCount = 0;
  tmp_y.resize(N);
  tmp_y_.resize(N);
}

// ------------------------------------------------------------------------------------------------
void TSmoothEvolvent::GetImage(const Extended& x, double* _y, int EvolventNum)
{
  if ((x.toDouble() < 0) || (x.toDouble() > 1))
  {
    throw EXCEPTION("x is out of range");
  }
  // x ---> y
  GetYOnXSmooth(x.toDouble(), tmp_y, tmp_y_);
  std::copy(tmp_y.begin(), tmp_y.end(), y);
  transform_P_to_D();

  memcpy(_y, y, N * sizeof(double));
}

// ------------------------------------------------------------------------------------------------
void TSmoothEvolvent::GetInverseImage(double* _y, Extended& x)
{
  throw EXCEPTION("This method is not implemented for the smooth evolvent");
}

// ------------------------------------------------------------------------------------------------
void TSmoothEvolvent::GetYOnXSmooth(double x, std::vector<double>&y, std::vector<double>&y_)
{
  y.assign(N, .0);
  y_.assign(N, .0);

  if (N == 1)
  {
    y[0] = x - 0.5;
    return;
  }

  int l = 0, iq = 0;
  std::vector<int> iu(N, 0);
  std::vector<int> iv(N, 0);
  int nexp = _Pow_int(2, N); // nexp=2**n
  double mnexp = pow(nexp, m); // mnexp=2**(nm)
  double d = 1.0 / mnexp;
  double dh = d * h;
  double xc = 0;
  while (xc < x)
    xc += d;
  if (((h > 0) && (h <= .5) && (x > dh) && (x < 1 - d - dh) && ((xc - x < dh) || (xc - x > d - dh))) && continuously)
  {
    smoothPointCount++;

    std::vector<double> y0(N, 0);
    std::vector<double> y1(N, 0);
    std::vector<double> y0_(N, 0);
    std::vector<double> y1_(N, 0);
    double xh = 0;
    if (xc - x < dh)
    {
      xh += x - xc;
      SmoothEvolventDer(xc - dh, N, m, y0, y0_, true);
      SmoothEvolventDer(xc + dh, N, m, y1, y1_, true);
    }
    else
    {
      xh += x - xc + d;
      SmoothEvolventDer(xc - d - dh, N, m, y0, y0_, true);
      SmoothEvolventDer(xc - d + dh, N, m, y1, y1_, true);
    }
    int i0 = -1, i1 = -1;
    for(int i = 0; i < N; i++)
    {
      if (y1[i] != y0[i])
      {
        if (i0 == -1)
          i0 = i;
        else
          i1 = i;
      }
      else
      {
        y[i] = y0[i];
        y_[i] = 0;
      }
    }
    if (i0 != -1 && i1 != -1)
    {
      y[i0] = Hermit(y0[i0], y0_[i0], y1[i0], y1_[i0], dh, xh);
      y_[i0] = HermitDer(y0[i0], y0_[i0], y1[i0], y1_[i0], dh, xh);
      y[i1] = Hermit(y0[i1], y0_[i1], y1[i1], y1_[i1], dh, xh);
      y_[i1] = HermitDer(y0[i1], y0_[i1], y1[i1], y1_[i1], dh, xh);
      return;
    }
  }
  SmoothEvolventDer(x, N, m, y, y_, continuously);
}

// ------------------------------------------------------------------------------------------------
TSmoothEvolvent::~TSmoothEvolvent()
{}

namespace
{
// ------------------------------------------------------------------------------------------------
  void numbr(Extended *iss, const int n1, const Extended& nexp,
             int& l, int* iu, int* iv,
             const Extended& extOne,
             const Extended& extZero,
             const Extended& extHalf)
  {
    /* calculate s(u)=is,l(u)=l,v(u)=iv by u=iu */

    Extended iff, is;
    int n, k1, k2, l1;

    n = n1 + 1;
    iff = nexp;
    is = extZero;
    k1 = -1;
    for (int i = 0; i < n; i++)
    {
      iff = iff * extHalf;
      k2 = -k1 * iu[i];
      iv[i] = iu[i];
      k1 = k2;
      if (k2 < 0) l1 = i;
      else { is += iff; l = i; }
    }
    if (is == extZero) l = n1;
    else
    {
      iv[n1] = -iv[n1];
      if (is == (nexp - extOne)) l = n1;
      else if (l1 == n1) iv[l] = -iv[l];
      else l = l1;
    }
    *iss = is;
  }

// ------------------------------------------------------------------------------------------------
  void xyd(Extended *xx, int m, double y[], int n,
           const Extended& nexp,
           const Extended& extOne,
           const Extended& extZero,
           const Extended& extHalf)
  {
    /* calculate preimage x  for nearest level  m center to y */
    /* (x - left boundary point of level m interval)          */
    int n1, l, iu[MaxDim], iv[MaxDim];

    Extended r1, x;
    double r;
    int iw[MaxDim + 1];
    int it;
    Extended is;

    n1 = n - 1;
    for (int i = 0; i < n; i++)
    {
      iw[i] = 1;
    }
    r = 0.5;
    r1 = extOne;
    x = extZero;
    it = 0;
    for (int j = 0; j < m; j++)
    {
      r *= 0.5;
      for (int i = 0; i < n; i++)
      {
        iu[i] = (y[i] < 0) ? -1 : 1;
        y[i] -= r * iu[i];
        iu[i] *= iw[i];
      }
      std::swap(iu[0], iu[it]);
      numbr(&is, n1, nexp, l, iu, iv, extOne, extZero, extHalf);
      std::swap(iv[0], iv[it]);
      for (int i = 0; i < n; i++)
        iw[i] = -iw[i] * iv[i];
      if (l == 0) l = it;
      else if (l == it) l = 0;
      it = l;
      r1 = r1 / nexp;
      x += r1 * is;
    }
    *xx = x;
  }

// ------------------------------------------------------------------------------------------------
  void invmad(int m, Extended xp[], int kp,
    int *kxx, double p[], int n, int incr,
    const Extended& nexp,
    const Extended& mne,
    const Extended& extOne,
    const Extended& extZero,
    const Extended& extHalf)
  {
    /* calculate kx preimage p node */
    /*   node type mapping m level  */

    Extended dr, dd, del, d1, x;
    double r, d, u[MaxDim], y[MaxDim];
    int i, k, kx;

    kx = 0;
    kp--;
    for (int i = 0; i < n; i++)
    {
      u[i] = -1.0;
    }
    dr = nexp;
    for (r = 0.5, i = 0; i < m; i++)
    {
      r *= 0.5;
    }
    dr = mne / nexp;

    dr = dr - fmod(dr, extOne);
    //dr = (dr>0) ? floor(dr) : ceil(dr);

    del = extOne / (mne - dr);
    d1 = del * (incr + 0.5);
    for (kx = -1; kx < kp;)
    {
      for (i = 0; i < n; i++)
      {       /* label 2 */
        d = p[i];
        y[i] = d - r * u[i];
      }
      for (i = 0; (i < n) && (fabs(y[i]) < 0.5); i++);
      if (i >= n)
      {
        xyd(&x, m, y, n, nexp, extOne, extZero, extHalf);
        dr = x * mne;
        dd = dr - fmod(dr, extOne);
        //dd = (dr>0) ? floor(dr) : ceil(dr);
        dr = dd / nexp;
        dd = dd - dr + fmod(dr, extOne);
        //dd = dd - ((dr>0) ? floor(dr) : ceil(dr));
        x = dd * del;
        if (kx > kp) break;
        k = kx++;                     /* label 9 */
        if (kx == 0) xp[0] = x;
        else
        {
          while (k >= 0)
          {
            dr = fabs(x - xp[k]);     /* label 11 */
            if (dr <= d1) {
              for (kx--; k < kx; k++, xp[k] = xp[k + 1]);
              goto m6;
            }
            else
              if (x <= xp[k])
              {
                xp[k + 1] = xp[k]; k--;
              }
              else break;
          }
          xp[k + 1] = x;
        }
      }
    m6: for (i = n - 1; (i >= 0) && (u[i] = (u[i] <= 0.0) ? 1 : -1) < 0; i--);
    if (i < 0) break;
    }
    *kxx = ++kx;

  }

// ------------------------------------------------------------------------------------------------
  void node(Extended is, int n1, Extended nexp, int& l, int& iq, int iu[], int iv[],
            const Extended& extOne,
            const Extended& extZero,
            const Extended& extHalf)
  {
    /* calculate iu=u[s], iv=v[s], l=l[s] by is=s */

    Extended iff;

    int n = n1 + 1;
    if (is == extZero)
    {
      l = n1;
      for (int i = 0; i < n; i++)
      {
        iu[i] = -1; iv[i] = -1;
      }
    }
    else if (is == (nexp - extOne))
    {
      l = n1;
      iu[0] = 1;
      iv[0] = 1;
      for (int i = 1; i < n; i++) {
        iu[i] = -1; iv[i] = -1;
      }
      iv[n1] = 1;
    }
    else
    {
      iff = nexp;
      int k1 = -1, k2;
      for (int i = 0; i < n; i++)
      {
        iff = iff / 2;
        if (is >= iff) {
          if ((is == iff) && (is != extOne)) { l = i; iq = -1; }
          is = is - iff;
          k2 = 1;
        }
        else
        {
          k2 = -1;
          if ((is == (iff - extOne)) && (is != extZero)) { l = i; iq = 1; }
        }
        int j = -k1 * k2;
        iv[i] = j;
        iu[i] = j;
        k1 = k2;
      }
      iv[l] = iv[l] * iq;
      iv[n1] = -iv[n1];
    }
  }

// ------------------------------------------------------------------------------------------------
  void mapd(Extended x, int m, double* y, int n, int key,
            const Extended& nexp,
            const Extended& mne,
            const Extended& extOne,
            const Extended& extZero,
            const Extended& extHalf)
  {
    /* mapping y(x) : 1 - center, 2 - line, 3 - node */
    // use key = 1

    int n1, l, iq, iu[MaxDim], iv[MaxDim];
    Extended d, is;
    double p, r;
    int iw[MaxDim];
    int it, k;

    p = 0.0;
    n1 = n - 1;
    d = x;
    r = 0.5;
    it = 0;
    for (int i = 0; i < n; i++)
    {
      iw[i] = 1; y[i] = 0.0;
    }

    if (key == 2)
    {
      d = d * (extOne - extOne / mne); k = 0;
    }
    else if (key > 2)
    {
      Extended dr = mne / nexp;
      dr = dr - fmod(dr, extOne);
      //dr=(dr>0)?floor(dr):ceil(dr);
      Extended dd = mne - dr;
      dr = d * dd;
      dd = dr - fmod(dr, extOne);
      //dd=(dr>0)?floor(dr):ceil(dr);
      dr = dd + (dd - extOne) / (nexp - extOne);
      dd = dr - fmod(dr, extOne);
      //dd=(dr>0)?floor(dr):ceil(dr);
      d = dd*(extOne / (mne - extOne));
    }

    for (int j = 0; j < m; j++)
    {
      iq = 0;
      if (x == extOne)
      {
        is = nexp - extOne; d = extZero;
      }
      else
      {
        d = d * nexp;
        is = floor(d);
        //is = (int)d.toDouble(); //опасное преобразование при n > 32
        d = d - is;
      }
      node(is, n1, nexp, l, iq, iu, iv, extOne, extZero, extHalf);
      std::swap(iu[0], iu[it]);
      std::swap(iv[0], iv[it]);
      if (l == 0)
        l = it;
      else if (l == it) l = 0;
      if ((iq > 0) || ((iq == 0) && (is == 0))) k = l;
      else if (iq<0) k = (it == n1) ? 0 : n1;
      r = r * 0.5;
      it = l;
      for (int i = 0; i<n; i++)
      {
        iu[i] = iu[i] * iw[i];
        iw[i] = -iv[i] * iw[i];
        p = r * iu[i];
        p = p + y[i];
        y[i] = p;
      }
    }
    if (key == 2)
    {
      int i;
      if (is == (nexp - extOne)) i = -1;
      else i = 1;
      p = 2 * i * iu[k] * r*d.toDouble();
      p = y[k] - p;
      y[k] = p;
    }
    else if (key == 3)
    {
      for (int i = 0; i < n; i++)
      {
        p = r * iu[i];
        p = p + y[i];
        y[i] = p;
      }
    }
  }
}

// ------------------------------------------------------------------------------------------------
TLinearEvolvent::TLinearEvolvent(int _N, int _m) :
  TEvolvent(_N, _m)
{
  nexpExtended = extOne;
  for (int i = 0; i < _N; nexpExtended *= 2, i++);
  mneExtended = extOne;
  for (int i = 0; i < _m; mneExtended *= nexpExtended, i++);
}

// ------------------------------------------------------------------------------------------------
TLinearEvolvent::~TLinearEvolvent()
{}

// ------------------------------------------------------------------------------------------------
void TLinearEvolvent::GetImage(const Extended& x, double* _y, int EvolventNum)
{
  if ((x.toDouble() < 0) || (x.toDouble() > 1))
  {
    throw EXCEPTION("x is out of range");
  }
  // x ---> y
  if (N == 1)
  {
    y[0] = x.toDouble() - 0.5;
    return;
  }
  mapd(x, m, y, N, 2, nexpExtended, mneExtended, extOne, extNull, extHalf);
  transform_P_to_D();
  memcpy(_y, y, N * sizeof(double));
}

// ------------------------------------------------------------------------------------------------
void TLinearEvolvent::GetInverseImage(double* _y, Extended& x)
{
  // y ---> x
  memcpy(y, _y, N * sizeof(double));
  transform_D_to_P();
  if (N == 1)
  {
    x = y[0] + 0.5;
    return;
  }
  xyd(&x, m, y, N, nexpExtended, extOne, extNull, extHalf);
}

// ------------------------------------------------------------------------------------------------
TNoninjectiveEvolvent::TNoninjectiveEvolvent(int _N, int _m, int _max_preimages) :
  TEvolvent(_N, _m)
{
  max_preimages = _max_preimages;
  nexpExtended = extOne;
  for (int i = 0; i < _N; nexpExtended *= 2, i++);
  mneExtended = extOne;
  for (int i = 0; i < _m; mneExtended *= nexpExtended, i++);
}

// ------------------------------------------------------------------------------------------------
TNoninjectiveEvolvent::~TNoninjectiveEvolvent()
{}

// ------------------------------------------------------------------------------------------------
void TNoninjectiveEvolvent::GetImage(const Extended& x, double* _y, int EvolventNum)
{
  if ((x.toDouble() < 0) || (x.toDouble() > 1))
  {
    throw EXCEPTION("x is out of range");
  }
  // x ---> y
  if (N == 1)
  {
    y[0] = x.toDouble() - 0.5;
    return;
  }
  mapd(x, m, y, N, 3, nexpExtended, mneExtended, extOne, extNull, extHalf);
  transform_P_to_D();
  memcpy(_y, y, N * sizeof(double));
}

// ------------------------------------------------------------------------------------------------
int TNoninjectiveEvolvent::GetNoninjectivePreimages(double* _y, Extended* x)
{
  memcpy(y, _y, N * sizeof(double));
  transform_D_to_P();
  int preimNumber = 1;
  if (N == 1)
    x[0] = y[0] + 0.5;
  else
    invmad(m, x, max_preimages, &preimNumber, y, N, 4, nexpExtended, mneExtended, extOne, extNull, extHalf);
  return preimNumber;
}
// - end of file ----------------------------------------------------------------------------------
