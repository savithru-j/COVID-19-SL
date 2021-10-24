#include "ModelParams.h"
#include "SurrealS.h"

template<class T>
constexpr double ModelParams<T>::dt;

template<class T>
void copyParam2FullVector(const ModelParams<T>& params, Vector<T>& v)
{
  const int nt = params.nt_hist + params.nt_pred;
  const int m = 6*nt + 10;
  if (v.m() != m)
    throwError("copyParam2FullVector - inconsistent dimensions!\n"
               + std::to_string(v.m()) + " != " + std::to_string(m));

  for (int i = 0; i < nt; ++i)
  {
    v[       i] = params.betaN[i];
    v[  nt + i] = params.c0[i];
    v[2*nt + i] = params.c1[i];
    v[3*nt + i] = params.c2[i];
    v[4*nt + i] = params.c3[i];
    v[5*nt + i] = params.IFR[i];
  }
  const int off = 6*nt;
  v[off  ] = params.T_incub0;
  v[off+1] = params.T_incub1;
  v[off+2] = params.T_asympt;
  v[off+3] = params.T_mild;
  v[off+4] = params.T_severe;
  v[off+5] = params.T_icu;
  v[off+6] = params.f;
  v[off+7] = params.frac_recover_I1;
  v[off+8] = params.frac_recover_I2;
  v[off+9] = params.vaccine_eff;
}

template<class T>
void copyFullVector2Param(const Vector<T>& v, ModelParams<T>& params)
{
  const int nt = params.nt_hist;
  const int m = 6*nt + 10;
  if (v.m() != m)
    throwError("copyFullVector2Param - inconsistent dimensions!");

  for (int i = 0; i < nt; ++i)
  {
    params.betaN[i] = v[       i];
    params.ce[i]    = v[  nt + i];
    params.c0[i]    = v[  nt + i]; //ce = c0
    params.c1[i]    = v[2*nt + i];
    params.c2[i]    = v[3*nt + i];
    params.c3[i]    = v[4*nt + i];
    params.IFR[i]   = v[5*nt + i];
  }

  const int N = params.nt_hist + params.nt_pred;
  for (int i = nt; i < N; ++i) //Copy last "history" value to prediction section
  {
    params.betaN[i] = params.betaN[nt-1];
    params.ce[i]    = params.ce[nt-1];
    params.c0[i]    = params.c0[nt-1];
    params.c1[i]    = params.c1[nt-1];
    params.c2[i]    = params.c2[nt-1];
    params.c3[i]    = params.c3[nt-1];
    params.IFR[i]   = params.IFR[nt-1];
  }

  const int off = 6*nt;
  params.T_incub0        = v[off  ];
  params.T_incub1        = v[off+1];
  params.T_asympt        = v[off+2];
  params.T_mild          = v[off+3];
  params.T_severe        = v[off+4];
  params.T_icu           = v[off+5];
  params.f               = v[off+6];
  params.frac_recover_I1 = v[off+7];
  params.frac_recover_I2 = v[off+8];
  params.vaccine_eff     = v[off+9];
}

//Explicit instantiations
template class ModelParams<double>;
template class ModelParams<SurrealS<1,double>>;

template void copyParam2FullVector(const ModelParams<double>& params, Vector<double>& v);
template void copyFullVector2Param(const Vector<double>& v, ModelParams<double>& params);
