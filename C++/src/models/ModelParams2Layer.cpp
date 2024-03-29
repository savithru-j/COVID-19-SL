#include "ModelParams2Layer.h"
#include "linearalgebra/SurrealS.h"
#include "utils/ErrorHandler.h"

template<class T>
constexpr double ModelParams2Layer<T>::dt;

template<class T>
void copyParam2FullVector(const ModelParams2Layer<T>& params, Vector<T>& v)
{
  const int nt = params.nt_hist + params.nt_pred;
  const int m = 2*nt + 6;
  if (v.m() != m)
    throwError("copyParam2FullVector - inconsistent dimensions!\n"
               + std::to_string(v.m()) + " != " + std::to_string(m));

  for (int i = 0; i < nt; ++i)
  {
    v[     i] = params.beta[i];
    v[nt + i] = params.IFR[i];
  }
  const int off = 2*nt;
  v[off  ] = params.T_incub;
  v[off+1] = params.T_recov;
  v[off+2] = params.beta_test_scaling;
  v[off+3] = params.beta_vac_scaling;
  v[off+4] = params.vaccine_alpha;
  v[off+5] = params.IFR_vac_scaling;
}

template<class T>
void copyFullVector2Param(const Vector<T>& v, ModelParams2Layer<T>& params)
{
  const int nt = params.nt_hist;
  const int m = 2*nt + 6;
  if (v.m() != m)
    throwError("copyFullVector2Param - inconsistent dimensions!");

  for (int i = 0; i < nt; ++i)
  {
    params.beta[i] = v[     i];
    params.IFR[i]  = v[nt + i];
  }

  const int N = params.nt_hist + params.nt_pred;
  for (int i = nt; i < N; ++i) //Copy last "history" value to prediction section
  {
    params.beta[i] = params.beta[nt-1];
    params.IFR[i]  = params.IFR[nt-1];
  }

  const int off = 2*nt;
  params.T_incub           = v[off  ];
  params.T_recov           = v[off+1];
  params.beta_test_scaling = v[off+2];
  params.beta_vac_scaling  = v[off+3];
  params.vaccine_alpha     = v[off+4];
  params.IFR_vac_scaling   = v[off+5];
}

//Explicit instantiations
template class ModelParams2Layer<double>;
template class ModelParams2Layer<SurrealS<1,double>>;

template void copyParam2FullVector(const ModelParams2Layer<double>& params, Vector<double>& v);
template void copyFullVector2Param(const Vector<double>& v, ModelParams2Layer<double>& params);
