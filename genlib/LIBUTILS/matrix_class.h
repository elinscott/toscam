
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/************************* MATRIX CLASS ****************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/


 ////////////////////////////////////////////////////
 // MATRIX UPDATES FOR CONTINUOUS TIME MONTE CARLO //
 ////////////////////////////////////////////////////


 // more generally, updates of determinant and inverse of a matrix
 // when columns are added, removed, or changed....


class MatrixM{

  double beta;
  const mesh1D* ptau;
  function2D<spline1D<double>* > vfun;
  const mesh1D* piom;
  function1D<double> d_ln, d_nl, d_lm, d_ml;
  double d_nn, ratio;
  function1D<double> Lt, Rt, Lt2, Rt2;
  function1D<double> temp;
  function2D<double> C, mD;
  function1D<int> ipiv;
  function2D<dcomplex> Gf;
  function2D<dcomplex> sum1, sum2, sum3;
  function1D<dcomplex> sum0;
  function2D<dcomplex> dG;
  const function2D<int>* ptfl_index;

public:
  void copy_data(MatrixM& source);

     // small inlined functions to access members

  int fullsize(){return Lt.fullsize();}
  const function2D<dcomplex>& gf() const { return Gf;}
  double iom_size(){return piom->size();}
  
     //  void SetUp(int N_max, const mesh1D& tau_, const vector<vector<spline1D<double> > >& fun_, const function2D<int>& tfl_index_, const mesh1D& iom, double beta_);

  void SetUp(int N_max, const mesh1D& tau_, const function2D<spline1D<double>* >& fun_, const function2D<int>& tfl_index_, const mesh1D& iom, double beta_);
  
  double AddDetRatio(function2D<double>& MD, int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval);
  double AddUpdateMatrix(function2D<double>& MD, int Nk, double t_start, int is, int btype_s, double t_end, int ie, int btype_e, const nIntervals& interval);
  void AddUpdate_Gf(const nIntervals& interval);
  
  double RemoveDetRatio(function2D<double>& MD, int Nk, int is, int ie);
  double RemoveDetRatio(function2D<double>& MD, int Nk, int is1, int ie1, int is2, int ie2);
  double RemoveUpdateMatrix(function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval);
  void RemoveUpdate_Gf(const function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval);
  
  double Move_start_DetRatio(function2D<double>& MD, int Nk, double ts_old, double ts_new, int btype_s, int is_old, const nIntervals& interval);
  double Move_end_DetRatio(function2D<double>& MD, int Nk, double te_old, double te_new, int btype_e, int ie_old, const nIntervals& interval);
  void Move_start_UpdateMatrix(function2D<double>& MD, int Nk, double ts_old, double ts_new, int btype_s, int is_old, int is_new, const nIntervals& interval);
  void Move_end_UpdateMatrix(function2D<double>& MD, int Nk, double te_old, double te_new, int btype_e, int ie_old, int ie_new, const nIntervals& interval);
  void MoveUpdate_Gf(const nIntervals& interval);
  
  
  void CleanUpdateMatrix(function2D<double>& MD, int Nk, const nIntervals& interval);
  void CleanUpdateGf(function2D<double>& MD, int Nk, const nIntervals& interval);
  

  double AddDetRatio(function2D<double>& MD, int Nk, int btype_s1, double t_start1, int btype_e1, double t_end1, 
                     int btype_s2, double t_start2, int btype_e2, double t_end2, const nIntervals& interval);
  double GlobalFlipDetRatio(const nIntervals& interval, const function2D<double>& MD, const function2D<spline1D<double>* >& Delta_flipped);

  void ComputeGtau(function2D<double>& MD, const nIntervals& interval, function2D<double>& Gtau);
private:
  void AddComputeGs(int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval);
  void AddComputeGs(int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval, function1D<double>& d_ln, function1D<double>& d_nl, double& dnn);
  friend void swapGf(MatrixM& m1, MatrixM& m2);
  double Interp_Delta(int bs, int be, const intpar& p);
  double Interp_Antiperiodic_Delta(int bs, int be, double dt);
  void SetMatrix_(int Nk, const nIntervals& interval, function2D<double>& mD, const function2D<spline1D<double>* >& Delta_any);
};


/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/


 //inline void MatrixM::SetUp(int N_max, const mesh1D& tau_, const vector<vector<spline1D<double> > >& fun_, const function2D<int>& tfl_index_, const mesh1D& iom, double beta_)
inline void MatrixM::SetUp(int N_max, const mesh1D& tau_, const function2D<spline1D<double>* >& fun_, const function2D<int>& tfl_index_, const mesh1D& iom, double beta_)
{
  beta = beta_;
  ptau = &tau_;
  ptfl_index = &tfl_index_;
  vfun = fun_;
  d_ln.resize(N_max);
  d_nl.resize(N_max);
  d_lm.resize(N_max);
  d_ml.resize(N_max);
  Lt.resize(N_max);
  Rt.resize(N_max);
  Lt2.resize(N_max);
  Rt2.resize(N_max);
  C.resize(N_max,N_max);
  mD.resize(N_max,N_max);
  ipiv.resize(N_max);
  piom = &iom;
  int dim = ptfl_index->size_N();
  Gf.resize(dim*dim,iom.size());
  sum1.resize(dim,iom.size());
  sum2.resize(dim,iom.size());
  sum3.resize(dim*dim,iom.size());
  sum0.resize(dim*dim);
  temp.resize(N_max);
  dG.resize(dim*dim,iom.size());
}




inline double Interp_Any_Delta(const spline1D<double>& Delta_any, const intpar& p, bool bs_equal_be)
{
  // Interpolates Delta(tau) and makes sure that Delta is causal
  double Delta_tau = Delta_any(p); // interpolation of Delta(tau) in point p.x==tau

 ///////////////////////////////////////////////////////////////////////
 //  if (bs_equal_be && Delta_tau>-common::minDeltat)
 //   Delta_tau = -common::minDeltat; // Fix any causality problem
 ///////////////////////////////////////////////////////////////////////
 
   if (bs_equal_be && Delta_tau>-0.0000001)
    Delta_tau = -0.0000001; // Fix any causality problem



  return Delta_tau;
}




inline double MatrixM::Interp_Delta(int bs, int be, const intpar& p)
{ // Interpolates Delta stored in MatrixM class. This Interpolation is for speed only.
  // It could be enough to have only Interp_Antiperiodic_Delta. In the latter case, 
  // program would not be optimized for correlated lookups.
  return Interp_Any_Delta(*vfun[bs][be], p, bs==be);
}



inline double MatrixM::Interp_Antiperiodic_Delta(int bs, int be, double dt)
{  // Interpolates Delta antiperiodically (because Delta is fermionic)
  if (dt>0)
    return Interp_Delta(bs, be, ptau->Interp(dt));
  else
    return -Interp_Delta(bs, be, ptau->Interp(dt+beta));
}



inline void MatrixM::AddComputeGs(int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval)
{
  tint pos1 = ptau->InitInterpLeft();
  tint pos2 = ptau->InitInterpLeft();
  for (int i=0; i<Nk; i++){
    int bs = interval.btype_s(i);
    int be = btype_e;
    if (interval.time_s(i)<t_end)
      d_ln[i] = -Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-t_end+beta,pos1));
    else
      d_ln[i] = Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-t_end,pos2));
  }
  pos1 = ptau->InitInterpRight();
  pos2 = ptau->InitInterpRight();
  for (int i=0; i<Nk; i++){
    int bs = btype_s; 
    int be = interval.btype_e(i);
    if (t_start<interval.time_e(i))
      d_nl[i] = -Interp_Delta(bs, be, ptau->InterpRight(t_start-interval.time_e(i)+beta,pos1));
    else
      d_nl[i] = Interp_Delta(bs, be, ptau->InterpRight(t_start-interval.time_e(i),pos2));
  }
  d_nn = Interp_Antiperiodic_Delta(btype_s, btype_e, t_start-t_end);
}





inline double MatrixM::AddDetRatio(function2D<double>& MD, int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval)
{
  AddComputeGs(Nk, t_start, btype_s, t_end, btype_e, interval);
  for (int i=0; i<Nk; i++){
    double sumL=0;
    double sumR=0;
    for (int l=0; l<Nk; l++){// should be optimised by dgemv
      sumL += MD(i,l)*d_ln[l];
      sumR += d_nl[l]*MD(l,i);
    }
    Lt[i] = sumL;
    Rt[i] = sumR;
  }
  double qsum=0;
  for (int l=0; l<Nk; l++) qsum += d_nl[l]*Lt[l]; // should be optimized by ddot
  ratio = d_nn - qsum;
  return ratio;
}




inline double MatrixM::Move_start_DetRatio(function2D<double>& MD, int Nk, double ts_old, double ts_new, int btype_s, int is_old, const nIntervals& interval)
{
  tint pos1 = ptau->InitInterpRight();
  tint pos2 = ptau->InitInterpRight();
  for (int i=0; i<Nk; i++){
    int bs = btype_s;
    int be = interval.btype_e(i);
    if (ts_old<interval.time_e(i))
      d_nl[i] = -Interp_Delta(bs, be, ptau->InterpRight(ts_old-interval.time_e(i)+beta,pos1));
    else
      d_nl[i] = Interp_Delta(bs, be, ptau->InterpRight(ts_old-interval.time_e(i),pos2));
  }
  pos1 = ptau->InitInterpRight();
  pos2 = ptau->InitInterpRight();
  for (int i=0; i<Nk; i++){
    int bs = btype_s;
    int be = interval.btype_e(i);
    if (ts_new<interval.time_e(i))
      d_ln[i] = -Interp_Delta(bs, be, ptau->InterpRight(ts_new-interval.time_e(i)+beta,pos1));
    else
      d_ln[i] = Interp_Delta(bs, be, ptau->InterpRight(ts_new-interval.time_e(i),pos2));
  }
  for (int i=0; i<Nk; i++) d_nl[i] = d_ln[i]-d_nl[i]; // new-old
  
  double sumR=0;
  for (int l=0; l<Nk; l++) sumR += d_nl[l]*MD(l,is_old);
  return 1. + sumR;
}




inline double MatrixM::Move_end_DetRatio(function2D<double>& MD, int Nk, double te_old, double te_new, int btype_e, int ie_old, const nIntervals& interval)
{
  tint pos1 = ptau->InitInterpLeft();
  tint pos2 = ptau->InitInterpLeft();
  for (int i=0; i<Nk; i++){
    int bs = interval.btype_s(i);
    int be = btype_e;
    if (interval.time_s(i)<te_old)
      d_ln[i] = -Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-te_old+beta,pos1));
    else
      d_ln[i] = Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-te_old,pos2));
  }
  pos1 = ptau->InitInterpLeft();
  pos2 = ptau->InitInterpLeft();
  for (int i=0; i<Nk; i++){
    int bs = interval.btype_s(i);
    int be = btype_e;
    if (interval.time_s(i)<te_new)
      d_nl[i] = -Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-te_new+beta,pos1));
    else
      d_nl[i] =  Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-te_new,pos2));
  }
  for (int i=0; i<Nk; i++) d_ln[i] = d_nl[i]-d_ln[i]; // new-old
  
  double sumL=0;
  for (int l=0; l<Nk; l++) sumL += MD(ie_old,l)*d_ln[l];
  return 1. + sumL;
}




inline void MatrixM::Move_start_UpdateMatrix(function2D<double>& MD, int Nk, double ts_old, double ts_new, int btype_s, int is_old, int is_new, const nIntervals& interval)
{
  double q = Move_start_DetRatio(MD, Nk, ts_old, ts_new, btype_s, is_old, interval);


  // Correction to the Green's function will be needed
  for (int im=0; im<piom->size(); im++){
    double xs = ts_new*(*piom)[im];
    int bs = btype_s;
    dcomplex dexp = dcomplex(cos(xs),-sin(xs))-interval.exp_s(is_old)[im].conj();
    sum0=0;
    for (int i=0; i<Nk; i++){
      int be = interval.btype_e(i);
      int ind = (*ptfl_index)[be][bs];
      sum0[ind] += interval.exp_e(i)[im]*MD(i,is_old)*dexp;
    }
    for (int ind=0; ind<sum0.size(); ind++)
      for (int i=0; i<sum0.size(); i++) dG[ind][im] = -sum0[ind]/beta;
  }
  
    
  Lt.resize(Nk); Rt.resize(Nk);
  for (int i=0; i<Nk; i++){
    double sumR=0;
    for (int l=0; l<Nk; l++){// should be optimised by dgemv
      sumR += d_nl[l]*MD(l,i);
    }
    Rt[i] = sumR/q;
  }
  for (int i=0; i<Nk; i++) Lt[i] = -MD(i,is_old);
  for (int i=0; i<Nk; i++)
    for (int j=0; j<Nk; j++)
      MD(i,j) += Lt[i]*Rt[j];

  
  if (is_new!=is_old){
    double rt=0;
    for (int i=0; i<Nk; i++){ // remember the old row
      temp[i] = MD(i,is_old);
      rt = Rt[is_old];
    }
    if (is_new<is_old){// shuffle rows
      for (int j=is_old; j>is_new; j--){
        for (int i=0; i<Nk; i++) MD(i,j) = MD(i,j-1);
        Rt[j] = Rt[j-1];
      }
      for (int i=0; i<Nk; i++) MD(i,is_new) = temp[i];
      Rt[is_new] = rt;
    }else{
      for (int j=is_old; j<is_new; j++){
        for (int i=0; i<Nk; i++) MD(i,j) = MD(i,j+1);
        Rt[j] = Rt[j+1];
      }
      for (int i=0; i<Nk; i++) MD(i,is_new) = temp[i];
      Rt[is_new] = rt;
    }
  }
}




inline void MatrixM::Move_end_UpdateMatrix(function2D<double>& MD, int Nk, double te_old, double te_new, int btype_e, int ie_old, int ie_new, const nIntervals& interval)
{
  double q = Move_end_DetRatio(MD, Nk, te_old, te_new, btype_e, ie_old, interval);


  // Correction to the Green's function will be needed
  for (int im=0; im<piom->size(); im++){
    double xe = te_new*(*piom)[im];
    int be = btype_e;
    dcomplex dexp = dcomplex(cos(xe),sin(xe))-interval.exp_e(ie_old)[im];
    sum0=0;
    for (int i=0; i<Nk; i++){
      int bs = interval.btype_s(i);
      int ind = (*ptfl_index)[be][bs];
      sum0[ind] += dexp*MD(ie_old,i)*interval.exp_s(i)[im].conj();
    }
    for (int ind=0; ind<sum0.size(); ind++)
      for (int i=0; i<sum0.size(); i++) dG[ind][im] = -sum0[ind]/beta;
  }
  
  Lt.resize(Nk); Rt.resize(Nk);
  for (int i=0; i<Nk; i++){
    double sumL=0;
    for (int l=0; l<Nk; l++){// should be optimised by dgemv
      sumL += MD(i,l)*d_ln[l];
    }
    Lt[i] = sumL/q;
  }
  for (int i=0; i<Nk; i++) Rt[i] = -MD(ie_old,i);
  for (int i=0; i<Nk; i++)
    for (int j=0; j<Nk; j++)
      MD(i,j) += Lt[i]*Rt[j];
  
  if (ie_new!=ie_old){
    double lt = Lt[ie_old];
    for (int i=0; i<Nk; i++) temp[i] = MD(ie_old,i); // remember the old row
    if (ie_new<ie_old){// shuffle columns
      for (int j=ie_old; j>ie_new; j--){
        for (int i=0; i<Nk; i++) MD(j,i) = MD(j-1,i);
        Lt[j] = Lt[j-1];
      }
      for (int i=0; i<Nk; i++) MD(ie_new,i) = temp[i];
      Lt[ie_new] = lt;
    }else{
      for (int j=ie_old; j<ie_new; j++){
        for (int i=0; i<Nk; i++) MD(j,i) = MD(j+1,i);
        Lt[j] = Lt[j+1];
      }
      for (int i=0; i<Nk; i++) MD(ie_new,i) = temp[i];
      Lt[ie_new] = lt;
    }
  }
}




inline double MatrixM::AddUpdateMatrix(function2D<double>& MD, int Nk, double t_start, int is, int btype_s, double t_end, int ie, int btype_e, const nIntervals& interval)
{
  AddDetRatio(MD, Nk, t_start, btype_s, t_end, btype_e, interval);

  Lt.resize(Nk+1); Rt.resize(Nk+1);
  
  for (int i=Nk; i>ie; i--) Lt[i] = Lt[i-1];
  Lt[ie] = -1;
  for (int i=Nk; i>is; i--) Rt[i] = Rt[i-1];
  Rt[is] = -1;
  double v = 1/ratio;
  for (int i=0; i<Nk+1; i++) Lt[i] *= v;
  
  for (int i=Nk; i>ie; i--)
    for (int j=Nk; j>is; j--)
      MD(i,j) = MD(i-1,j-1);
  
  for (int i=0; i<ie; i++)
    for (int j=Nk; j>is; j--)
      MD(i,j) = MD(i,j-1);
  
  for (int i=Nk; i>ie; i--)
    for (int j=0; j<is; j++)
      MD(i,j) = MD(i-1,j);
  
  for (int i=0; i<Nk+1; i++) MD(i,is)=0;
  for (int j=0; j<Nk+1; j++) MD(ie,j)=0;
  
  for (int i=0; i<Nk+1; i++)// should be optimized by dger
    for (int j=0; j<Nk+1; j++)
      MD(i,j) += Lt[i]*Rt[j];

  MD.resize(Nk+1,Nk+1);

  return ratio*(1-2*((is+ie)%2));
}





inline double MatrixM::RemoveDetRatio(function2D<double>& MD, int Nk, int is, int ie)
{
  ratio = MD(ie,is);
  return ratio*(1-2*((is+ie)%2));
}
inline double MatrixM::RemoveDetRatio(function2D<double>& MD, int Nk, int is1, int ie1, int is2, int ie2)
{
  ratio = MD(ie1,is1)*MD(ie2,is2)-MD(ie2,is1)*MD(ie1,is2);
  return ratio*(1-2*((is1+ie1+is2+ie2)%2));
}
inline double MatrixM::RemoveUpdateMatrix(function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval)
{
  RemoveDetRatio(MD,Nk,is,ie);
  
  for (int i=0; i<Nk; i++)
    for (int j=0; j<Nk; j++){
      if (i!=ie && j!=is)
        MD(i,j) -= MD(i,is)*MD(ie,j)/MD(ie,is);
    }

  for (int i=ie; i<Nk-1; i++)
    for (int j=is; j<Nk-1; j++)
      MD(i,j) = MD(i+1,j+1);
  
  for (int i=0; i<ie; i++)
    for (int j=is; j<Nk-1; j++)
      MD(i,j) = MD(i,j+1);
  
  for (int i=ie; i<Nk-1; i++)
    for (int j=0; j<is; j++)
      MD(i,j) = MD(i+1,j);

  MD.resize(Nk-1,Nk-1);
  return ratio*(1-2*((is+ie)%2));
}




inline void MatrixM::AddUpdate_Gf(const nIntervals& interval)
{
  // Updates Green's function
  sum1=0; sum2=0;
  for (int it=0; it<Lt.size(); it++){
    const funProxy<dcomplex>& exp_e = interval.exp_e(it);
    const funProxy<dcomplex>& exp_s = interval.exp_s(it);
    int be = interval.btype_e(it);
    int bs = interval.btype_s(it);
    funProxy<dcomplex>& _sum1 = sum1[be];
    funProxy<dcomplex>& _sum2 = sum2[bs];
    double lt = Lt[it];
    double rt = Rt[it];
    for (int im=0; im<piom->size(); im++){
      _sum1[im] += exp_e[im]*lt;
      _sum2[im] += exp_s[im]*rt;
    }
  }
  for (int be=0; be<ptfl_index->size_N(); be++){
    for (int bs=0; bs<ptfl_index->size_Nd(); bs++){
      int ind = (*ptfl_index)[be][bs];
      for (int im=0; im<piom->size(); im++){
        Gf[ind][im] -= sum1[be][im]*conj(sum2[bs][im])/beta;
      }
    }
  }
}




inline void MatrixM::MoveUpdate_Gf(const nIntervals& interval)
{
  AddUpdate_Gf(interval);
  for (int ind=0; ind<sum0.size(); ind++)
    for (int im=0; im<piom->size(); im++) Gf[ind][im] += dG[ind][im];
}




inline void MatrixM::RemoveUpdate_Gf(const function2D<double>& MD, int Nk, int is, int ie, const nIntervals& interval)
{
  // Updates Green's function
  sum1=0; sum2=0;
  for (int it=0; it<Nk; it++){
    const funProxy<dcomplex>& exp_e = interval.exp_e(it);
    const funProxy<dcomplex>& exp_s = interval.exp_s(it);
    int be = interval.btype_e(it);
    int bs = interval.btype_s(it);
    double me = MD(it,is);
    double ms = MD(ie,it);
    for (int im=0; im<piom->size(); im++){
      sum1[be][im] += exp_e[im]*me;
      sum2[bs][im] += exp_s[im]*ms;
    }
  }
  for (int be=0; be<ptfl_index->size_N(); be++){
    for (int bs=0; bs<ptfl_index->size_Nd(); bs++){
      int ind = (*ptfl_index)[be][bs];  
      for (int im=0; im<piom->size(); im++)
        Gf[ind][im] += sum1[be][im]*conj(sum2[bs][im])/MD(ie,is)/beta;
    }
  }
}



void MatrixM::SetMatrix_(int Nk, const nIntervals& interval, function2D<double>& mD, const function2D<spline1D<double>* >& Delta_any)
{
  mD.resize(Nk,Nk);
  for (int i=0; i<Nk; i++){
    double t_start = interval.time_s(i);
    int bs = interval.btype_s(i);
    tint pos1 = ptau->InitInterpRight(), pos2 = ptau->InitInterpRight();
    for (int j=0; j<Nk; j++){
      double dt = t_start-interval.time_e(j);
      int be = interval.btype_e(j);
      if (dt<0){
        mD(i,j) = -Interp_Any_Delta(*Delta_any[bs][be], ptau->InterpRight(dt+beta,pos1), bs==be);
      }else{
        mD(i,j) =  Interp_Any_Delta(*Delta_any[bs][be], ptau->InterpRight(dt,pos2), bs==be);
      }
    }
  }
}




inline void MatrixM::CleanUpdateMatrix(function2D<double>& MD, int Nk, const nIntervals& interval)
{
  SetMatrix_(Nk, interval, mD, vfun);
  //  SetMatrix_(Nk, interval, mD, vfun, *ptau, beta);
  Inverse(mD, MD, ipiv);
}




inline void MatrixM::CleanUpdateGf(function2D<double>& MD, int Nk, const nIntervals& interval)
{
  sum3=0;
  for (int ie=0; ie<Nk; ie++){
    const funProxy<dcomplex>& exp_e = interval.exp_e(ie);
    int be = interval.btype_e(ie);
    for (int is=0; is<Nk; is++){
      const funProxy<dcomplex>& exp_s = interval.exp_s(is);
      int bs = interval.btype_s(is);      
      double md = MD(ie,is);      
      int ind = (*ptfl_index)[be][bs];
      for (int im=0; im<piom->size(); im++){
        sum3[ind][im] += exp_e[im]*md*conj(exp_s[im]);
      }
    }
  }
  for (int ind=0; ind<sum3.size_N(); ind++)
    for (int im=0; im<piom->size(); im++) Gf[ind][im] = -sum3[ind][im]/beta;
}



inline void MatrixM::AddComputeGs(int Nk, double t_start, int btype_s, double t_end, int btype_e, const nIntervals& interval,
                                  function1D<double>& d_ln, function1D<double>& d_nl, double& d_nn)
{
  tint pos1 = ptau->InitInterpLeft();
  tint pos2 = ptau->InitInterpLeft();
  for (int i=0; i<Nk; i++){
    int bs = interval.btype_s(i);
    int be = btype_e;
    if (interval.time_s(i)<t_end)
      d_ln[i] = -Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-t_end+beta,pos1));
    else
      d_ln[i] = Interp_Delta(bs, be, ptau->InterpLeft(interval.time_s(i)-t_end,pos2));
  }
  pos1 = ptau->InitInterpRight();
  pos2 = ptau->InitInterpRight();
  for (int i=0; i<Nk; i++){
    int bs = btype_s; 
    int be = interval.btype_e(i);
    if (t_start<interval.time_e(i))
      d_nl[i] = -Interp_Delta(bs, be, ptau->InterpRight(t_start-interval.time_e(i)+beta,pos1));
    else
      d_nl[i] = Interp_Delta(bs, be, ptau->InterpRight(t_start-interval.time_e(i),pos2));
  }
  d_nn = Interp_Antiperiodic_Delta(btype_s, btype_e, t_start-t_end);
}



inline double MatrixM::AddDetRatio(function2D<double>& MD, int Nk,
                                   int btype_s1, double t_start1, int btype_e1, double t_end1, 
                                   int btype_s2, double t_start2, int btype_e2, double t_end2, 
                                   const nIntervals& interval)
{
  double d_nn;
  double d_mm;
  AddComputeGs(Nk, t_start1, btype_s1, t_end1, btype_e1, interval, d_ln, d_nl, d_nn);
  AddComputeGs(Nk, t_start2, btype_s2, t_end2, btype_e2, interval, d_lm, d_ml, d_mm);
  double d_nm = Interp_Antiperiodic_Delta(btype_s1, btype_e2, t_start1-t_end2);
  double d_mn = Interp_Antiperiodic_Delta(btype_s2, btype_e1, t_start2-t_end1);
  
  for (int i=0; i<Nk; i++){
    double sumL1=0;
    double sumR1=0;
    double sumL2=0;
    //    double sumR2=0;
    for (int l=0; l<Nk; l++){
      sumL1 += MD(i,l)*d_ln[l];
      sumR1 += d_nl[l]*MD(l,i);
      sumL2 += MD(i,l)*d_lm[l];
      //      sumR2 += d_ml[l]*MD(l,i);
    }
    Lt[i] = sumL1;
    Rt[i] = sumR1;
    Lt2[i] = sumL2;
    //    Rt2[i] = sumR2;
  }
  double sumL_mn=0;
  double sumR_nm=0;
  double sumL_nn=0;
  double sumL_mm=0;
  for (int l=0; l<Nk; l++){
    sumL_mn += d_ml[l]*Lt[l];
    sumR_nm += Rt[l]*d_lm[l];
    sumL_nn += d_nl[l]*Lt[l];
    sumL_mm += d_ml[l]*Lt2[l];
  }
  ratio = (d_nn - sumL_nn)*(d_mm - sumL_mm) - (sumL_mn - d_mn)*(sumR_nm - d_nm);
  return ratio;
}



double MatrixM::GlobalFlipDetRatio(const nIntervals& interval, const function2D<double>& MD, const function2D<spline1D<double>* >& Delta_flipped)
{
  int Nk = interval.size()/2;
  if (Nk<=0) return 1;
  SetMatrix_(Nk, interval, mD, Delta_flipped);
  //  SetMatrix_(Nk, interval, mD, Delta_flipped, *ptau, beta);
  C.Product("N","N",mD,MD);
  return Det(C);
}



void MatrixM::ComputeGtau(function2D<double>& MD, const nIntervals& interval, function2D<double>& Gtau)
{
  // Gtau=0;
  int Nk = interval.size()/2;
  int Gtsize = Gtau[0].size();
  for (int ie=0; ie<Nk; ie++){
    double te = interval.time_e(ie);
    int be = interval.btype_e(ie);
    for (int is=0; is<Nk; is++){
      double ts = interval.time_s(is);
      int bs = interval.btype_s(is);
      double md = MD(ie,is);
      int idt;
      if (te-ts>0){
        idt = static_cast<int>((te-ts)/beta*Gtsize);
      }else{
        idt = static_cast<int>((beta+te-ts)/beta*Gtsize);
        md = -md;
      }
      if (idt>=Gtau[0].size()) idt = Gtau[0].size()-1;
      if (idt<0) idt=0;
      int ind = (*ptfl_index)[be][bs];
      Gtau[ind][idt] -= md;
    }
  }
}




void swapGf(MatrixM& m1, MatrixM& m2)
{
  static function2D<dcomplex> Gtmp;
  Gtmp = m1.Gf;
  m1.Gf = m2.Gf;
  m2.Gf = Gtmp;
}




