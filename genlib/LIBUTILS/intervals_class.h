
#include <map>
#include <deque>


   ///////////////////////////////////////////
   // MODULE TO HANDLE A TIME ORDERED CHAIN //
   ///////////////////////////////////////////

using namespace std;


typedef multimap<double,int> Tinterval;


class IntervalIndex
{
public: 
  int ifl;
  int type;
  int in;
  IntervalIndex(int ifl_, int type_, int in_) : ifl(ifl_), type(type_), in(in_){};
  IntervalIndex(int x=-1) : ifl(x), type(x), in(x){}; 

  // This is for conversion from integer to default state of the class
  // This state of the index should never be used. The last slot in Operators (at t=beta) does not contain this index. This
  // default value -1 will be set in this case.
};



inline std::ostream& operator<< (std::ostream& stream, const IntervalIndex& r)
{
  int width = stream.width();
  stream << std::setw(width)<< r.ifl << " " << std::setw(width) << r.type << " " << std::setw(width) << r.in;
  return stream;
}



class nIntervals{
  double beta;
  function1D<double> _time_e, _time_s;

  function2D<dcomplex> _exp_e, _exp_s;

  function2D<dcomplex> b_exp_e, b_exp_s; 
  // bosonic frequencies in case needed

  function1D<int> index_e, index_s;
  function1D<int> _btype_e, _btype_s;
  deque<int> Empty_e, Empty_s;

  const mesh1D* piom;  
   // fermionic frequencies

  const mesh1D* biom;  
   // bosonic frequencies

  function1D<int> Nbtyp_s, Nbtyp_e;

public: 
  // members
  static const int cd = 0;
  static const int c = 1;

public:  
  // methods  
  nIntervals() {};
  void SetUp(int N_max, const mesh1D& iom, const mesh1D& iomb, double beta_, int dim);


  // inlined short methods to access class members
  int size() const {return _time_e.size()+_time_s.size();}
  int fullsize() const {return _time_e.fullsize()+_time_s.fullsize();}
  const funProxy<dcomplex>& exp_e(int i) const {return _exp_e[index_e[i]];}
  const funProxy<dcomplex>& exp_s(int i) const {return _exp_s[index_s[i]];}
  const funProxy<dcomplex>& exp(int type, int i) const {return type==cd ? _exp_s[index_s[i]] : _exp_e[index_e[i]];}
  double time_s(int i) const {return _time_s[index_s[i]];}
  double time_e(int i) const {return _time_e[index_e[i]];}
  double time(int type, int i) const {return type==cd ? _time_s[index_s[i]] : _time_e[index_e[i]];}
  
  double time_direct(int type, int i) const {return type==cd ? _time_s[i] : _time_e[i];}
  
  template <int boson_fermion>
  const funProxy<dcomplex>& exp_direct(int type, int i) const {cerr<<"Should not happen"<<endl; return NULL;}
  
  int btype_s(int i) const {return _btype_s[index_s[i]];}
  int btype_e(int i) const {return _btype_e[index_e[i]];}
  int btype(int type, int i) const {return type==cd ? _btype_s[index_s[i]] : _btype_e[index_e[i]];}
  int Nbtype(int type, int b) const {return type==cd ? Nbtyp_s[b] : Nbtyp_e[b];}
  void Find_is_ie(double t_start, int& is, double t_end, int& ie);
  pair<int,int> InsertExponents(double t_start, int is, int btyp_s, double t_end, int ie, int btyp_e);
  void RemoveExponents(int is, int ie);
  void MoveExponents(int type, double t_old, int i_old, double t_new, int i_new);
  int FindIndex(double time, double t_old, int type);
  void print(ostream& stream);
  void TestOrder();
  int FindSuccessive(int type, int ie, int bfle);
  double PreviousTime(int type, int to_move);
  double NextTime(int type, int to_move);
  void copy_data(const nIntervals& source);
};


template <>
const funProxy<dcomplex>& nIntervals::exp_direct<0>(int type, int i) const
{return type==cd ? b_exp_s[i] : b_exp_e[i];}

template <>
const funProxy<dcomplex>& nIntervals::exp_direct<1>(int type, int i)  const
{return type==cd ? _exp_s[i] : _exp_e[i];}

 


 

double nIntervals::PreviousTime(int type, int to_move)
{
  if (type==c){
    int bfl = _btype_e[index_e[to_move]];
    for (int ii=to_move-1; ii>=0; ii--)
      if (_btype_e[index_e[ii]]==bfl) return _time_e[index_e[ii]];
    for (int ii=_time_e.size()-1; ii>=to_move; ii--)
      if (_btype_e[index_e[ii]]==bfl) return _time_e[index_e[ii]];
  }else{
    int bfl = _btype_s[index_s[to_move]];
    for (int ii=to_move-1; ii>=0; ii--)
      if (_btype_s[index_s[ii]]==bfl) return _time_s[index_s[ii]];
    for (int ii=_time_s.size()-1; ii>=to_move; ii--)
      if (_btype_s[index_s[ii]]==bfl) return _time_s[index_s[ii]];
  }
  return 0;
}



double nIntervals::NextTime(int type, int to_move)
{
  if (type==c){
    int bfl = _btype_e[index_e[to_move]];
    for (int ii=to_move+1; ii<_time_e.size(); ii++)
      if (_btype_e[index_e[ii]]==bfl) return _time_e[index_e[ii]];
    for (int ii=0; ii<=to_move; ii++)
      if (_btype_e[index_e[ii]]==bfl) return _time_e[index_e[ii]];
  }else{
    int bfl = _btype_s[index_s[to_move]];
    for (int ii=to_move+1; ii<_time_s.size(); ii++)
      if (_btype_s[index_s[ii]]==bfl) return _time_s[index_s[ii]];
    for (int ii=0; ii<=to_move; ii++)
      if (_btype_s[index_s[ii]]==bfl) return _time_s[index_s[ii]];
  }
  return 0;
}



int nIntervals::FindSuccessive(int type, int ix, int bflx)
{
  int ii=0;
  if (type==c){
    for (int i=0; i<ix; i++) if (_btype_e[index_e[i]]==bflx) ii++;
  } else{
    for (int i=0; i<ix; i++) if (_btype_s[index_s[i]]==bflx) ii++;
  }
  return ii;
}



void nIntervals::SetUp(int N_max, const mesh1D& iom, const mesh1D& iomb, double beta_, int dim=1)
{
  beta = beta_;
  int nom = iom.size();
  int nomb = iomb.size();
  _time_e.resize(N_max); _time_e.resize_virtual(0);
  _time_s.resize(N_max); _time_s.resize_virtual(0);
  _exp_e.resize(N_max,nom);
  _exp_s.resize(N_max,nom);
  b_exp_e.resize(N_max,nomb);
  b_exp_s.resize(N_max,nomb);
  index_e.resize(N_max);
  index_s.resize(N_max);
  Empty_e.resize(N_max);
  Empty_s.resize(N_max);
  _btype_s.resize(N_max);
  _btype_e.resize(N_max);
  Nbtyp_s.resize(dim);
  Nbtyp_s=0;
  Nbtyp_e.resize(dim);
  Nbtyp_e=0;
  piom = &iom;
  biom = &iomb;
  for (int i=0; i<N_max; i++) Empty_e[i]=i;
  for (int i=0; i<N_max; i++) Empty_s[i]=i;
};




inline pair<int,int> nIntervals::InsertExponents(double t_start, int is, int btyp_s, double t_end, int ie, int btyp_e)
{
  int nsize = index_e.size()-Empty_e.size();
  
  int to_insert_c = Empty_e.front();
  Empty_e.pop_front();
  for (int i=nsize; i>ie; i--) index_e[i] = index_e[i-1];
  index_e[ie] = to_insert_c;
  
  int to_insert_cd = Empty_s.front();
  Empty_s.pop_front();
  for (int i=nsize; i>is; i--) index_s[i] = index_s[i-1];    
  index_s[is] = to_insert_cd;
  
  _time_e[to_insert_c] = t_end;
  _time_s[to_insert_cd] = t_start;

  _time_e.resize_virtual(_time_e.size()+1);
  _time_s.resize_virtual(_time_s.size()+1);
  
  _btype_s[to_insert_cd] = btyp_s;
  _btype_e[to_insert_c ] = btyp_e;
  Nbtyp_e[btyp_e]++;
  Nbtyp_s[btyp_s]++;
    
  for (int i=0; i<piom->size(); i++){
    double xe = t_end*(*piom)[i];
    _exp_e[to_insert_c][i] = dcomplex(cos(xe),sin(xe));
    double xs = t_start*(*piom)[i];
    _exp_s[to_insert_cd][i] = dcomplex(cos(xs),sin(xs));
  }

  //////////////////////////////////////////////////////////
  //if (common::SampleSusc){
    for (int i=0; i<biom->size(); i++){
      double xe = t_end*(*biom)[i];
      b_exp_e[to_insert_c][i] = dcomplex(cos(xe),sin(xe));
      double xs = t_start*(*biom)[i];
      b_exp_s[to_insert_cd][i] = dcomplex(cos(xs),sin(xs));
    }
  //}
  //////////////////////////////////////////////////////////

  return make_pair(to_insert_cd, to_insert_c);
}



inline void nIntervals::RemoveExponents(int is, int ie)
{
  int nsize = index_e.size()-Empty_e.size();

  Nbtyp_e[_btype_e[ie]]--;
  Nbtyp_s[_btype_s[is]]--;
  
  Empty_e.push_back(index_e[ie]);
  Empty_s.push_back(index_s[is]);
  for (int i=is; i<nsize-1; i++) index_s[i] = index_s[i+1];
  for (int i=ie; i<nsize-1; i++) index_e[i] = index_e[i+1];
  
  _time_e.resize_virtual(_time_e.size()-1);
  _time_s.resize_virtual(_time_s.size()-1);
}





inline void nIntervals::MoveExponents(int type, double t_old, int i_old, double t_new, int i_new)
{
  int nsize = index_e.size()-Empty_e.size();
  
  if (type==c){
    int to_insert = index_e[i_old];
    _time_e[to_insert] = t_new;
    for (int i=i_old; i<nsize-1; i++) index_e[i] = index_e[i+1];
    for (int i=nsize-1; i>i_new; i--) index_e[i] = index_e[i-1];
    index_e[i_new] = to_insert;
    for (int i=0; i<piom->size(); i++){
      double xe = t_new*(*piom)[i];
      _exp_e[to_insert][i] = dcomplex(cos(xe),sin(xe));
    }


   //////////////////////////////////////////////////////
   //    if (common::SampleSusc){
      for (int i=0; i<biom->size(); i++){
        double xe = t_new*(*biom)[i];
        b_exp_e[to_insert][i].Set(cos(xe),sin(xe));
      }
   //    }
   //////////////////////////////////////////////////////


  }else{
    int to_insert = index_s[i_old];
    _time_s[to_insert] = t_new;
    for (int i=i_old; i<nsize-1; i++) index_s[i] = index_s[i+1];
    for (int i=nsize-1; i>i_new; i--) index_s[i] = index_s[i-1];
    index_s[i_new] = to_insert;
    for (int i=0; i<piom->size(); i++){
      double xe = t_new*(*piom)[i];
      _exp_s[to_insert][i] = dcomplex(cos(xe),sin(xe));
    }
    for (int i=0; i<biom->size(); i++){
      double xe = t_new*(*biom)[i];
      b_exp_s[to_insert][i].Set(cos(xe),sin(xe));
    }
  }
}




inline int nIntervals::FindIndex(double t_new, double t_old, int type)
{
  int is=0;
  if (type==cd){
    while (is<_time_s.size() && time_s(is)<t_new) is++;
  }else{
    while (is<_time_e.size() && time_e(is)<t_new) is++;
  }

  if (t_new>t_old) is--;  
  return is;
}




inline void nIntervals::Find_is_ie(double t_start, int& is, double t_end, int& ie)
{
  is=0;
  while (is<_time_s.size() && time_s(is)<t_start) is++;
  ie=0;
  while (ie<_time_e.size() && time_e(ie)<t_end) ie++;
}


  
inline void nIntervals::TestOrder()
{
  int nsize = index_e.size()-Empty_e.size();
  for (int i=0; i<nsize-1; i++)
    if (time_s(i)>time_s(i+1)) cout<<"Times are not ordered! "<<time_s(i)<<" "<<time_s(i+1)<<" "<<i<<" "<<i+1<<endl;
  for (int i=0; i<nsize-1; i++)
    if (time_e(i)>time_e(i+1)) cout<<"Times are not ordered! "<<time_e(i)<<" "<<time_e(i+1)<<" "<<i<<" "<<i+1<<endl;
}



inline void nIntervals::print(ostream& stream)
{
  stream<<setw(4)<<_time_s.size()<<":  s= ";
  for (int is=0; is<_time_s.size(); is++) stream<<setw(4)<<time_s(is)<<" ";
  cout<<" e= ";
  for (int ie=0; ie<_time_e.size(); ie++) stream<<setw(4)<<time_e(ie)<<" ";
  stream<<endl;
}



void nIntervals::copy_data(const nIntervals& s)
{
  //  beta= s.beta;
  _time_e.copy_full(s._time_e);
  _time_s.copy_full(s._time_s);
  _exp_e = s._exp_e;
  _exp_s = s._exp_s;
  b_exp_e = s.b_exp_e;
  b_exp_s = s.b_exp_s;
  index_e = s.index_e;
  index_s = s.index_s;
  _btype_e.copy_full(s._btype_e);
  _btype_s.copy_full(s._btype_s);
  Empty_e = s.Empty_e;
  Empty_s = s.Empty_s;
  Nbtyp_s.copy_full(s.Nbtyp_s);
  Nbtyp_e.copy_full(s.Nbtyp_e);
}




