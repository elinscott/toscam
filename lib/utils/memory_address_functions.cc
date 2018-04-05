#include <iostream>
using namespace std;

#ifdef LOGICAL1   
#define bool int
#endif

extern "C" int same_address_(double &i1,double &i2)
{
if(&i1 == &i2){
return 1;
}else{
return 0;
};
}

extern "C" int same_address_r_(float &i1,float &i2)
{
if(&i1 == &i2){
return 1;
}else{
return 0;
};
}


extern "C" long int getaddress_char_(char &i1)
{
return (long int) &i1;
}


extern "C" long int getaddress_r_(float &i1)
{
return (long int) &i1;
}


extern "C" long int getaddress_(double &i1)
{
return (long int) &i1;
}

extern "C" long int getaddress_vec_(double i1[])
{
return (long int) &i1[0];
}


extern "C" long int getaddress_i_(int &i1)
{
return (long int) &i1;
}

extern "C" long int getaddress_l_(bool &i1)
{
return (long int) &i1;
}

extern "C" int putaddress_char_(long int &j, char &i1)
{
 char *p;
  p = reinterpret_cast<char*>(j);
 *p=i1;
 return 0;
}


extern "C" int putaddress_r_(long int &j, float &i1)
{

 float *p;
  p = reinterpret_cast<float*>(j);
 *p=i1;
 return 0;
}


extern "C" int putaddress_(long int &j, double &i1)
{
 double *p;
  p = reinterpret_cast<double*>(j);
 *p=i1;
 return 0;
}


extern "C" int putaddress_i_(long int &j, int &i1)
{
  int *p;
  p = reinterpret_cast<int*>(j);
 *p=i1;
 return 0;
}

extern "C" int putaddress_l_(long int &j, bool &i1)
{
  bool *p;
  p = reinterpret_cast<bool*>(j);
 *p=i1;
 return 0;
}


