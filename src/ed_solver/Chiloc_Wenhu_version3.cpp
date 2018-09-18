#include<complex>
#include<iostream>
#include<fstream>
#include<iomanip>
using namespace std;

// PhiM {{{
/**
 * @brief Fourier transform of  \f$\langle Tcccc \rangle\f$ WITHOUT \f$\exp^{-\beta E_i} \f$ for Matsubara frequencies
 *
 * The number of a term denotes the number of the summand in the following expression (see output of doxygen).
 *\f$
\Phi  \left( E_i, E_j, E_k, E_l; \omega_{1}, \omega_{2}, \omega_{3} \right) =
 \frac{e^{-\beta E_i}}{i \omega_3 + E_{kl}} \left\{
\frac{1-\delta_{\omega_2, -\omega_3}\delta_{E_{jl}, 0}}{i \left( \omega_2 + \omega_3 \right) + E_{jl}}
\left[ \frac{e^{-\beta E_{ji}}+1}{i \omega_1 + E_{ij}} - \frac{e^{-\beta
E_{li} }+1}{i (\omega_1 + \omega_2 + \omega_3) + E_{il}} \right] \right. \\
  \left. +
\delta_{\omega_2, - \omega_3} \delta_{E_{jl}, 0} \left[ \frac{e^{-\beta
E_{ji}}+1}{(i \omega_1 + E_{ij})^2} - \beta \frac{e^{-\beta E_{ji} } }{(i
\omega_1 + E_{ij})} \right] \right. \left. - \frac{1}{i \omega_2 + E_{jk}} \left[
\frac{e^{-\beta E_{ji}}+1}{i\omega_1 + E_{ij}} + \left( 1 -
\delta_{\omega_1, -\omega_2} \delta_{E_{ik}, 0} \right) \frac{e^{-\beta
E_{ki}}-1}{i(\omega_1 + \omega_2) + E_{ik}} + \beta \delta_{\omega_1, -
\omega_2} \delta_{E_{ik}, 0} \right]\right\}
 \f$
 */

    inline static complex<double>
        PhiM_ii(double beta, double Ei, double Ej, double Ek, double El,
            complex<double>  w1,complex<double>  w2,complex<double> w3, double PHI_EPS)
    {
        complex<double> I_w3_Ekl    = 1./(w3+(Ek-El));

        // ii  Term 2,4 or 6
        complex<double> w23_Ejl     = w2+w3+(Ej-El);
        complex<double> I_w1_Eij    = 1./(w1+(Ei-Ej));
        complex<double> T24o6;

        if(abs(w23_Ejl) > PHI_EPS)
            T24o6 = 1./w23_Ejl*(I_w1_Eij-1./(w1+w2+w3+Ei-El));
        else
            T24o6 = I_w1_Eij*I_w1_Eij;


        // ii  Term 11 or 12
        complex<double> w12_Eik   = (w1+w2+(Ei-Ek));
        complex<double> T11o12;
        if(abs(w12_Eik) > PHI_EPS)
            T11o12 = -1./w12_Eik;
        else
            T11o12 = beta;

        return I_w3_Ekl*(T24o6-1./(w2+Ej-Ek)*(I_w1_Eij+T11o12));
    };

    inline static complex<double>
        PhiM_ji(double beta, double Ei, double Ej, double Ek, double El,
            complex<double>  w1,complex<double>  w2,complex<double> w3, double PHI_EPS)
    {
        complex<double> I_w3_Ekl    = 1./(w3+(Ek-El));
        complex<double> w23_Eil     = (w2+w3+Ei-El);
        complex<double> I_w1_Eji    = 1./(w1+Ej-Ei);

        complex<double> T1o57;
        if(abs(w23_Eil)>PHI_EPS)
            T1o57 = 1./w23_Eil* I_w1_Eji;
        else
            T1o57 = I_w1_Eji*I_w1_Eji - beta*I_w1_Eji;

        //                     T8
        return I_w3_Ekl*(T1o57 - 1./(w2+Ei-Ek)*I_w1_Eji);
    };

    inline static complex<double>
        PhiM_ki(double beta, double Ei, double Ej, double Ek, double El,
            complex<double>  w1,complex<double>  w2,complex<double> w3, double PHI_EPS)
    {
        complex<double> w12_Eki   = (w1+w2+Ek-Ei);
        if(abs(w12_Eki)>PHI_EPS)
            return -1./(w3+Ei-El)*1./w12_Eki*1./(w2+Ej-Ei);
        return 0.0;
    };

    inline static complex<double>
        PhiM_li(double beta, double Ei, double Ej, double Ek, double El,
            complex<double>  w1,complex<double>  w2,complex<double> w3, double PHI_EPS)
    {
        complex<double> w23_Eji   = (w2+w3+Ej-Ei);
        if(abs(w23_Eji)>PHI_EPS)
            return -1./( (w3+Ek-Ei) * w23_Eji * (w1+w2+w3+El-Ei) );
        return 0.0;
    };
// }}}


// **** Chi_loc ***************************************************************************
#ifdef __cplusplus
extern "C"{
            void   chi_tilde_loc_(int &op, // which one is calculated. If(op == 1 ): DNDN; If(op == 2):UPDN
                   complex<double> &w1, complex<double> &w2, complex<double> &w3, // external frequencies.
                   double &PHI_EPS, double &beta,
                   double &Z, // partition function.
                   double &gsE,

                   int &sites, // number of sites.
                   int &nup, int &ndn, // nup and ndn specifying the sector for state i.

                   double *cp_i_E,      int &dim_E_i,      // energy list of sector (nup,ndn) and length of this list.
                   double *cp_pup_E,    int &dim_E_pup,    // energy list of sector (nup+1, ndn) and its length.
                   double *cp_pdn_E,    int &dim_E_pdn,    // energy list of sector (nup, ndn+1) and its length.
                   double *cp_mup_E,    int &dim_E_mup,    // energy list of sector (nup-1, ndn) and its length.
                   double *cp_mdn_E,    int &dim_E_mdn,    // energy list of sec        tor (nup, ndn-1) and its length.
                   double *cp_p2dn_E,   int &dim_E_p2dn,   // energy list of sector (nup, ndn+2) and its length.
                   double *cp_m2dn_E,   int &dim_E_m2dn,   // energy list of sector (nup, ndn-2) and its length.
                   double *cp_puppdn_E, int &dim_E_puppdn, // energy list of sector (nup+1, ndn+1) and its length.
                   double *cp_muppdn_E, int &dim_E_muppdn, // energy list of sector (nup-1, ndn+1) and its length.
                   double *cp_pupmdn_E, int &dim_E_pupmdn, // energy list of sector (nup+1, ndn-1) and its length.
                   double *cp_mupmdn_E, int &dim_E_mupmdn, // energy list of sector (nup-1, ndn-1) and its length.

                   int &n1, int &n2, int &n3, int &n4, int &n5, int &n6, int &n7, int &n8, int &n9, int &n10, int &n11,
                   int &n12, int &n13, int &n14, int &n15, int &n16, int &n17, int &n18, int &n19, int &n20, int &n21, int &n22,
                   int &n23, int &n24, int &n25, int &n26, int &n27, int &n28,

                   double *cp_i_cdup,          // < nup+1,ndn   | c^dag_up | nup, ndn >
                   double *cp_i_cddn,          // < nup,ndn+1   | c^dag_dn | nup, ndn >
                   double *cp_pup_cddn,        // < nup+1,ndn+1 | c^dag_dn | nup+1, ndn >
                   double *cp_pdn_cdup,        // < nup+1,ndn+1 | c^dag_up | nup, ndn+1 >
                   double *cp_pdn_cddn,        // < nup,ndn+2   | c^dag_dn | nup, ndn+1 >
                   double *cp_mup_cdup,        // < nup,ndn     | c^dag_up | nup-1, ndn >
                   double *cp_mup_cddn,        // < nup-1,ndn+1 | c^dag_dn | nup-1, ndn >
                   double *cp_mdn_cddn,        // < nup,ndn     | c^dag_dn | nup, ndn-1 >
                   double *cp_mdn_cdup,        // < nup+1,ndn-1 | c^dag_up | nup, ndn-1 >
                   double *cp_muppdn_cdup,     // < nup,ndn+1   | c^dag_up | nup-1, ndn+1 >
                   double *cp_pupmdn_cddn,     // < nup+1,ndn   | c^dag_dn | nup+1, ndn-1 >
                   double *cp_mupmdn_cdup,     // < nup,ndn-1   | c^dag_up | nup-1, ndn-1 >
                   double *cp_mupmdn_cddn,     // < nup-1,ndn   | c^dag_dn | nup-1, ndn-1 >
                   double *cp_m2dn_cddn,       // < nup,ndn-1   | c^dag_dn | nup, ndn-2 >

                   complex<double> &pDNDN,complex<double> &pUPDN // to output result. 
                   );

}
#endif

void chi_tilde_loc_(int &op, // which one is calculated. If(op &1): DNDN; If(op &2):UPDN
                   complex<double> &w1, complex<double> &w2, complex<double> &w3, // external frequencies.
                   double &PHI_EPS, double &beta,
		   double &Z, // partition function.
                   double &gsE,

                   int &sites, // number of sites.
                   int &nup, int &ndn, // nup and ndn specifying the sector for state i.

                   double *cp_i_E,      int &dim_E_i,      // energy list of sector (nup,ndn) and length of this list.
                   double *cp_pup_E,    int &dim_E_pup,    // energy list of sector (nup+1, ndn) and its length.
                   double *cp_pdn_E,    int &dim_E_pdn,    // energy list of sector (nup, ndn+1) and its length.
                   double *cp_mup_E,    int &dim_E_mup,    // energy list of sector (nup-1, ndn) and its length.
                   double *cp_mdn_E,    int &dim_E_mdn,    // energy list of sec	tor (nup, ndn-1) and its length.
                   double *cp_p2dn_E,   int &dim_E_p2dn,   // energy list of sector (nup, ndn+2) and its length.
                   double *cp_m2dn_E,   int &dim_E_m2dn,   // energy list of sector (nup, ndn-2) and its length.
		   double *cp_puppdn_E, int &dim_E_puppdn, // energy list of sector (nup+1, ndn+1) and its length.
                   double *cp_muppdn_E, int &dim_E_muppdn, // energy list of sector (nup-1, ndn+1) and its length.
                   double *cp_pupmdn_E, int &dim_E_pupmdn, // energy list of sector (nup+1, ndn-1) and its length.
                   double *cp_mupmdn_E, int &dim_E_mupmdn, // energy list of sector (nup-1, ndn-1) and its length.

                   int &n1, int &n2, int &n3, int &n4, int &n5, int &n6, int &n7, int &n8, int &n9, int &n10, int &n11, 
                   int &n12, int &n13, int &n14, int &n15, int &n16, int &n17, int &n18, int &n19, int &n20, int &n21, int &n22, 
                   int &n23, int &n24, int &n25, int &n26, int &n27, int &n28,

                   double *cp_i_cdup,          // < nup+1,ndn   | c^dag_up | nup, ndn >
                   double *cp_i_cddn,          // < nup,ndn+1   | c^dag_dn | nup, ndn >
                   double *cp_pup_cddn,        // < nup+1,ndn+1 | c^dag_dn | nup+1, ndn >
                   double *cp_pdn_cdup,        // < nup+1,ndn+1 | c^dag_up | nup, ndn+1 >
                   double *cp_pdn_cddn,        // < nup,ndn+2   | c^dag_dn | nup, ndn+1 >
                   double *cp_mup_cdup,        // < nup,ndn     | c^dag_up | nup-1, ndn >
                   double *cp_mup_cddn,        // < nup-1,ndn+1 | c^dag_dn | nup-1, ndn >
                   double *cp_mdn_cddn,        // < nup,ndn     | c^dag_dn | nup, ndn-1 >
                   double *cp_mdn_cdup,        // < nup+1,ndn-1 | c^dag_up | nup, ndn-1 >
                   double *cp_muppdn_cdup,     // < nup,ndn+1   | c^dag_up | nup-1, ndn+1 >
                   double *cp_pupmdn_cddn,     // < nup+1,ndn   | c^dag_dn | nup+1, ndn-1 >
                   double *cp_mupmdn_cdup,     // < nup,ndn-1   | c^dag_up | nup-1, ndn-1 >
                   double *cp_mupmdn_cddn,     // < nup-1,ndn   | c^dag_dn | nup-1, ndn-1 >
                   double *cp_m2dn_cddn,       // < nup,ndn-1   | c^dag_dn | nup, ndn-2 >

                   /*   The indices of the matrix should be consistent with those of lists of energy eigenvalues. 
                        
                        Thus only the matrix of creation operators, c^dag_up and c^dag_dn are used. 
			The matrix elements are stored in a 1D array, the convention is as the following examples:
                        Eg. 1: if | j > is a state in sector (nup, ndn), 
                               then < i | c^dag_up | j > =  *cp_i_cdup[dimj*i+j], where dimj=dim_E_i is the number of states in sector (nup, ndn). 
                        Eg. 2: if | j > is a state in sector (nup-1, ndn), 
                               then < i | c^dag_dn | j > =  cp_mup_cddn[dimj*i+j], where dimj=dim_E_mup is the number of states in sector (nup-1, ndn).
                        Eg. 3: if | j > is a state in sector (nup+1, ndn-1), 
			       then < i | c^dag_dn | j > =  cp_pupmdn_cddn[dimj*i+j], where dimj=dim_E_pupmdn is the number of states in sector (nup+1, ndn-1).

                   */

                   complex<double> &pDNDN,complex<double> &pUPDN // to output result. 
                   )
{
//******* check scalar passing. **************
//    cout<<"op in C++ = "<<op<<endl;
//    cout<<"w1 w2 w3 in C++ = "<<w1<<w2<<w3<<endl;
//    cout<<"PHI_EPS, beta, Z in C++ = "<<PHI_EPS<<" "<<beta<<" "<<" "<<Z<<endl;
//    cout<<"sites in C++ = "<<sites<<endl;
//    cout<<"w1, w2, w3 in C++ = "<<w1<<" "<<w2<<" "<<w3<<endl;
//    cout<<"dim_E's in C++ = "<<dim_E_i<<" "<<dim_E_pup<<" "<<dim_E_pdn<<" "<<dim_E_mup<<" "<<dim_E_mdn<<" "<<dim_E_p2dn<<" "<<dim_E_m2dn<<" "
//                             <<dim_E_puppdn<<" "<<dim_E_muppdn<<" "<<dim_E_pupmdn<<" "<<dim_E_mupmdn<<" "<<endl;
//    cout<<"cp_pup_cddn in C++ = "<<cp_pup_cddn[1]<<endl;

//******* print out eigenvalues and matrix to check if they are correctly passed. *****
//    ofstream file("cpp_eigens.dat", ios_base::app);
//    file<<"partition function: "<<Z<<endl;
//    file<<"beta = "<<beta<<endl;
//    file<<"nsites:"<<sites<<endl;
//    file<<"##########################################"<<endl;
//    file<<"nsites="<<sites<<", nup="<<nup<<", ndn="<<ndn<<endl;

//    file<<showpoint<<fixed<<right<<setprecision(4);
//    file<<"cp_i_E="<<endl;    for (int i=0; i<n2; i++){       file<<setw(8)<<cp_i_E[i]<<" "; if ((i+1)%6==0) file<<endl;}    file<<endl;
//    file<<"cp_pup_E="<<endl;    for (int i=0; i<n1; i++){     file<<setw(8)<<cp_pup_E[i]<<" "; if ((i+1)%6==0) file<<endl;}    file<<endl;
//    file<<"cp_pdn_E="<<endl;    for (int i=0; i<n3; i++){     file<<setw(8)<<cp_pdn_E[i]<<" "; if ((i+1)%6==0) file<<endl;}   file<<endl;
//    file<<"cp_mup_E="<<endl;    for (int i=0; i<n12; i++){    file<<setw(8)<<cp_mup_E[i]<<" "; if ((i+1)%6==0) file<<endl;}    file<<endl;
//    file<<"cp_mdn_E="<<endl;    for (int i=0; i<n16; i++){    file<<setw(8)<<cp_mdn_E[i]<<" "; if ((i+1)%6==0) file<<endl;}    file<<endl;
//    file<<"cp_p2dn_E="<<endl;    for (int i=0; i<n9; i++){    file<<setw(8)<<cp_p2dn_E[i]<<" "; if ((i+1)%6==0) file<<endl;}    file<<endl;
//    file<<"cp_m2dn_E="<<endl;    for (int i=0; i<n28; i++){   file<<setw(8)<<cp_m2dn_E[i]<<" "; if ((i+1)%6==0) file<<endl;}    file<<endl;
//    file<<"cp_puppdn_E="<<endl;    for (int i=0; i<n5; i++){  file<<setw(8)<<cp_puppdn_E[i]<<" "; if ((i+1)%6==0) file<<endl;}    file<<endl;
//    file<<"cp_muppdn_E="<<endl;    for (int i=0; i<n13; i++){ file<<setw(8)<<cp_muppdn_E[i]<<" "; if ((i+1)%6==0) file<<endl;}    file<<endl;
//    file<<"cp_pupmdn_E="<<endl;    for (int i=0; i<n17; i++){ file<<setw(8)<<cp_pupmdn_E[i]<<" "; if ((i+1)%6==0) file<<endl;}   file<<endl;
//    file<<"cp_mupmdn_E="<<endl;  for (int i=0; i<n24; i++){   file<<setw(8)<<cp_mupmdn_E[i]<<" "; if ((i+1)%6==0) file<<endl;}   file<<endl;

//    file<<"cp_i_cdup="<<endl; for (int i=0; i<n1*n2; i++ ){        file<<setw(8)<<cp_i_cdup[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_i_cddn="<<endl; for (int i=0; i<n3*n4; i++ ){        file<<setw(8)<<cp_i_cddn[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_pup_cddn="<<endl; for (int i=0; i<n5*n6; i++ ){      file<<setw(8)<<cp_pup_cddn[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_pdn_cdup="<<endl; for (int i=0; i<n7*n8; i++ ){      file<<setw(8)<<cp_pdn_cdup[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_pdn_cddn="<<endl; for (int i=0; i<n9*n10; i++ ){     file<<setw(8)<<cp_pdn_cddn[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_mup_cdup="<<endl; for (int i=0; i<n11*n12; i++ ){    file<<setw(8)<<cp_mup_cdup[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_mup_cddn="<<endl; for (int i=0; i<n13*n14; i++ ){    file<<setw(8)<<cp_mup_cddn[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_mdn_cddn="<<endl; for (int i=0; i<n15*n16; i++ ){    file<<setw(8)<<cp_mdn_cddn[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_mdn_cdup="<<endl; for (int i=0; i<n17*n18; i++ ){    file<<setw(8)<<cp_mdn_cdup[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_muppdn_cdup="<<endl; for (int i=0; i<n19*n20; i++ ){ file<<setw(8)<<cp_muppdn_cdup[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_pupmdn_cddn="<<endl; for (int i=0; i<n21*n22; i++ ){ file<<setw(8)<<cp_pupmdn_cddn[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_mupmdn_cdup="<<endl; for (int i=0; i<n23*n24; i++ ){ file<<setw(8)<<cp_mupmdn_cdup[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_mupmdn_cddn="<<endl; for (int i=0; i<n25*n26; i++ ){ file<<setw(8)<<cp_mupmdn_cddn[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;
//    file<<"cp_m2dn_cddn="<<endl; for (int i=0; i<n27*n28; i++ ){   file<<setw(8)<<cp_m2dn_cddn[i]<<" "; if ((i+1)%7==0) file<<endl;}  file<<endl;

//    file.close();

//*************************************************************************************
    pDNDN = complex<double>(0.,0.);
    pUPDN = complex<double>(0.,0.);

    // temporary variables
    complex<double> xi1 = complex<double>(0.,0.); complex<double> yi1 = complex<double>(0.,0.);
    complex<double> xi2 = complex<double>(0.,0.); complex<double> yi2 = complex<double>(0.,0.);

    for(int stati=0; stati < dim_E_i ; stati++) // summation over i.
    {
        double boltzZ = exp(-(cp_i_E[stati]-gsE)*beta)/Z;
//        cout<<"boltzZ="<<boltzZ<<endl;

        int dimi = dim_E_i;
        int diml, dimk, dimj;

// **** DNDN ******************************************************************************
// ****************************************************************************************

// For DNDN the following terms are important:
// tot  perm    w   w   w         matrix elements      fermi       //  Indices:
//   0   0      1   2   3         -d  +d  -d  +d   ii    1         //    ii :   i j k l i
//   1   0      1   2   3         +d  -d  +d  -d   ji    1         //    ji :   i k l j i
//   2   0      1   2   3         -d  +d  -d  +d   ki    1         //    ki :   i l k j i
//   3   0      1   2   3         +d  -d  +d  -d   li    1         //    li :   i l j k i

//   4   1      2   1   3         +d  -d  -d  +d   ii   -1
//   5   1      2   1   3         -d  -d  +d  +d   ji   -1
//   6   1      2   1   3         -d  +d  +d  -d   ki   -1
//   7   1      2   1   3         +d  +d  -d  -d   li   -1

//   8   2      3   1   2         -d  -d  +d  +d   ii    1
//   9   2      3   1   2         -d  +d  +d  -d   ji    1
//  10   2      3   1   2         +d  +d  -d  -d   ki    1
//  11   2      3   1   2         +d  -d  -d  +d   li    1

//  12   3      1   3   2         -d  -d  +d  +d   ii   -1
//  13   3      1   3   2         -d  +d  +d  -d   ji   -1
//  14   3      1   3   2         +d  +d  -d  -d   ki   -1
//  15   3      1   3   2         +d  -d  -d  +d   li   -1

//  16   4      2   3   1         +d  -d  -d  +d   ii    1
//  17   4      2   3   1         -d  -d  +d  +d   ji    1
//  18   4      2   3   1         -d  +d  +d  -d   ki    1
//  19   4      2   3   1         +d  +d  -d  -d   li    1

//  20   5      3   2   1         -d  +d  -d  +d   ii   -1
//  21   5      3   2   1         +d  -d  +d  -d   ji   -1
//  22   5      3   2   1         -d  +d  -d  +d   ki   -1
//  23   5      3   2   1         +d  -d  +d  -d   li   -1

    if(op == 1)
    {
        // initialize
        complex<double> cDNDN = complex<double>(0.,0.);

        if(ndn!=sites)
        {   
            diml = dim_E_pdn ;
            dimk = dim_E_i;
            dimj = dim_E_pdn;

//            cout<<"ndn!=nsites\n";
//            cout<<"case 0/20 and 2/22 invoked: \n";
//            cout<<"diml=dim_E_pdn, dimk=dim_E_i, dimj=dim_E_pdn\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ii(beta,cp_i_E[stati], cp_pdn_E[j],cp_i_E[k],cp_pdn_E[l],w1,w2,w3,PHI_EPS);
                xi2 =-PhiM_ii(beta,cp_i_E[stati], cp_pdn_E[j],cp_i_E[k],cp_pdn_E[l],w3,w2,w1,PHI_EPS);
                // tot =0/20          <stati|cdn|j>          *   <j|c+dn|k>        * <k|cdn|l>            * <l|c+dn|stati>
                cDNDN +=  (xi1+xi2)* cp_i_cddn[dimi*j+stati]* cp_i_cddn[dimk*j+k] * cp_i_cddn[dimk*l+k] * cp_i_cddn[dimi*l+stati];
            }


            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                yi1 = PhiM_ki(beta,cp_i_E[stati], cp_pdn_E[j],cp_i_E[k],cp_pdn_E[l],w1,w2,w3,PHI_EPS);
                yi2 =-PhiM_ki(beta,cp_i_E[stati], cp_pdn_E[j],cp_i_E[k],cp_pdn_E[l],w3,w2,w1,PHI_EPS);
                // tot = 2/22       <stati|cdn|l>          *   <l|c+dn|k>        * <k|cdn|j>            * <j|c+dn|stati>
                cDNDN +=  (yi1+yi2)* cp_i_cddn[dimi*l+stati]* cp_i_cddn[dimk*l+k] * cp_i_cddn[dimk*j+k] * cp_i_cddn[dimi*j+stati];
            }
//           cout<<"cDNDN="<<cDNDN<<endl;
        }


        if(ndn!=sites && ndn!=0)
        {
            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;
 
//            cout<<"ndn!=sites && ndn!=0\n";
//            cout<<"case 4/16 invoked:\n";
//            cout<<"diml=dim_E_pdn, dimk=dim_E_i, dimj=dim_E_mdn\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 =-PhiM_ii(beta,cp_i_E[stati], cp_mdn_E[j],cp_i_E[k],cp_pdn_E[l],w2,w1,w3,PHI_EPS);
                xi2 = PhiM_ii(beta,cp_i_E[stati], cp_mdn_E[j],cp_i_E[k],cp_pdn_E[l],w2,w3,w1,PHI_EPS);
                // tot =4/16          <stati|c+dn|j>          *   <j|cdn|k>        * <k|cdn|l>            * <l|c+dn|stati>
                cDNDN +=  (xi1+xi2)* cp_mdn_cddn[dimj*stati+j]* cp_mdn_cddn[dimj*k+j] * cp_i_cddn[dimk*l+k] * cp_i_cddn[dimi*l+stati];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;
             
            diml = dim_E_mdn;
            dimk = dim_E_pdn;
            dimj = dim_E_i;

//            cout<<"case 11/15 invoked:\n";
//            cout<<"diml=dim_E_mdn, dimk=dim_pdn, dimj=dim_E_i\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;
//            cout<<"\n";
 
            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                yi1 = PhiM_li(beta,cp_i_E[stati], cp_i_E[j],cp_pdn_E[k],cp_mdn_E[l],w3,w1,w2,PHI_EPS);
                yi2 =-PhiM_li(beta,cp_i_E[stati], cp_i_E[j],cp_pdn_E[k],cp_mdn_E[l],w1,w3,w2,PHI_EPS);
                // tot = 11/15       <stati|c+dn|l>          *   <l|cdn|j>        * <j|cdn|k>            * <k|c+dn|stati>
                cDNDN +=  (yi1+yi2)* cp_mdn_cddn[diml*stati+l]* cp_mdn_cddn[diml*j+l] * cp_i_cddn[dimj*k+j] * cp_i_cddn[dimi*k+stati];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;

            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;

//            cout<<"case 6/18 invoked:\n";
//            cout<<"diml=dim_E_pdn, dimk=dim_E_i, dimj=dim_E_mdn\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;
//            cout<<"\n";

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 =-PhiM_ki(beta,cp_i_E[stati], cp_mdn_E[j],cp_i_E[k],cp_pdn_E[l],w2,w1,w3,PHI_EPS);
                xi2 = PhiM_ki(beta,cp_i_E[stati], cp_mdn_E[j],cp_i_E[k],cp_pdn_E[l],w2,w3,w1,PHI_EPS);
                // tot =6/18          <stati|cdn|l>          *   <l|c+dn|k>        * <k|c+dn|j>            * <j|cdn|stati>
                cDNDN +=  (xi1+xi2)* cp_i_cddn[dimi*l+stati]* cp_i_cddn[dimk*l+k] * cp_mdn_cddn[dimj*k+j] * cp_mdn_cddn[dimj*stati+j];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;

            diml = dim_E_i;
            dimk = dim_E_pdn;
            dimj = dim_E_mdn;

//            cout<<"case 9/13 invoked:\n";
//            cout<<"diml=dim_E_i, dimk=dim_E_pdn, dimj=dim_E_mdn\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;
//            cout<<"\n";

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                yi1 = PhiM_ji(beta,cp_i_E[stati], cp_mdn_E[j],cp_pdn_E[k],cp_i_E[l],w3,w1,w2,PHI_EPS);
                yi2 =-PhiM_ji(beta,cp_i_E[stati], cp_mdn_E[j],cp_pdn_E[k],cp_i_E[l],w1,w3,w2,PHI_EPS);
                // tot = 9/13       <stati|cdn|k>          *   <k|c+dn|l>        * <l|c+dn|j>            * <j|cdn|stati>
                cDNDN +=  (yi1+yi2)* cp_i_cddn[dimi*k+stati]* cp_i_cddn[diml*k+l] * cp_mdn_cddn[dimj*l+j] * cp_mdn_cddn[dimj*stati+j];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;
        }

        if(ndn<sites-1) // 5,8,12,17
        {
            diml = dim_E_p2dn;
            dimk = dim_E_pdn;
            dimj = dim_E_pdn;

//            cout<<"ndn<nsites-1\n";
//            cout<<"case 5/17 invoked: \n";
//            cout<<"diml=dim_E_p2dn, dimk=dim_E_pdn, dimj=dim_E_pdn\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                yi1 =-PhiM_ji(beta,cp_i_E[stati], cp_pdn_E[j],cp_pdn_E[k],cp_p2dn_E[l],w2,w1,w3,PHI_EPS);
                yi2 = PhiM_ji(beta,cp_i_E[stati], cp_pdn_E[j],cp_pdn_E[k],cp_p2dn_E[l],w2,w3,w1,PHI_EPS);
                // tot = 5/17       <stati|cdn|k>          *   <k|cdn|l>        * <l|c+dn|j>            * <j|c+dn|stati>
                cDNDN +=  (yi1+yi2)* cp_i_cddn[dimi*k+stati]* cp_pdn_cddn[dimk*l+k] * cp_pdn_cddn[dimj*l+j] * cp_i_cddn[dimi*j+stati];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;

            diml = dim_E_pdn;
            dimk = dim_E_p2dn;
            dimj = dim_E_pdn;

//            cout<<"case 8/12 invoked: \n";
//            cout<<"diml=dim_E_pdn, dimk=dim_E_p2dn, dimj=dim_E_pdn\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;
 
            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ii(beta,cp_i_E[stati], cp_pdn_E[j],cp_p2dn_E[k],cp_pdn_E[l],w3,w1,w2,PHI_EPS);
                xi2 =-PhiM_ii(beta,cp_i_E[stati], cp_pdn_E[j],cp_p2dn_E[k],cp_pdn_E[l],w1,w3,w2,PHI_EPS);
                // tot =8/12          <stati|cdn|j>          *   <j|cdn|k>        * <k|c+dn|l>            * <l|c+dn|stati>
                cDNDN +=  (xi1+xi2)* cp_i_cddn[dimi*j+stati]* cp_pdn_cddn[dimj*k+j] * cp_pdn_cddn[diml*k+l] * cp_i_cddn[dimi*l+stati];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;
        }

        if(ndn>1) // 7,10,14,19
        {
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_m2dn;

//            cout<<"ndn>1\n";
//            cout<<"case 7/19 invoked: \n";
//            cout<<"diml=dim_E_mdn, dimk=dim_E_mdn, dimj=dim_E_m2dn\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                yi1 =-PhiM_li(beta,cp_i_E[stati], cp_m2dn_E[j],cp_mdn_E[k],cp_mdn_E[l],w2,w1,w3,PHI_EPS);
                yi2 = PhiM_li(beta,cp_i_E[stati], cp_m2dn_E[j],cp_mdn_E[k],cp_mdn_E[l],w2,w3,w1,PHI_EPS);
                // tot = 7/19       <stati|c+dn|l>          *   <l|c+dn|j>        * <j|cdn|k>            * <k|cdn|stati>
                cDNDN +=  (yi1+yi2)* cp_mdn_cddn[diml*stati+l]* cp_m2dn_cddn[dimj*l+j] * cp_m2dn_cddn[dimj*k+j] * cp_mdn_cddn[dimk*stati+k];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;

            diml = dim_E_mdn;
            dimk = dim_E_m2dn;
            dimj = dim_E_mdn;

//            cout<<"case 10/14 invoked: \n";
//            cout<<"diml=dim_E_mdn, dimk=dim_E_m2dn, dimj=dim_E_mdn\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                yi1 = PhiM_ki(beta,cp_i_E[stati], cp_mdn_E[j],cp_m2dn_E[k],cp_mdn_E[l],w3,w1,w2,PHI_EPS);
                yi2 =-PhiM_ki(beta,cp_i_E[stati], cp_mdn_E[j],cp_m2dn_E[k],cp_mdn_E[l],w1,w3,w2,PHI_EPS);
                // tot = 10/14       <stati|c+dn|l>          *   <l|c+dn|k>        * <k|cdn|j>            * <j|cdn|stati>
                cDNDN +=  (yi1+yi2)* cp_mdn_cddn[diml*stati+l]* cp_m2dn_cddn[dimk*l+k] * cp_m2dn_cddn[dimk*j+k] * cp_mdn_cddn[dimj*stati+j];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;
        }

        if(ndn!=0)
        {
            diml = dim_E_i;
            dimk = dim_E_mdn;
            dimj = dim_E_mdn;

//            cout<<"ndn!=0\n";
//            cout<<"case 1/21 invoked: \n";
//            cout<<"diml=dim_E_i, dimk=dim_E_mdn, dimj=dim_E_mdn\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ji(beta,cp_i_E[stati], cp_mdn_E[j],cp_mdn_E[k],cp_i_E[l],w1,w2,w3,PHI_EPS);
                xi2 =-PhiM_ji(beta,cp_i_E[stati], cp_mdn_E[j],cp_mdn_E[k],cp_i_E[l],w3,w2,w1,PHI_EPS);
                // tot=1/21             <stati|c+dn|k>         *  <k|cdn|l>             * <l|c+dn|j>             * <j|cdn|stati>
                cDNDN +=  (xi1+xi2)* cp_mdn_cddn[dimk*stati+k] * cp_mdn_cddn[dimk*l+k] * cp_mdn_cddn[dimj*l+j] * cp_mdn_cddn[dimj*stati+j];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;

            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_i;

//            cout<<"case 3/12 invoked: \n";
//            cout<<"diml=dim_E_mdn, dimk=dim_E_mdn, dimj=dim_E_i\n";
//            cout<<"diml="<<diml<<", dimk="<<dimk<<", dimj="<<dimj<<endl;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                yi1 = PhiM_li(beta,cp_i_E[stati], cp_i_E[j],cp_mdn_E[k],cp_mdn_E[l],w1,w2,w3,PHI_EPS);
                yi2 =-PhiM_li(beta,cp_i_E[stati], cp_i_E[j],cp_mdn_E[k],cp_mdn_E[l],w3,w2,w1,PHI_EPS);
                // tot=3/23           <stati|c+dn|l>           *  <l|cdn|j>             * <j|c+dn|k>             * <k|cdn|stati>
                cDNDN +=  (yi1+yi2)* cp_mdn_cddn[diml*stati+l] * cp_mdn_cddn[diml*j+l] * cp_mdn_cddn[dimk*j+k] * cp_mdn_cddn[dimk*stati+k];
            }
//            cout<<"cDNDN="<<cDNDN<<endl;
        }
        pDNDN += cDNDN*boltzZ;
//        cout<<"pDNDN = "<<pDNDN<<endl;
    }

// **** UPDN *****************************************************************************
// ***************************************************************************************

// For UPDN *all* terms are important:
// tot  perm    w   w   w         matrix elements      fermi       //  Indices:
//   0   0      1   2   3         -u  +u  -d  +d   ii    1         //    ii :   i j k l i
//   1   0      1   2   3         +u  -d  +d  -u   ji    1         //    ji :   i k l j i
//   2   0      1   2   3         -d  +d  -u  +u   ki    1         //    ki :   i l k j i
//   3   0      1   2   3         +d  -u  +u  -d   li    1         //    li :   i l j k i

//   4   1      2   1   3         +u  -u  -d  +d   ii   -1
//   5   1      2   1   3         -u  -d  +d  +u   ji   -1
//   6   1      2   1   3         -d  +d  +u  -u   ki   -1
//   7   1      2   1   3         +d  +u  -u  -d   li   -1

//   8   2      3   1   2         -d  -u  +u  +d   ii    1
//   9   2      3   1   2         -u  +u  +d  -d   ji    1
//  10   2      3   1   2         +u  +d  -d  -u   ki    1
//  11   2      3   1   2         +d  -d  -u  +u   li    1

//  12   3      1   3   2         -u  -d  +u  +d   ii   -1
//  13   3      1   3   2         -d  +u  +d  -u   ji   -1
//  14   3      1   3   2         +u  +d  -u  -d   ki   -1
//  15   3      1   3   2         +d  -u  -d  +u   li   -1

//  16   4      2   3   1         +u  -d  -u  +d   ii    1
//  17   4      2   3   1         -d  -u  +d  +u   ji    1
//  18   4      2   3   1         -u  +d  +u  -d   ki    1
//  19   4      2   3   1         +d  +u  -d  -u   li    1

//  20   5      3   2   1         -d  +u  -u  +d   ii   -1
//  21   5      3   2   1         +u  -u  +d  -d   ji   -1
//  22   5      3   2   1         -u  +d  -d  +u   ki   -1
//  23   5      3   2   1         +d  -d  +u  -u   li   -1

    if(op == 2)
    {
//        cout<<"compute UPDN case"<<endl;
        // initialize
        complex<double> cUPDN = complex<double>(0.,0.);

        if(ndn!=sites && nup != sites)  // tot: 0,2,5,8,12,17
        {
            // tot 0,2
            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_pup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ii(beta,cp_i_E[stati], cp_pup_E[j],cp_i_E[k],cp_pdn_E[l],w1,w2,w3,PHI_EPS);
                xi2 = PhiM_ki(beta,cp_i_E[stati], cp_pup_E[j],cp_i_E[k],cp_pdn_E[l],w1,w2,w3,PHI_EPS);

// tot=0              +            <stati|cup|j>          *  <j|c+up|k>         * <k|cdn|l>           * <l|c+dn|stati>
// tot=2              +            <stati|cdn|l>          *  <l|c+dn|k>         * <k|cup|j>           * <j|c+up|stati>
                cUPDN +=(xi1+xi2) * cp_i_cdup[dimi*j+stati]* cp_i_cdup[dimk*j+k] * cp_i_cddn[dimk*l+k] * cp_i_cddn[dimi*l+stati];
            }

            // tot 5
            diml = dim_E_puppdn;
            dimk = dim_E_pup;
            dimj = dim_E_pup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ji(beta,cp_i_E[stati], cp_pup_E[j],cp_pup_E[k],cp_puppdn_E[l],w2,w1,w3,PHI_EPS);

// tot=5              -            <stati|cup|k>          *  <k|cdn|l>         * <l|c+dn|j>           * <j|c+up|stati>
                cUPDN -=xi1 * cp_i_cdup[dimi*k+stati]* cp_pup_cddn[dimk*l+k] * cp_pup_cddn[dimj*l+j] * cp_i_cdup[dimi*j+stati];
            }

            // tot 8
            diml = dim_E_pdn;
            dimk = dim_E_puppdn;
            dimj = dim_E_pdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ii(beta,cp_i_E[stati], cp_pdn_E[j],cp_puppdn_E[k],cp_pdn_E[l],w3,w1,w2,PHI_EPS);

// tot=8              +      <stati|cdn|j>          *  <j|cup|k>            * <k|c+up|l>            * <l|c+dn|stati>
                cUPDN +=xi1 * cp_i_cddn[dimi*j+stati]* cp_pdn_cdup[dimj*k+j] * cp_pdn_cdup[diml*k+l] * cp_i_cddn[dimi*l+stati];
            }


            // tot 12
            diml = dim_E_pdn;
            dimk = dim_E_puppdn;
            dimj = dim_E_pup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ii(beta,cp_i_E[stati], cp_pup_E[j],cp_puppdn_E[k],cp_pdn_E[l],w1,w3,w2,PHI_EPS);

// tot=12             -            <stati|cup|j>          *  <j|cdn|k>         * <k|c+up|l>           * <l|c+dn|stati>
                cUPDN -=xi1 * cp_i_cdup[dimi*j+stati]* cp_pup_cddn[dimj*k+j] * cp_pdn_cdup[diml*k+l] * cp_i_cddn[dimi*l+stati];
            }


            // tot 17
            diml = dim_E_puppdn;
            dimk = dim_E_pdn;
            dimj = dim_E_pup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ji(beta,cp_i_E[stati], cp_pup_E[j],cp_pdn_E[k],cp_puppdn_E[l],w2,w3,w1,PHI_EPS);

// tot=17             +            <stati|cdn|k>          *  <k|cup|l>         * <l|c+dn|j>           * <j|c+up|stati>
                cUPDN +=xi1 * cp_i_cddn[dimi*k+stati]* cp_pdn_cdup[dimk*l+k] * cp_pup_cddn[dimj*l+j] * cp_i_cdup[dimi*j+stati];
            }
        }


        if(ndn!=sites && nup != 0)  // tot: 1,4,6,13,16,20
        {
            // tot 1
            diml = dim_E_muppdn;
            dimk = dim_E_mup;
            dimj = dim_E_mup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ji(beta,cp_i_E[stati], cp_mup_E[j],cp_mup_E[k],cp_muppdn_E[l],w1,w2,w3,PHI_EPS);

// tot=1              +         <stati|c+up|k>      *  <k|cdn|l>            * <l|c+dn|j>            * <j|cup|stati>
                cUPDN +=xi1 * cp_mup_cdup[dimk*stati+k]* cp_mup_cddn[dimk*l+k] * cp_mup_cddn[dimj*l+j] * cp_mup_cdup[dimj*stati+j];
            }


            // tot 4
            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ii(beta,cp_i_E[stati], cp_mup_E[j],cp_i_E[k],cp_pdn_E[l],w2,w1,w3,PHI_EPS);

// tot=4              -         <stati|c+up|j>      *  <j|cup|k>            * <k|cdn|l>           * <l|c+dn|stati>
                cUPDN -=xi1 * cp_mup_cdup[dimj*stati+j]* cp_mup_cdup[dimj*k+j] * cp_i_cddn[dimk*l+k] * cp_i_cddn[dimi*l+stati];
            }


            // tot 6
            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ki(beta,cp_i_E[stati], cp_mup_E[j],cp_i_E[k],cp_pdn_E[l],w2,w1,w3,PHI_EPS);

// tot=6              -         <stati|cdn|l>      *  <l|c+dn|k>            * <k|c+up|j>            * <j|cup|stati>
                cUPDN -=xi1 * cp_i_cddn[dimi*l+stati]* cp_i_cddn[dimk*l+k] * cp_mup_cdup[dimj*k+j] * cp_mup_cdup[dimj*stati+j];
            }

            // tot 13
            diml = dim_E_muppdn;
            dimk = dim_E_pdn;
            dimj = dim_E_mup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ji(beta,cp_i_E[stati], cp_mup_E[j],cp_pdn_E[k],cp_muppdn_E[l],w1,w3,w2,PHI_EPS);

// tot=13             -         <stati|cdn|k>      *  <k|c+up|l>            * <l|c+dn|j>            * <j|cup|stati>
                cUPDN -=xi1 * cp_i_cddn[dimi*k+stati]* cp_muppdn_cdup[diml*k+l] * cp_mup_cddn[dimj*l+j] * cp_mup_cdup[dimj*stati+j];
            }


            // tot 16
            diml = dim_E_pdn;
            dimk = dim_E_muppdn;
            dimj = dim_E_mup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ii(beta,cp_i_E[stati], cp_mup_E[j],cp_muppdn_E[k],cp_pdn_E[l],w2,w3,w1,PHI_EPS);

// tot=16             +         <stati|c+up|j>      *  <j|cdn|k>            * <k|cup|l>            * <l|c+dn|stati>
                cUPDN +=xi1 * cp_mup_cdup[dimj*stati+j]* cp_mup_cddn[dimj*k+j] * cp_muppdn_cdup[dimk*l+k] * cp_i_cddn[dimi*l+stati];
            }


            // tot 20
            diml = dim_E_pdn;
            dimk = dim_E_muppdn;
            dimj = dim_E_pdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ii(beta,cp_i_E[stati], cp_pdn_E[j],cp_muppdn_E[k],cp_pdn_E[l],w3,w2,w1,PHI_EPS);

// tot=20             -         <stati|cdn|j>      *  <j|c+up|k>            * <k|cup|l>          * <l|c+dn|stati>
                cUPDN -=xi1 * cp_i_cddn[dimi*j+stati]* cp_muppdn_cdup[dimk*j+k] * cp_muppdn_cdup[dimk*l+k] * cp_i_cddn[dimi*l+stati];
            }
        } // 1,4,6,13,16,20






        if(ndn!=0 && nup != sites)  // tot: 3,9,11,15,18,22
        {
            // tot 3
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_pupmdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_li(beta,cp_i_E[stati], cp_pupmdn_E[j],cp_mdn_E[k],cp_mdn_E[l],w1,w2,w3,PHI_EPS);

// tot=3              +         <stati|c+dn|l>      *  <l|cup|j>            * <j|c+up|k>            * <k|cdn|stati>
                cUPDN +=xi1 * cp_mdn_cddn[diml*stati+l]* cp_mdn_cdup[diml*j+l] * cp_mdn_cdup[dimk*j+k] * cp_mdn_cddn[dimk*stati+k];
            }


            // tot 9
            diml = dim_E_i;
            dimk = dim_E_pup;
            dimj = dim_E_mdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ji(beta,cp_i_E[stati], cp_mdn_E[j],cp_pup_E[k],cp_i_E[l],w3,w1,w2,PHI_EPS);

// tot=16             +         <stati|cup|k>      *  <k|c+up|l>            * <l|c+dn|j>            * <j|cdn|stati>
                cUPDN +=xi1 * cp_i_cdup[dimi*k+stati]* cp_i_cdup[diml*k+l] * cp_mdn_cddn[dimj*l+j] * cp_mdn_cddn[dimj*stati+j];
            }


            // tot 11
            diml = dim_E_mdn;
            dimk = dim_E_pup;
            dimj = dim_E_i;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_li(beta,cp_i_E[stati], cp_i_E[j],cp_pup_E[k],cp_mdn_E[l],w3,w1,w2,PHI_EPS);

// tot=11             +         <stati|c+dn|l>      *  <l|cdn|j>            * <j|cup|k>          * <k|c+up|stati>
                cUPDN +=xi1 * cp_mdn_cddn[diml*stati+l]* cp_mdn_cddn[diml*j+l] * cp_i_cdup[dimj*k+j] * cp_i_cdup[dimi*k+stati];
            }


            // tot 15
            diml = dim_E_mdn;
            dimk = dim_E_pup;
            dimj = dim_E_pupmdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_li(beta,cp_i_E[stati], cp_pupmdn_E[j],cp_pup_E[k],cp_mdn_E[l],w1,w3,w2,PHI_EPS);

// tot=15             -         <stati|c+dn|l>      *  <l|cup|j>            * <j|cdn|k>            * <k|c+up|stati>
                cUPDN -=xi1 * cp_mdn_cddn[diml*stati+l]* cp_mdn_cdup[diml*j+l] * cp_pupmdn_cddn[dimj*k+j] * cp_i_cdup[dimi*k+stati];
            }


            // tot 18
            diml = dim_E_pup;
            dimk = dim_E_pupmdn;
            dimj = dim_E_mdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ki(beta,cp_i_E[stati], cp_mdn_E[j],cp_pupmdn_E[k],cp_pup_E[l],w2,w3,w1,PHI_EPS);

// tot=18             +         <stati|cup|l>      *  <l|c+dn|k>            * <k|c+up|j>            * <j|cdn|stati>
                cUPDN +=xi1 * cp_i_cdup[dimi*l+stati]* cp_pupmdn_cddn[dimk*l+k] * cp_mdn_cdup[dimj*k+j] * cp_mdn_cddn[dimj*stati+j];
            }


            // tot 22
            diml = dim_E_pup;
            dimk = dim_E_pupmdn;
            dimj = dim_E_pup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ki(beta,cp_i_E[stati], cp_pup_E[j],cp_pupmdn_E[k],cp_pup_E[l],w3,w2,w1,PHI_EPS);

// tot=22             -         <stati|cup|l>         *  <l|c+dn|k>            * <k|cdn|j>          * <j|c+up|stati>
                cUPDN -=xi1 * cp_i_cdup[dimi*l+stati]* cp_pupmdn_cddn[dimk*l+k] * cp_pupmdn_cddn[dimk*j+k] * cp_i_cdup[dimi*j+stati];
            }

        } // 3,9,11,15,18,22





        if(ndn!=0 && nup != 0)      // tot: 7,10,14,19,21,23
        {
            // tot 7
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_mupmdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_li(beta,cp_i_E[stati], cp_mupmdn_E[j],cp_mdn_E[k],cp_mdn_E[l],w2,w1,w3,PHI_EPS);

// tot=7              -         <stati|c+dn|l>      *  <l|c+up|j>            * <j|cup|k>            * <k|cdn|stati>
                cUPDN -=xi1 * cp_mdn_cddn[diml*stati+l]* cp_mupmdn_cdup[dimj*l+j] * cp_mupmdn_cdup[dimj*k+j] * cp_mdn_cddn[dimk*stati+k];
            }


            // tot 10
            diml = dim_E_mup;
            dimk = dim_E_mupmdn;
            dimj = dim_E_mup;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ki(beta,cp_i_E[stati], cp_mup_E[j],cp_mupmdn_E[k],cp_mup_E[l],w3,w1,w2,PHI_EPS);

// tot=10             +         <stati|c+up|l>      *  <l|c+dn|k>            * <k|cdn|j>            * <j|cup|stati>
                cUPDN +=xi1 * cp_mup_cdup[diml*stati+l]* cp_mupmdn_cddn[dimk*l+k] * cp_mupmdn_cddn[dimk*j+k] * cp_mup_cdup[dimj*stati+j];
            }


            // tot 14
            diml = dim_E_mup;
            dimk = dim_E_mupmdn;
            dimj = dim_E_mdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ki(beta,cp_i_E[stati], cp_mdn_E[j],cp_mupmdn_E[k],cp_mup_E[l],w1,w3,w2,PHI_EPS);

// tot=14             -         <stati|c+up|l>      *  <l|c+dn|k>            * <k|cup|j>          * <j|cdn|stati>
                cUPDN -=xi1 * cp_mup_cdup[diml*stati+l]* cp_mupmdn_cddn[dimk*l+k] * cp_mupmdn_cdup[dimk*j+k] * cp_mdn_cddn[dimj*stati+j];
            }

            // tot 19
            diml = dim_E_mdn;
            dimk = dim_E_mup;
            dimj = dim_E_mupmdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_li(beta,cp_i_E[stati], cp_mupmdn_E[j],cp_mup_E[k],cp_mdn_E[l],w2,w3,w1,PHI_EPS);

// tot=19             +         <stati|c+dn|l>      *  <l|c+up|j>            * <j|cdn|k>            * <k|cup|stati>
                cUPDN +=xi1 * cp_mdn_cddn[diml*stati+l]* cp_mupmdn_cdup[dimj*l+j] * cp_mupmdn_cddn[dimj*k+j] * cp_mup_cdup[dimk*stati+k];
            }


            // tot 21
            diml = dim_E_i;
            dimk = dim_E_mup;
            dimj = dim_E_mdn;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_ji(beta,cp_i_E[stati], cp_mdn_E[j],cp_mup_E[k],cp_i_E[l],w3,w2,w1,PHI_EPS);

// tot=21             -         <stati|c+up|k>      *  <k|cup|l>               * <l|c+dn|j>            * <j|cdn|stati>
                cUPDN -=xi1 * cp_mup_cdup[dimk*stati+k]* cp_mup_cdup[dimk*l+k] * cp_mdn_cddn[dimj*l+j] * cp_mdn_cddn[dimj*stati+j];
            }


            // tot 23
            diml = dim_E_mdn;
            dimk = dim_E_mup;
            dimj = dim_E_i;

            for(int j=0; j < dimj; j++)
            for(int k=0; k < dimk; k++)
            for(int l=0; l < diml; l++)
            {
                xi1 = PhiM_li(beta,cp_i_E[stati], cp_i_E[j],cp_mup_E[k],cp_mdn_E[l],w3,w2,w1,PHI_EPS);

// tot=23             -         <stati|c+dn|l>         *  <l|cdn|j>            * <j|c+up|k>          * <k|cup|stati>
                cUPDN -=xi1 * cp_mdn_cddn[diml*stati+l]* cp_mdn_cddn[diml*j+l] * cp_mup_cdup[dimk*j+k] * cp_mup_cdup[dimk*stati+k];
            }
        }

        pUPDN += cUPDN*boltzZ;
//        cout<<"pUPDN = "<<pUPDN<<endl;
    }
    }// end loop impstates
// ******************************************************************************************************************************
//    cout<<"pDNDN = "<<pDNDN<<endl;
//    cout<<"pUPDN = "<<pUPDN<<endl;
}

//int main(){return 0;};
// vim:foldmethod=marker:tw=0
