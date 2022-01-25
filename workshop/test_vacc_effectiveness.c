/* pomp C snippet file: test_vacc_effectiveness */
/* Time: 2022-01-25 06:01:09.974 -0500 */
/* Salt: 462EFC33CA9C69F616761294 */

#include <pomp.h>
#include <R_ext/Rdynload.h>

 


/* C snippet: 'rinit' */
#define S_0		(__p[__parindex[0]])
#define I_0		(__p[__parindex[1]])
#define q		(__p[__parindex[2]])
#define q_age_1		(__p[__parindex[3]])
#define q_age_2		(__p[__parindex[4]])
#define q_age_3		(__p[__parindex[5]])
#define q_age_4		(__p[__parindex[6]])
#define q_age_5		(__p[__parindex[7]])
#define gamma		(__p[__parindex[8]])
#define sigma		(__p[__parindex[9]])
#define beta1		(__p[__parindex[10]])
#define rho		(__p[__parindex[11]])
#define rho_age_1		(__p[__parindex[12]])
#define rho_age_2		(__p[__parindex[13]])
#define rho_age_3		(__p[__parindex[14]])
#define rho_age_4		(__p[__parindex[15]])
#define rho_age_5		(__p[__parindex[16]])
#define rho_age_u		(__p[__parindex[17]])
#define psi		(__p[__parindex[18]])
#define psi_1		(__p[__parindex[19]])
#define psi_2		(__p[__parindex[20]])
#define psi_3		(__p[__parindex[21]])
#define psi_4		(__p[__parindex[22]])
#define psi_5		(__p[__parindex[23]])
#define psi_u		(__p[__parindex[24]])
#define Cv1		(__p[__parindex[25]])
#define Cv2		(__p[__parindex[26]])
#define Cv3		(__p[__parindex[27]])
#define Cv4		(__p[__parindex[28]])
#define Cv5		(__p[__parindex[29]])
#define Cv6		(__p[__parindex[30]])
#define Cv7		(__p[__parindex[31]])
#define Cv8		(__p[__parindex[32]])
#define Cv9		(__p[__parindex[33]])
#define Cv10		(__p[__parindex[34]])
#define Cv11		(__p[__parindex[35]])
#define Cv12		(__p[__parindex[36]])
#define Cv13		(__p[__parindex[37]])
#define Cv14		(__p[__parindex[38]])
#define Cv15		(__p[__parindex[39]])
#define Cv16		(__p[__parindex[40]])
#define Cv17		(__p[__parindex[41]])
#define Cv18		(__p[__parindex[42]])
#define Cv19		(__p[__parindex[43]])
#define Cv20		(__p[__parindex[44]])
#define Cv21		(__p[__parindex[45]])
#define Cv22		(__p[__parindex[46]])
#define Cv23		(__p[__parindex[47]])
#define Cv24		(__p[__parindex[48]])
#define Cv25		(__p[__parindex[49]])
#define lcv1		(__p[__parindex[50]])
#define lcv2		(__p[__parindex[51]])
#define lcv3		(__p[__parindex[52]])
#define lcv4		(__p[__parindex[53]])
#define lcv5		(__p[__parindex[54]])
#define alpha		(__p[__parindex[55]])
#define dwan		(__p[__parindex[56]])
#define agegroups		(__p[__parindex[57]])
#define N_1		(__covars[__covindex[0]])
#define N_2		(__covars[__covindex[1]])
#define N_3		(__covars[__covindex[2]])
#define N_4		(__covars[__covindex[3]])
#define N_5		(__covars[__covindex[4]])
#define MU_1		(__covars[__covindex[5]])
#define MU_2		(__covars[__covindex[6]])
#define MU_3		(__covars[__covindex[7]])
#define MU_4		(__covars[__covindex[8]])
#define MU_5		(__covars[__covindex[9]])
#define IN_1		(__covars[__covindex[10]])
#define IN_2		(__covars[__covindex[11]])
#define IN_3		(__covars[__covindex[12]])
#define IN_4		(__covars[__covindex[13]])
#define IN_5		(__covars[__covindex[14]])
#define OUT_1		(__covars[__covindex[15]])
#define OUT_2		(__covars[__covindex[16]])
#define OUT_3		(__covars[__covindex[17]])
#define OUT_4		(__covars[__covindex[18]])
#define OUT_5		(__covars[__covindex[19]])
#define p1		(__covars[__covindex[20]])
#define p2		(__covars[__covindex[21]])
#define Births		(__covars[__covindex[22]])
#define eta_a		(__covars[__covindex[23]])
#define Ss_1		(__x[__stateindex[0]])
#define Ss_2		(__x[__stateindex[1]])
#define Ss_3		(__x[__stateindex[2]])
#define Ss_4		(__x[__stateindex[3]])
#define Ss_5		(__x[__stateindex[4]])
#define Sw_1		(__x[__stateindex[5]])
#define Sw_2		(__x[__stateindex[6]])
#define Sw_3		(__x[__stateindex[7]])
#define Sw_4		(__x[__stateindex[8]])
#define Sw_5		(__x[__stateindex[9]])
#define St_1		(__x[__stateindex[10]])
#define St_2		(__x[__stateindex[11]])
#define St_3		(__x[__stateindex[12]])
#define St_4		(__x[__stateindex[13]])
#define St_5		(__x[__stateindex[14]])
#define Es_1		(__x[__stateindex[15]])
#define Es_2		(__x[__stateindex[16]])
#define Es_3		(__x[__stateindex[17]])
#define Es_4		(__x[__stateindex[18]])
#define Es_5		(__x[__stateindex[19]])
#define Ew_1		(__x[__stateindex[20]])
#define Ew_2		(__x[__stateindex[21]])
#define Ew_3		(__x[__stateindex[22]])
#define Ew_4		(__x[__stateindex[23]])
#define Ew_5		(__x[__stateindex[24]])
#define Et_1		(__x[__stateindex[25]])
#define Et_2		(__x[__stateindex[26]])
#define Et_3		(__x[__stateindex[27]])
#define Et_4		(__x[__stateindex[28]])
#define Et_5		(__x[__stateindex[29]])
#define Is_1		(__x[__stateindex[30]])
#define Is_2		(__x[__stateindex[31]])
#define Is_3		(__x[__stateindex[32]])
#define Is_4		(__x[__stateindex[33]])
#define Is_5		(__x[__stateindex[34]])
#define Iw_1		(__x[__stateindex[35]])
#define Iw_2		(__x[__stateindex[36]])
#define Iw_3		(__x[__stateindex[37]])
#define Iw_4		(__x[__stateindex[38]])
#define Iw_5		(__x[__stateindex[39]])
#define It_1		(__x[__stateindex[40]])
#define It_2		(__x[__stateindex[41]])
#define It_3		(__x[__stateindex[42]])
#define It_4		(__x[__stateindex[43]])
#define It_5		(__x[__stateindex[44]])
#define V_1		(__x[__stateindex[45]])
#define V_2		(__x[__stateindex[46]])
#define V_3		(__x[__stateindex[47]])
#define V_4		(__x[__stateindex[48]])
#define V_5		(__x[__stateindex[49]])
#define C_1		(__x[__stateindex[50]])
#define C_2		(__x[__stateindex[51]])
#define C_3		(__x[__stateindex[52]])
#define C_4		(__x[__stateindex[53]])
#define C_5		(__x[__stateindex[54]])
#define Vs_1		(__x[__stateindex[55]])
#define Vs_2		(__x[__stateindex[56]])
#define Vs_3		(__x[__stateindex[57]])
#define Vs_4		(__x[__stateindex[58]])
#define Vs_5		(__x[__stateindex[59]])

void __pomp_rinit (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{
 
  double *Ss = &Ss_1;
  double *Sw = &Sw_1;
  double *St = &St_1;
  double *Es = &Es_1;
  double *Ew = &Ew_1;
  double *Et = &Et_1;
  double *Is = &Is_1;
  double *Iw = &Iw_1;
  double *It = &It_1;
  
  double *V = &V_1;
  double *C = &C_1;
  double *Vs = &Vs_1;
  
  int n = agegroups;
  
  #define N_AGE(J)  (__covars[__covindex[0]  + (J)])
  
  for (int i = 0; i < n; i++) {
    Ss[i]  = nearbyint(S_0*N_AGE(i));
    Sw[i]  = 0;
    St[i]  = Ss[i] + Sw[i];
    
    Es[i] = 0;
    Ew[i] = 0;
    Et[i] = Es[i] + Ew[i];
    
    Is[i] = nearbyint(I_0*N_AGE(i));
    Iw[i] = 0;
    It[i] = Is[i] + Iw[i];
    
    V[i]  = 0;
    C[i]  = 0;
    Vs[i] = 0;
  }
 
}

#undef S_0
#undef I_0
#undef q
#undef q_age_1
#undef q_age_2
#undef q_age_3
#undef q_age_4
#undef q_age_5
#undef gamma
#undef sigma
#undef beta1
#undef rho
#undef rho_age_1
#undef rho_age_2
#undef rho_age_3
#undef rho_age_4
#undef rho_age_5
#undef rho_age_u
#undef psi
#undef psi_1
#undef psi_2
#undef psi_3
#undef psi_4
#undef psi_5
#undef psi_u
#undef Cv1
#undef Cv2
#undef Cv3
#undef Cv4
#undef Cv5
#undef Cv6
#undef Cv7
#undef Cv8
#undef Cv9
#undef Cv10
#undef Cv11
#undef Cv12
#undef Cv13
#undef Cv14
#undef Cv15
#undef Cv16
#undef Cv17
#undef Cv18
#undef Cv19
#undef Cv20
#undef Cv21
#undef Cv22
#undef Cv23
#undef Cv24
#undef Cv25
#undef lcv1
#undef lcv2
#undef lcv3
#undef lcv4
#undef lcv5
#undef alpha
#undef dwan
#undef agegroups
#undef N_1
#undef N_2
#undef N_3
#undef N_4
#undef N_5
#undef MU_1
#undef MU_2
#undef MU_3
#undef MU_4
#undef MU_5
#undef IN_1
#undef IN_2
#undef IN_3
#undef IN_4
#undef IN_5
#undef OUT_1
#undef OUT_2
#undef OUT_3
#undef OUT_4
#undef OUT_5
#undef p1
#undef p2
#undef Births
#undef eta_a
#undef Ss_1
#undef Ss_2
#undef Ss_3
#undef Ss_4
#undef Ss_5
#undef Sw_1
#undef Sw_2
#undef Sw_3
#undef Sw_4
#undef Sw_5
#undef St_1
#undef St_2
#undef St_3
#undef St_4
#undef St_5
#undef Es_1
#undef Es_2
#undef Es_3
#undef Es_4
#undef Es_5
#undef Ew_1
#undef Ew_2
#undef Ew_3
#undef Ew_4
#undef Ew_5
#undef Et_1
#undef Et_2
#undef Et_3
#undef Et_4
#undef Et_5
#undef Is_1
#undef Is_2
#undef Is_3
#undef Is_4
#undef Is_5
#undef Iw_1
#undef Iw_2
#undef Iw_3
#undef Iw_4
#undef Iw_5
#undef It_1
#undef It_2
#undef It_3
#undef It_4
#undef It_5
#undef V_1
#undef V_2
#undef V_3
#undef V_4
#undef V_5
#undef C_1
#undef C_2
#undef C_3
#undef C_4
#undef C_5
#undef Vs_1
#undef Vs_2
#undef Vs_3
#undef Vs_4
#undef Vs_5

/* C snippet: 'skeleton' */
#define S_0		(__p[__parindex[0]])
#define I_0		(__p[__parindex[1]])
#define q		(__p[__parindex[2]])
#define q_age_1		(__p[__parindex[3]])
#define q_age_2		(__p[__parindex[4]])
#define q_age_3		(__p[__parindex[5]])
#define q_age_4		(__p[__parindex[6]])
#define q_age_5		(__p[__parindex[7]])
#define gamma		(__p[__parindex[8]])
#define sigma		(__p[__parindex[9]])
#define beta1		(__p[__parindex[10]])
#define rho		(__p[__parindex[11]])
#define rho_age_1		(__p[__parindex[12]])
#define rho_age_2		(__p[__parindex[13]])
#define rho_age_3		(__p[__parindex[14]])
#define rho_age_4		(__p[__parindex[15]])
#define rho_age_5		(__p[__parindex[16]])
#define rho_age_u		(__p[__parindex[17]])
#define psi		(__p[__parindex[18]])
#define psi_1		(__p[__parindex[19]])
#define psi_2		(__p[__parindex[20]])
#define psi_3		(__p[__parindex[21]])
#define psi_4		(__p[__parindex[22]])
#define psi_5		(__p[__parindex[23]])
#define psi_u		(__p[__parindex[24]])
#define Cv1		(__p[__parindex[25]])
#define Cv2		(__p[__parindex[26]])
#define Cv3		(__p[__parindex[27]])
#define Cv4		(__p[__parindex[28]])
#define Cv5		(__p[__parindex[29]])
#define Cv6		(__p[__parindex[30]])
#define Cv7		(__p[__parindex[31]])
#define Cv8		(__p[__parindex[32]])
#define Cv9		(__p[__parindex[33]])
#define Cv10		(__p[__parindex[34]])
#define Cv11		(__p[__parindex[35]])
#define Cv12		(__p[__parindex[36]])
#define Cv13		(__p[__parindex[37]])
#define Cv14		(__p[__parindex[38]])
#define Cv15		(__p[__parindex[39]])
#define Cv16		(__p[__parindex[40]])
#define Cv17		(__p[__parindex[41]])
#define Cv18		(__p[__parindex[42]])
#define Cv19		(__p[__parindex[43]])
#define Cv20		(__p[__parindex[44]])
#define Cv21		(__p[__parindex[45]])
#define Cv22		(__p[__parindex[46]])
#define Cv23		(__p[__parindex[47]])
#define Cv24		(__p[__parindex[48]])
#define Cv25		(__p[__parindex[49]])
#define lcv1		(__p[__parindex[50]])
#define lcv2		(__p[__parindex[51]])
#define lcv3		(__p[__parindex[52]])
#define lcv4		(__p[__parindex[53]])
#define lcv5		(__p[__parindex[54]])
#define alpha		(__p[__parindex[55]])
#define dwan		(__p[__parindex[56]])
#define agegroups		(__p[__parindex[57]])
#define N_1		(__covars[__covindex[0]])
#define N_2		(__covars[__covindex[1]])
#define N_3		(__covars[__covindex[2]])
#define N_4		(__covars[__covindex[3]])
#define N_5		(__covars[__covindex[4]])
#define MU_1		(__covars[__covindex[5]])
#define MU_2		(__covars[__covindex[6]])
#define MU_3		(__covars[__covindex[7]])
#define MU_4		(__covars[__covindex[8]])
#define MU_5		(__covars[__covindex[9]])
#define IN_1		(__covars[__covindex[10]])
#define IN_2		(__covars[__covindex[11]])
#define IN_3		(__covars[__covindex[12]])
#define IN_4		(__covars[__covindex[13]])
#define IN_5		(__covars[__covindex[14]])
#define OUT_1		(__covars[__covindex[15]])
#define OUT_2		(__covars[__covindex[16]])
#define OUT_3		(__covars[__covindex[17]])
#define OUT_4		(__covars[__covindex[18]])
#define OUT_5		(__covars[__covindex[19]])
#define p1		(__covars[__covindex[20]])
#define p2		(__covars[__covindex[21]])
#define Births		(__covars[__covindex[22]])
#define eta_a		(__covars[__covindex[23]])
#define Ss_1		(__x[__stateindex[0]])
#define Ss_2		(__x[__stateindex[1]])
#define Ss_3		(__x[__stateindex[2]])
#define Ss_4		(__x[__stateindex[3]])
#define Ss_5		(__x[__stateindex[4]])
#define Sw_1		(__x[__stateindex[5]])
#define Sw_2		(__x[__stateindex[6]])
#define Sw_3		(__x[__stateindex[7]])
#define Sw_4		(__x[__stateindex[8]])
#define Sw_5		(__x[__stateindex[9]])
#define St_1		(__x[__stateindex[10]])
#define St_2		(__x[__stateindex[11]])
#define St_3		(__x[__stateindex[12]])
#define St_4		(__x[__stateindex[13]])
#define St_5		(__x[__stateindex[14]])
#define Es_1		(__x[__stateindex[15]])
#define Es_2		(__x[__stateindex[16]])
#define Es_3		(__x[__stateindex[17]])
#define Es_4		(__x[__stateindex[18]])
#define Es_5		(__x[__stateindex[19]])
#define Ew_1		(__x[__stateindex[20]])
#define Ew_2		(__x[__stateindex[21]])
#define Ew_3		(__x[__stateindex[22]])
#define Ew_4		(__x[__stateindex[23]])
#define Ew_5		(__x[__stateindex[24]])
#define Et_1		(__x[__stateindex[25]])
#define Et_2		(__x[__stateindex[26]])
#define Et_3		(__x[__stateindex[27]])
#define Et_4		(__x[__stateindex[28]])
#define Et_5		(__x[__stateindex[29]])
#define Is_1		(__x[__stateindex[30]])
#define Is_2		(__x[__stateindex[31]])
#define Is_3		(__x[__stateindex[32]])
#define Is_4		(__x[__stateindex[33]])
#define Is_5		(__x[__stateindex[34]])
#define Iw_1		(__x[__stateindex[35]])
#define Iw_2		(__x[__stateindex[36]])
#define Iw_3		(__x[__stateindex[37]])
#define Iw_4		(__x[__stateindex[38]])
#define Iw_5		(__x[__stateindex[39]])
#define It_1		(__x[__stateindex[40]])
#define It_2		(__x[__stateindex[41]])
#define It_3		(__x[__stateindex[42]])
#define It_4		(__x[__stateindex[43]])
#define It_5		(__x[__stateindex[44]])
#define V_1		(__x[__stateindex[45]])
#define V_2		(__x[__stateindex[46]])
#define V_3		(__x[__stateindex[47]])
#define V_4		(__x[__stateindex[48]])
#define V_5		(__x[__stateindex[49]])
#define C_1		(__x[__stateindex[50]])
#define C_2		(__x[__stateindex[51]])
#define C_3		(__x[__stateindex[52]])
#define C_4		(__x[__stateindex[53]])
#define C_5		(__x[__stateindex[54]])
#define Vs_1		(__x[__stateindex[55]])
#define Vs_2		(__x[__stateindex[56]])
#define Vs_3		(__x[__stateindex[57]])
#define Vs_4		(__x[__stateindex[58]])
#define Vs_5		(__x[__stateindex[59]])
#define DSs_1		(__f[__stateindex[0]])
#define DSs_2		(__f[__stateindex[1]])
#define DSs_3		(__f[__stateindex[2]])
#define DSs_4		(__f[__stateindex[3]])
#define DSs_5		(__f[__stateindex[4]])
#define DSw_1		(__f[__stateindex[5]])
#define DSw_2		(__f[__stateindex[6]])
#define DSw_3		(__f[__stateindex[7]])
#define DSw_4		(__f[__stateindex[8]])
#define DSw_5		(__f[__stateindex[9]])
#define DSt_1		(__f[__stateindex[10]])
#define DSt_2		(__f[__stateindex[11]])
#define DSt_3		(__f[__stateindex[12]])
#define DSt_4		(__f[__stateindex[13]])
#define DSt_5		(__f[__stateindex[14]])
#define DEs_1		(__f[__stateindex[15]])
#define DEs_2		(__f[__stateindex[16]])
#define DEs_3		(__f[__stateindex[17]])
#define DEs_4		(__f[__stateindex[18]])
#define DEs_5		(__f[__stateindex[19]])
#define DEw_1		(__f[__stateindex[20]])
#define DEw_2		(__f[__stateindex[21]])
#define DEw_3		(__f[__stateindex[22]])
#define DEw_4		(__f[__stateindex[23]])
#define DEw_5		(__f[__stateindex[24]])
#define DEt_1		(__f[__stateindex[25]])
#define DEt_2		(__f[__stateindex[26]])
#define DEt_3		(__f[__stateindex[27]])
#define DEt_4		(__f[__stateindex[28]])
#define DEt_5		(__f[__stateindex[29]])
#define DIs_1		(__f[__stateindex[30]])
#define DIs_2		(__f[__stateindex[31]])
#define DIs_3		(__f[__stateindex[32]])
#define DIs_4		(__f[__stateindex[33]])
#define DIs_5		(__f[__stateindex[34]])
#define DIw_1		(__f[__stateindex[35]])
#define DIw_2		(__f[__stateindex[36]])
#define DIw_3		(__f[__stateindex[37]])
#define DIw_4		(__f[__stateindex[38]])
#define DIw_5		(__f[__stateindex[39]])
#define DIt_1		(__f[__stateindex[40]])
#define DIt_2		(__f[__stateindex[41]])
#define DIt_3		(__f[__stateindex[42]])
#define DIt_4		(__f[__stateindex[43]])
#define DIt_5		(__f[__stateindex[44]])
#define DV_1		(__f[__stateindex[45]])
#define DV_2		(__f[__stateindex[46]])
#define DV_3		(__f[__stateindex[47]])
#define DV_4		(__f[__stateindex[48]])
#define DV_5		(__f[__stateindex[49]])
#define DC_1		(__f[__stateindex[50]])
#define DC_2		(__f[__stateindex[51]])
#define DC_3		(__f[__stateindex[52]])
#define DC_4		(__f[__stateindex[53]])
#define DC_5		(__f[__stateindex[54]])
#define DVs_1		(__f[__stateindex[55]])
#define DVs_2		(__f[__stateindex[56]])
#define DVs_3		(__f[__stateindex[57]])
#define DVs_4		(__f[__stateindex[58]])
#define DVs_5		(__f[__stateindex[59]])

void __pomp_skelfn (double *__f, const double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
  
    // define  the rate of waning of vaccine derived immunity
    double delta = 1/dwan; 

    int i,j,k,l;
    int n = agegroups; 
    
    double lambda1;
    double lambda2;
    
    const double *Cvlocal = &Cv1;
    const double *lvlocal = &lcv1;
    const double *qvlocal = &q_age_1;
    
    const double *Sslocal = &Ss_1;
    const double *Swlocal = Sslocal+n;
    const double *Stlocal = Swlocal+n;
    const double *Eslocal = Stlocal+n;
    const double *Ewlocal = Eslocal+n;
    const double *Etlocal = Ewlocal+n;
    const double *Islocal = Etlocal+n;
    const double *Iwlocal = Islocal+n;
    const double *Itlocal = Iwlocal+n;
    const double *Vlocal  = Itlocal+n;
    const double *Clocal  = Vlocal+n;
    const double *Vslocal = Clocal+n; 
    
    
    double *DSslocal = &DSs_1;
    double *DSwlocal = DSslocal+n;
    double *DStlocal = DSwlocal+n;
    double *DEslocal = DStlocal+n;
    double *DEwlocal = DEslocal+n;
    double *DEtlocal = DEwlocal+n;
    double *DIslocal = DEtlocal+n;
    double *DIwlocal = DIslocal+n;
    double *DItlocal = DIwlocal+n;
    double *DVlocal  = DItlocal+n;
    double *DClocal  = DVlocal+n;
    double *DVslocal = DClocal+n; 
    
    #define CM(J,K) Cvlocal[(J)+n*(K)]
    #define LV(K) lvlocal[(K)]
    #define QAGE(K) qvlocal[(K)]
    
    // Special Age-specific Covariates macro definitions  
    
    //#define N_AGE(J)   (__covars[__covindex[0]  + (J)])
    #define MU_AGE(J)  (__covars[__covindex[5]  + (J)])
    
    #define Ss(K) Sslocal[(K)]
    #define Sw(K) Swlocal[(K)]
    #define St(K) Stlocal[(K)]
    #define Es(K) Eslocal[(K)]
    #define Ew(K) Ewlocal[(K)]
    #define Et(K) Etlocal[(K)]
    #define Is(K) Islocal[(K)]
    #define Iw(K) Iwlocal[(K)]
    #define It(K) Itlocal[(K)]
    #define V(K) Vlocal[(K)]
    #define C(K) Clocal[(K)]
    #define Vs(K) Vslocal[(K)]
    
    #define DSs(K) DSslocal[(K)]
    #define DSw(K) DSwlocal[(K)]
    #define DSt(K) DStlocal[(K)]
    #define DEs(K) DEslocal[(K)]
    #define DEw(K) DEwlocal[(K)]
    #define DEt(K) DEtlocal[(K)]
    #define DIs(K) DIslocal[(K)]
    #define DIw(K) DIwlocal[(K)]
    #define DIt(K) DItlocal[(K)]
    #define DV(K)  DVlocal[(K)]
    #define DC(K)  DClocal[(K)]
    #define DVs(K) DVslocal[(K)]
    
    
    
    /* ============================================================================== */
    /* augment the contact matrix to incorporate seasonality in the [5,15) sge cohort */
    /* ============================================================================== */      
    // declare and initialize
    double augCM[n][n];
    
    // modify to incorporate seasonality  
    for(k=0; k<n; k++) {
      for(l=0; l<n; l++) {
        if(k == 1 && l == 1) {
          augCM[k][l] = CM(k,l)*(1-beta1*sin(2*M_PI*t));
        } else {
          augCM[k][l] = CM(k,l);
        }
      }
    }
  
    
    /* ============================================================================= */
    /* ========================== defining epidemic process ======================== */
    /* ============================ with type 1 mortality ========================== */
    /* ============================================================================= */
    
    for (i = 0; i < n; i++) {
  
      lambda1 = 0.0;
      
      if(i == 0) {
        
        // births occur only in the first age class - specially treat the [0,5) age cohort
        
        for (j = 0; j < n; j++) {
          lambda1 += q*QAGE(i)*augCM[i][j]*(It(j))/N_AGE(j);
          }
            
        /* ============================ */
        /* balance transition equations */
        /* ============================ */
        
        // vaccinated and unvaccinated births - calculated using covariate population  
        
        double uv_births = (1-(1-alpha)*p1)*Births;
        double  v_births = (1-alpha)*p1*Births;
        
        DSs(i) = uv_births  - (lambda1 + LV(i) - MU_AGE(i))*Ss(i);
        DSw(i) = delta*V(i) - (lambda1 + LV(i) - MU_AGE(i))*Sw(i);
        DSt(i) = Ss(i) + Sw(i);  
        
        DEs(i) = lambda1*Ss(i) - (sigma + LV(i) - MU_AGE(i))*Es(i);   
        DEw(i) = lambda1*Sw(i) - (sigma + LV(i) - MU_AGE(i))*Ew(i);
        DEt(i) = Es(i) + Ew(i); 
        
        
        DIs(i) = sigma*Es(i) - (gamma + LV(i) - MU_AGE(i))*Is(i);
        DIw(i) = sigma*Ew(i) - (gamma + LV(i) - MU_AGE(i))*Iw(i);
        DIt(i) = Is(i) + Iw(i);
        
        DV(i)  = v_births - (delta + LV(i) - MU_AGE(i))*V(i);
        
        DVs(i) = delta*V(i); 
      
      } else if (i == 1) {
        
          // individuals age in and out of age cohort after the first
        
          for (j = 0; j < n; j++) {
            lambda1 += q*QAGE(i)*augCM[i][j]*(It(j))/N_AGE(j);
          }
            
          /* ============================ */
          /* balance tranistion equations */
          /* ============================ */
          
          double uv_grads_Ss = LV(i-1)*(1-(1-alpha)*p2)*Ss(i-1);
          double uv_grads_Sw = LV(i-1)*(1-(1-alpha)*p2)*Sw(i-1);
          double     v_grads = LV(i-1)*(1-alpha)*p2*St(i-1);
          
          DSs(i) =  uv_grads_Ss - (lambda1 + LV(i) - MU_AGE(i))*Ss(i);               
          DSw(i) =  uv_grads_Sw + delta*V(i) - (lambda1 + LV(i) - MU_AGE(i))*Sw(i);               
          DSt(i) =  Ss(i) + Sw(i); 
          
          DEs(i) = LV(i-1)*Es(i-1) + lambda1*Ss(i) - (sigma + LV(i) - MU_AGE(i))*Es(i); 
          DEw(i) = LV(i-1)*Ew(i-1) + lambda1*Sw(i) - (sigma + LV(i) - MU_AGE(i))*Ew(i); 
          DEt(i) = Es(i) + Ew(i); 
          
          DIs(i) = LV(i-1)*Is(i-1) + sigma*Es(i) - (gamma + LV(i) - MU_AGE(i))*Is(i);
          DIw(i) = LV(i-1)*Iw(i-1) + sigma*Ew(i) - (gamma + LV(i) - MU_AGE(i))*Iw(i);
          DIt(i) = Is(i) + Iw(i);
          
          DV(i)  = v_grads + LV(i-1)*V(i-1) - (delta + LV(i) - MU_AGE(i))*V(i);
          
          DVs(i) = delta*V(i);
      
      } else {
        
          // individuals age in and out of age cohort after the first
          for (j = 0; j < n; j++) {
            lambda1 += q*QAGE(i)*augCM[i][j]*(It(j))/N_AGE(j);
          }
          
          /* ============================ */
          /* balance transition equations */
          /* ============================ */
          
          DSs(i) = LV(i-1)*Ss(i-1) - (lambda1 + LV(i) - MU_AGE(i))*Ss(i);
          DSw(i) = LV(i-1)*Sw(i-1) + delta*V(i) - (lambda1 + LV(i) - MU_AGE(i))*Sw(i);
          DSt(i) = Ss(i) + Sw(i);
          
          DEs(i) = LV(i-1)*Es(i-1) + lambda1*Ss(i) - (sigma + LV(i) - MU_AGE(i))*Es(i); 
          DEw(i) = LV(i-1)*Ew(i-1) + lambda1*Sw(i) - (sigma + LV(i) - MU_AGE(i))*Ew(i); 
          DEt(i) = Es(i) + Ew(i); 
          
          DIs(i) = LV(i-1)*Is(i-1) + sigma*Es(i) - (gamma + LV(i) - MU_AGE(i))*Is(i);
          DIw(i) = LV(i-1)*Iw(i-1) + sigma*Ew(i) - (gamma + LV(i) - MU_AGE(i))*Iw(i);
          DIt(i) = Is(i) + Iw(i); 
          
          DV(i)  = LV(i-1)*V(i-1)  - (delta + LV(i) - MU_AGE(i))*V(i);
          
          DVs(i) = delta*V(i);
          
      }
      
      /* ========================= */
      /* Keep account of new cases */
      /* ========================= */
      DC(i) = gamma*It(i); 
    }
 
}

#undef S_0
#undef I_0
#undef q
#undef q_age_1
#undef q_age_2
#undef q_age_3
#undef q_age_4
#undef q_age_5
#undef gamma
#undef sigma
#undef beta1
#undef rho
#undef rho_age_1
#undef rho_age_2
#undef rho_age_3
#undef rho_age_4
#undef rho_age_5
#undef rho_age_u
#undef psi
#undef psi_1
#undef psi_2
#undef psi_3
#undef psi_4
#undef psi_5
#undef psi_u
#undef Cv1
#undef Cv2
#undef Cv3
#undef Cv4
#undef Cv5
#undef Cv6
#undef Cv7
#undef Cv8
#undef Cv9
#undef Cv10
#undef Cv11
#undef Cv12
#undef Cv13
#undef Cv14
#undef Cv15
#undef Cv16
#undef Cv17
#undef Cv18
#undef Cv19
#undef Cv20
#undef Cv21
#undef Cv22
#undef Cv23
#undef Cv24
#undef Cv25
#undef lcv1
#undef lcv2
#undef lcv3
#undef lcv4
#undef lcv5
#undef alpha
#undef dwan
#undef agegroups
#undef N_1
#undef N_2
#undef N_3
#undef N_4
#undef N_5
#undef MU_1
#undef MU_2
#undef MU_3
#undef MU_4
#undef MU_5
#undef IN_1
#undef IN_2
#undef IN_3
#undef IN_4
#undef IN_5
#undef OUT_1
#undef OUT_2
#undef OUT_3
#undef OUT_4
#undef OUT_5
#undef p1
#undef p2
#undef Births
#undef eta_a
#undef Ss_1
#undef Ss_2
#undef Ss_3
#undef Ss_4
#undef Ss_5
#undef Sw_1
#undef Sw_2
#undef Sw_3
#undef Sw_4
#undef Sw_5
#undef St_1
#undef St_2
#undef St_3
#undef St_4
#undef St_5
#undef Es_1
#undef Es_2
#undef Es_3
#undef Es_4
#undef Es_5
#undef Ew_1
#undef Ew_2
#undef Ew_3
#undef Ew_4
#undef Ew_5
#undef Et_1
#undef Et_2
#undef Et_3
#undef Et_4
#undef Et_5
#undef Is_1
#undef Is_2
#undef Is_3
#undef Is_4
#undef Is_5
#undef Iw_1
#undef Iw_2
#undef Iw_3
#undef Iw_4
#undef Iw_5
#undef It_1
#undef It_2
#undef It_3
#undef It_4
#undef It_5
#undef V_1
#undef V_2
#undef V_3
#undef V_4
#undef V_5
#undef C_1
#undef C_2
#undef C_3
#undef C_4
#undef C_5
#undef Vs_1
#undef Vs_2
#undef Vs_3
#undef Vs_4
#undef Vs_5
#undef DSs_1
#undef DSs_2
#undef DSs_3
#undef DSs_4
#undef DSs_5
#undef DSw_1
#undef DSw_2
#undef DSw_3
#undef DSw_4
#undef DSw_5
#undef DSt_1
#undef DSt_2
#undef DSt_3
#undef DSt_4
#undef DSt_5
#undef DEs_1
#undef DEs_2
#undef DEs_3
#undef DEs_4
#undef DEs_5
#undef DEw_1
#undef DEw_2
#undef DEw_3
#undef DEw_4
#undef DEw_5
#undef DEt_1
#undef DEt_2
#undef DEt_3
#undef DEt_4
#undef DEt_5
#undef DIs_1
#undef DIs_2
#undef DIs_3
#undef DIs_4
#undef DIs_5
#undef DIw_1
#undef DIw_2
#undef DIw_3
#undef DIw_4
#undef DIw_5
#undef DIt_1
#undef DIt_2
#undef DIt_3
#undef DIt_4
#undef DIt_5
#undef DV_1
#undef DV_2
#undef DV_3
#undef DV_4
#undef DV_5
#undef DC_1
#undef DC_2
#undef DC_3
#undef DC_4
#undef DC_5
#undef DVs_1
#undef DVs_2
#undef DVs_3
#undef DVs_4
#undef DVs_5

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {++__pomp_load_stack;}

void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}

void R_init_test_vacc_effectiveness (DllInfo *info)
{
R_RegisterCCallable("test_vacc_effectiveness", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("test_vacc_effectiveness", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("test_vacc_effectiveness", "__pomp_rinit", (DL_FUNC) __pomp_rinit);
R_RegisterCCallable("test_vacc_effectiveness", "__pomp_skelfn", (DL_FUNC) __pomp_skelfn);
}
