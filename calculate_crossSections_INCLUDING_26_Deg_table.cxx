
#include <TGraphErrors.h>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph2DErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include <fstream>
#include <iostream>
#include <TSystem.h>
#include <vector>
#include <cmath>
#include <algorithm>

using std::vector;
using std::ifstream;
using namespace std; 

//GLOBAL VARIABLES
const float  N_A            = 6.02214*1e+23;         // number of particles per one gram
const float  Q_E            = 1.60217663*1e-19;         //1C = 6.24e18 electrons
const int    phi_hi         =  100;
const int    phi_lo         = -100;
const int    theta_hi       =  65; //100; dave uses 65
const int    theta_lo       = -65; //-100; dave uses -65
const double Einitial       = 10.551;
const double Mp             = 0.9382731;
const float  density_LT     = 0.16743  ;       //Has g/cm3 units 
const double length_LT      = 10.0;          //Has cm units  
const double sim_charge_ST  = 1.;            //Has uC units 
const double sim_charge_LT  = 1.;            //Has uC units 
const double AM_LT          = 2.0;
const int    N_BINS         = 16;
const double TOLERANCE      = 1e-8;
const bool   F2_option      = false;        // it should be false always, calculating F2 event by event is not the right approach
const char* var_dep         = "delta";      // delta or xb
const char* trigger_select  = "ELREAL";     // ELREAL or anything else to include 3/4 triggers (not reason to include those in physics analysis)

struct Variables {
  // Variables for SIMULATION
  /*  float hsdp;
  float hspxp;
  float hspyp;
  float stop_id;
  float hpsxfp;    // SIM focal plane x
  float hpsyfp;    // SIM focal plane y
  float hpsxptar;  // SIM target xp      
  float hpsyptar;  // SIM target yp      
  float hpsytar;   // SIM target y       
  float hpsyptari;
  float hpsxptari;
  float hsdpi;*/

  double hsdp;
  double hspxp;
  double hspyp;
  double stop_id;
  double hpsxfp;    // SIM focal plane x                                                                                    
  double hpsyfp;    // SIM focal plane y                                                                
  double hpsxptar;  // SIM target xp                                                                                                                            
  double hpsyptar;  // SIM target yp                                                                                                                            
  double hpsytar;   // SIM target y                                                              
  double hpsyptari;
  double hpsxptari;
  double hsdpi;




  // Variables for ACTUAL ROOTFILE (Actual Data)
  double ptardp;   // data delta         
  double ptarth;   // data target theta  
  double ptarph;   // data target phi    
  double ptarx;    // data target x      
  double ptary;    // data target y      
  double dpsxfp;   // data focal plane x 
  double dpsyfp;   // data focal plane y 
  double dpsxptar; // data focal plane xp
  double dpsyptar; // data focal plane yp
  double dW;       // data W     
  double dQ2;      // data Q2    
  double dnu;      // data omega 
  double dtheta;   // data theta         
  double decal;
  double npeSum;
  double AvgCurr;
  double CurrentFlag;
};

struct ScaleFactors {
  int    run_number;
  double charge;
  int    PS;
  double Eff_Fid ;
  double Eff_Elec_LiveTime;
  double Eff_Comp_LiveTime;
  double trigger ; 
  double cal_eff;
  double cer_eff;
};

void fillScaleFactors(const TString& table_file, 
                      const TString& target_input, 
                      const TString& angle_input, 
                      const TString& spectrometer, 
                      const TString& momentum_spec_input, 
                      double p_spec,
                      std::vector<ScaleFactors>& run_data_vector);

void ImportRadcor(vector<double> &v1,vector<double> &v2,vector<double> &v3,vector<double> &v4,vector<double> &v5,vector<double> &v6,vector<double> &v7,vector<double> &v8,vector<double> &v9, vector<double> &v10,vector<double> &v11,vector<double> &v12,vector<double> &v13, const char * filename);
void setBranchAddresses_mc(TTree *tree, Variables &vars, TString spectrometer_option); 
void setBranchAddresses_data(TTree *tree, Variables &vars, TString spectrometer_option) ;
double calculate_mc_scale_factor(int nentries_mc , double p_spec, double density , double length, double AM , double sim_charge , int delta_hi, int delta_lo , TString targType ) ;
double calculate_BjorkenX ( double p_spec , double hsdp , double hpsyptar , double hpsxptar , double& Eprime , double& thetadeg, double theta_central , TString spectrometer_option);
TGraph2DErrors* createGraph2D( const std::vector<double>& V2, const std::vector<double>& V3, const std::vector<double>& VData, int size);
void printHistogramContents(TH1F *hist);
void PlotWithRatio(TCanvas *&canvas, TH1F *Data_ST, TH1F *Data_LT, TH1F *MC_ST, TH1F *MC_LT, TH1F *MC_ST_UNWeigthed,TH1F *MC_LT_UNWeigthed , TString out_name);
double deltaToX(double delta , double p_spec , double Einitial , double angle);
double XTodelta(double X_Bjorken, double p_spec, double Einitial, double angle);
void convertHistogramToGraph(TH1F* h_delta, TGraphErrors*& graph,  double p_spec , double Einitial , double angle , TString out_name);
void printGraphContents(TGraphErrors *graph);
void writeHistogramToTextFile(TH1F* hist, double p_spec , double Einitial , double angle, const char* filename) ;
void scaleGraph(TGraph* graph, double scaleFactor);
void plottingWeightFactors( std::vector<double> hsdp_values_ST , std::vector<double> weight_values_ST, std::vector<double> hsdp_values_LT , std::vector<double> weight_values_LT, TString out_name);
void plottingWeightFactors_Born( std::vector<double> hsdp_values_ST , std::vector<double> weight_values_ST, std::vector<double> hsdp_values_LT , std::vector<double> weight_values_LT, TString out_name);
void mergePDFs( TString out_name );
double linearInterpolation(double x0, double y0, double x1, double y1, double x) ;
bool areEqual(double a, double b, double tolerance);
double Interpolate1DFromGraph2DErrors(TGraph2DErrors* graph, double Eprime, double angle);
double findNearestValue(double input, double angle); 
double weight_calculateF2(double theta_deg , double Eprime , double Eb);


double weight_CSB_corr(TString targType, TString spectrometer_option, double angle, double Ep) {
  double a, b, R, finalFactor;

  // Define a tolerance for angle comparison
  const double tolerance = 0.5;

  if (spectrometer_option == "HMS") {

    if (targType == "He03" && fabs(angle - 20) <= tolerance) { a = 2.770; b = -2.543; }
    if (targType == "He04" && fabs(angle - 20) <= tolerance) { a = 2.573; b = -2.377; }
    if (targType == "Li06" && fabs(angle - 20) <= tolerance) { a = 1.973; b = -2.353; }
    if (targType == "Li07" && fabs(angle - 20) <= tolerance) { a = 2.118; b = -2.434; }
    if (targType == "Be09" && fabs(angle - 20) <= tolerance) { a = 2.207; b = -2.319; }
    if (targType == "BC10" && fabs(angle - 20) <= tolerance) { a = 1.942; b = -2.269; }
    if (targType == "BC11" && fabs(angle - 20) <= tolerance) { a = 1.971; b = -2.277; }
    if (targType == "Al27" && fabs(angle - 20) <= tolerance) { a = 1.949; b = -2.225; }
    if (targType == "Fe54" && fabs(angle - 20) <= tolerance) { a = 1.877; b = -2.136; }
    if (targType == "Ni58" && fabs(angle - 20) <= tolerance) { a = 1.679; b = -2.146; }
    if (targType == "Cu63" && fabs(angle - 20) <= tolerance) { a = 2.453; b = -2.065; }
    if (targType == "Ni64" && fabs(angle - 20) <= tolerance) { a = 1.652; b = -2.136; }
    if (targType == "Ag10" && fabs(angle - 20) <= tolerance) { a = 2.215; b = -2.060; }
    if (targType == "Sn11" && fabs(angle - 20) <= tolerance) { a = 1.195; b = -1.671; }
    if (targType == "Au19" && fabs(angle - 20) <= tolerance) { a = 2.009; b = -2.051; }
    if (targType == "Th23" && fabs(angle - 20) <= tolerance) { a = 2.143; b = -2.045; }
    if (targType == "Ca12" && fabs(angle - 20) <= tolerance) { a = 1.907; b = -2.242; }
    if (targType == "Ca12" && fabs(angle - 26) <= tolerance) { a = 2.628; b = -2.820; }
    if (targType == "Ca12" && fabs(angle - 35) <= tolerance) { a = 3.185; b = -3.452; }
    if (targType == "Ca40" && fabs(angle - 20) <= tolerance) { a = 2.378; b = -2.154; }
    if (targType == "Ca40" && fabs(angle - 26) <= tolerance) { a = 2.833; b = -2.535; }
    if (targType == "Ca40" && fabs(angle - 35) <= tolerance) { a = 3.443; b = -3.054; }
    if (targType == "Ca48" && fabs(angle - 20) <= tolerance) { a = 2.450; b = -2.150; }
    if (targType == "Ca48" && fabs(angle - 35) <= tolerance) { a = 3.391; b = -2.949; }
    if (targType == "LD2" && fabs(angle - 20) <= tolerance) { a = 2.983; b = -2.508; }
    if (targType == "LD2" && fabs(angle - 26) <= tolerance) { a = 3.711; b = -3.221; }
    if (targType == "LD2" && fabs(angle - 35) <= tolerance) { a = 4.257; b = -4.120; }
    if (targType == "Dummy" && fabs(angle - 20) <= tolerance) { a = 1.768; b = -2.183; }
    if (targType == "Dummy" && fabs(angle - 26) <= tolerance) { a = 2.480; b = -2.739; }
    if (targType == "Dummy" && fabs(angle - 35) <= tolerance) { a = 3.104; b = -3.353; }
  }

  R = TMath::Exp(a + b * Ep);
  finalFactor = (1.0 / (1.0 + R));
  std::cout << targType << " , " << spectrometer_option << " , " << angle
	    << " CSB! Ep:" << Ep << " and FinalFactor is : " << finalFactor
	    << " , and R is " << R << " , and a and b are: " << a << " , " << b << std::endl;
  return finalFactor;
}





/*
double weight_CSB_corr(   TString targType ,  TString spectrometer_option , double angle , double Ep ){
  double a, b , R, finalFactor; 

  if ( spectrometer_option=="HMS"  ) {

          if (targType=="He03"  && angle==20) { a = 2.770 ;  b = -2.543 ;   }
          if (targType=="He04"  && angle==20) { a = 2.573 ;  b = -2.377 ;   }
	  if (targType=="Li06"  && angle==20) { a = 1.973 ;  b = -2.353 ;   }
	  if (targType=="Li07"  && angle==20) { a = 2.118 ;  b = -2.434 ;   }
	  if (targType=="Be09"  && angle==20) { a = 2.207 ;  b = -2.319 ;   }
	  if (targType=="BC10"  && angle==20) { a = 1.942 ;  b = -2.269 ;   }
	  if (targType=="BC11"  && angle==20) { a = 1.971 ;  b = -2.277 ;   }
	  if (targType=="Al27"  && angle==20) { a = 1.949 ;  b = -2.225 ;   }
	  if (targType=="Fe54"  && angle==20) { a = 1.877 ;  b = -2.136 ;   }
	  if (targType=="Ni58"  && angle==20) { a = 1.679 ;  b = -2.146 ;   }
	  if (targType=="Cu63"  && angle==20) { a = 2.453 ;  b = -2.065 ;   }
	  if (targType=="Ni64"  && angle==20) { a = 1.652 ;  b = -2.136 ;   }
	  if (targType=="Ag10"  && angle==20) { a = 2.215 ;  b = -2.060 ;   }
	  if (targType=="Sn11"  && angle==20) { a = 1.195 ;  b = -1.671 ;   }
	  if (targType=="Au19"  && angle==20) { a = 2.009 ;  b = -2.051 ;   }
	  if (targType=="Th23"  && angle==20) { a = 2.143 ;  b = -2.045 ;   }
          if (targType=="Ca12"  && angle==20) { a = 1.907 ;  b = -2.242 ;   }
	  if (targType=="Ca12"  && angle==26) { a = 2.628 ;  b = -2.820 ;   }
	  if (targType=="Ca12"  && angle==35) { a = 3.185 ;  b = -3.452 ;   }
          if (targType=="Ca40"  && angle==20) { a = 2.378 ;  b = -2.154 ;   }
          if (targType=="Ca40"  && angle==26) { a = 2.833 ;  b = -2.535 ;   }
          if (targType=="Ca40"  && angle==35) { a = 3.443 ;  b = -3.054 ;   }
          if (targType=="Ca48"  && angle==20) { a = 2.450 ;  b = -2.150 ;   }
          if (targType=="Ca48"  && angle==35) { a = 3.391 ;  b = -2.949 ;   }
	  if (targType=="LD2"   && angle==20) { a = 2.983 ;  b = -2.508 ;   }
	  if (targType=="LD2"   && angle==26) { a = 3.711 ;  b = -3.221 ;   }
	  if (targType=="LD2"   && angle==35) { a = 4.257 ;  b = -4.120 ;   }
	  if (targType=="Dummy" && angle==20) { a = 1.768 ;  b = -2.183 ;   }
	  if (targType=="Dummy" && angle==26) { a = 2.480 ;  b = -2.739 ;   }
	  if (targType=="Dummy" && angle==35) { a = 3.104 ;  b = -3.353 ;   }
  }
  R = TMath::Exp(a+b*Ep);
  finalFactor=(1.0/(1.0+R)) ; 
  cout<< targType << " , " << spectrometer_option << " , " << angle << " CSB! Ep:" << Ep << " and FinalFactor is : "<< finalFactor <<  " , and R is " << R<<  " , and a and b are: "<< a << " , "<< b << endl;
  return finalFactor;
}
*/
double weight_ytar_corr(double ytar , TString spectrometer_option){
double weight;
if ( spectrometer_option=="HMS") { weight =  -0.00812174*ytar*ytar - 0.0000415678*ytar+1.00021;}
else if ( spectrometer_option=="SHMS") { weight = 1.0;}
return weight;
}
double weight_MC_Jacobian_corr(double xptar , double yptar, TString spectrometer_option){

double weight;
if ( spectrometer_option=="HMS") { weight =  TMath::Sqrt(1+xptar*xptar + yptar*yptar);  }
else if ( spectrometer_option=="SHMS") { weight = 1.0;}
 double final_weight = (1.0/weight) * (1.0/weight) * (1.0/weight)  ;
return final_weight;

}


double weight_Isoscalar_corr(double Xb , double Q2 , int N , int Z) {

    double xtab[] = {
      0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090,
      0.100, 0.110, 0.120, 0.130, 0.140, 0.150, 0.160, 0.170, 0.180,
      0.190, 0.200, 0.210, 0.220, 0.230, 0.240, 0.250, 0.260, 0.270,
      0.280, 0.290, 0.300, 0.310, 0.320, 0.330, 0.340, 0.350, 0.360,
      0.370, 0.380, 0.390, 0.400, 0.410, 0.420, 0.430, 0.440, 0.450,
      0.460, 0.470, 0.480, 0.490, 0.500, 0.510, 0.520, 0.530, 0.540,
      0.550, 0.560, 0.570, 0.580, 0.590, 0.600, 0.610, 0.620, 0.630,
      0.640, 0.650, 0.660, 0.670, 0.680, 0.690, 0.700, 0.710, 0.720,
      0.730, 0.740, 0.750, 0.760, 0.770, 0.780, 0.790, 0.800, 0.810,
      0.820, 0.830, 0.840, 0.850, 0.860, 0.870, 0.880, 0.890, 0.900,
      0.910, 0.920, 0.930, 0.940, 0.950, 0.960, 0.970, 0.980, 0.990, 1.000
    };

    double q2tab[] = {
      1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 15.0
    };

    double nonp[100][10] = { 
      {0.983509,0.967562,0.952133,0.937198,0.922738,0.908734,0.895167,0.882020,0.869275,0.856917},
      {0.844931,0.833300,0.822010,0.811048,0.800401,0.790055,0.779999,0.770222,0.760711,0.751457},
      {0.742449,0.733679,0.725136,0.716812,0.708699,0.700790,0.693076,0.685551,0.678208,0.671042},
      {0.664041,0.657206,0.650530,0.644006,0.637628,0.631395,0.625297,0.619336,0.613508,0.607788},
      {0.602222,0.596735,0.591388,0.586151,0.581032,0.575976,0.571070,0.566249,0.561505,0.556859},
      {0.552326,0.547884,0.543518,0.539216,0.535019,0.530865,0.526926,0.522888,0.518921,0.515048},
      {0.511480,0.507615,0.504017,0.500443,0.497119,0.493618,0.490066,0.486863,0.483462,0.480752},
      {0.477100,0.473876,0.471597,0.468325,0.465172,0.463005,0.459858,0.457916,0.452807,0.450750},
      {0.451604,0.446949,0.438826,0.443575,0.444461,0.436106,0.428253,0.432000,0.436110,0.426526},
      {0.417893,0.421976,0.433359,0.412051,0.425995,0.426235,0.414856,0.419595,0.414556,0.457513},
        
      {0.983509,0.967559,0.952119,0.937164,0.922671,0.908618,0.894987,0.881756,0.868907,0.856424},
      {0.844289,0.832487,0.821001,0.809819,0.798926,0.788308,0.777955,0.767855,0.757995,0.748366},
      {0.738959,0.729762,0.720769,0.711970,0.703358,0.694925,0.686665,0.678569,0.670634,0.662851},
      {0.655216,0.647724,0.640369,0.633148,0.626055,0.619087,0.612238,0.605507,0.598888,0.592381},
      {0.585978,0.579676,0.573483,0.567381,0.561373,0.555459,0.549643,0.543901,0.538247,0.532687},
      {0.527185,0.521772,0.516424,0.511157,0.505953,0.500826,0.495763,0.490746,0.485828,0.480940},
      {0.476151,0.471344,0.466696,0.462116,0.457552,0.452952,0.448466,0.444209,0.439813,0.435568},
      {0.431184,0.427269,0.423081,0.419041,0.415156,0.410801,0.407098,0.403425,0.400097,0.395570},
      {0.391607,0.389038,0.386063,0.381268,0.377802,0.374977,0.371842,0.368327,0.364838,0.362433},
      {0.360181,0.353288,0.355565,0.350272,0.353403,0.343797,0.356637,0.354206,0.348619,0.313281},
        
      {0.983510,0.967561,0.952118,0.937157,0.922654,0.908587,0.894934,0.881676,0.868794,0.856270},
      {0.844086,0.832226,0.820675,0.809418,0.798441,0.787732,0.777276,0.767064,0.757083,0.747323},
      {0.737774,0.728427,0.719272,0.710302,0.701508,0.692884,0.684421,0.676114,0.667956,0.659941},
      {0.652064,0.644319,0.636702,0.629208,0.621833,0.614572,0.607422,0.600379,0.593440,0.586600},
      {0.579859,0.573212,0.566655,0.560190,0.553808,0.547511,0.541293,0.535154,0.529089,0.523092},
      {0.517166,0.511309,0.505516,0.499776,0.494105,0.488500,0.482945,0.477457,0.472016,0.466636},
      {0.461332,0.456062,0.450841,0.445676,0.440586,0.435560,0.430543,0.425598,0.420714,0.415880},
      {0.411154,0.406366,0.401623,0.397016,0.392580,0.388002,0.383386,0.379056,0.374824,0.370799},
      {0.366380,0.361443,0.357684,0.354736,0.350574,0.345649,0.342257,0.339736,0.333568,0.331548},
      {0.328415,0.326729,0.320067,0.318466,0.316831,0.318958,0.307941,0.306496,0.314458,0.341841},
        
      {0.983511,0.967562,0.952118,0.937155,0.922648,0.908574,0.894911,0.881641,0.868742,0.856198},
      {0.843990,0.832102,0.820519,0.809225,0.798207,0.787451,0.776944,0.766676,0.756634,0.746807},
      {0.737187,0.727763,0.718526,0.709468,0.700581,0.691858,0.683291,0.674874,0.666600,0.658464},
      {0.650460,0.642582,0.634826,0.627188,0.619662,0.612245,0.604933,0.597723,0.590610,0.583591},
      {0.576664,0.569826,0.563074,0.556404,0.549815,0.543304,0.536865,0.530498,0.524198,0.517961},
      {0.511786,0.505672,0.499617,0.493619,0.487675,0.481789,0.475957,0.470175,0.464451,0.458773},
      {0.453154,0.447575,0.442062,0.436587,0.431165,0.425793,0.420489,0.415204,0.409975,0.404812},
      {0.399690,0.394627,0.389593,0.384685,0.379695,0.374906,0.370107,0.365316,0.360556,0.355935},
      {0.351579,0.346887,0.342370,0.337866,0.334286,0.329150,0.324774,0.321725,0.318475,0.312829},
      {0.310749,0.309096,0.300508,0.302876,0.296301,0.297655,0.289021,0.289589,0.293994,0.314936},
        
      {0.983511,0.967564,0.952119,0.937155,0.922645,0.908567,0.894900,0.881622,0.868714,0.856158},
      {0.843936,0.832031,0.820429,0.809113,0.798070,0.787287,0.776750,0.766448,0.756369,0.746503},
      {0.736840,0.727370,0.718083,0.708972,0.700029,0.691246,0.682616,0.674131,0.665787,0.657576},
      {0.649494,0.641535,0.633694,0.625966,0.618348,0.610834,0.603421,0.596106,0.588885,0.581754},
      {0.574711,0.567753,0.560876,0.554079,0.547357,0.540709,0.534130,0.527617,0.521165,0.514775},
      {0.508440,0.502162,0.495938,0.489768,0.483649,0.477581,0.471566,0.465599,0.459680,0.453810},
      {0.447989,0.442215,0.436491,0.430812,0.425180,0.419597,0.414061,0.408569,0.403131,0.397743},
      {0.392384,0.387107,0.381856,0.376672,0.371528,0.366410,0.361420,0.356439,0.351487,0.346606},
      {0.341924,0.337194,0.332407,0.327782,0.323699,0.319170,0.314455,0.310125,0.307215,0.301642},
      {0.298502,0.295424,0.289368,0.289307,0.284499,0.280933,0.277798,0.277803,0.280455,0.284791},
        
      {0.983512,0.967565,0.952120,0.937156,0.922644,0.908564,0.894893,0.881610,0.868696,0.856133},
      {0.843902,0.831987,0.820372,0.809042,0.797983,0.787181,0.776624,0.766300,0.756197,0.746305},
      {0.736613,0.727111,0.717792,0.708646,0.699665,0.690842,0.682169,0.673640,0.665248,0.656987},
      {0.648852,0.640838,0.632939,0.625151,0.617470,0.609890,0.602409,0.595023,0.587727,0.580520},
      {0.573398,0.566357,0.559396,0.552510,0.545697,0.538955,0.532279,0.525664,0.519108,0.512609},
      {0.506164,0.499771,0.493430,0.487140,0.480898,0.474704,0.468557,0.462460,0.456407,0.450400},
      {0.444440,0.438524,0.432654,0.426828,0.421050,0.415315,0.409625,0.403984,0.398382,0.392835},
      {0.387346,0.381878,0.376478,0.371131,0.365850,0.360577,0.355399,0.350270,0.345240,0.340166},
      {0.335218,0.330465,0.325688,0.320838,0.316225,0.311863,0.307204,0.302922,0.299101,0.294904},
      {0.290491,0.287610,0.284479,0.279501,0.279012,0.274414,0.273407,0.272938,0.272034,0.276103},
        
      {0.983512,0.967566,0.952121,0.937156,0.922644,0.908562,0.894889,0.881603,0.868685,0.856116},
      {0.843879,0.831956,0.820332,0.808992,0.797922,0.787107,0.776536,0.766196,0.756076,0.746165},
      {0.736453,0.726930,0.717587,0.708416,0.699408,0.690556,0.681852,0.673291,0.664865,0.656569},
      {0.648396,0.640342,0.632402,0.624571,0.616843,0.609217,0.601686,0.594248,0.586899,0.579636},
      {0.572456,0.565356,0.558332,0.551383,0.544504,0.537692,0.530944,0.524256,0.517624,0.511045},
      {0.504518,0.498042,0.491614,0.485235,0.478903,0.472616,0.466375,0.460178,0.454027,0.447919},
      {0.441853,0.435834,0.429856,0.423923,0.418031,0.412188,0.406386,0.400631,0.394916,0.389253},
      {0.383636,0.378070,0.372545,0.367084,0.361686,0.356316,0.351006,0.345767,0.340628,0.335499},
      {0.330442,0.325487,0.320655,0.315780,0.311117,0.306730,0.302236,0.297676,0.293601,0.289589},
      {0.285513,0.282426,0.279029,0.274611,0.273444,0.268930,0.269485,0.267128,0.266523,0.270773},
        
      {0.983512,0.967567,0.952122,0.937157,0.922644,0.908561,0.894886,0.881598,0.868677,0.856104},
      {0.843862,0.831934,0.820304,0.808956,0.797877,0.787053,0.776472,0.766120,0.755987,0.746062},
      {0.736335,0.726796,0.717435,0.708245,0.699217,0.690343,0.681617,0.673032,0.664580,0.656257},
      {0.648056,0.639972,0.632001,0.624137,0.616375,0.608712,0.601144,0.593668,0.586278,0.578973},
      {0.571749,0.564604,0.557533,0.550535,0.543605,0.536741,0.529939,0.523194,0.516503,0.509864},
      {0.503275,0.496734,0.490241,0.483794,0.477391,0.471033,0.464719,0.458448,0.452220,0.446034},
      {0.439890,0.433788,0.427729,0.421711,0.415736,0.409803,0.403917,0.398072,0.392272,0.386518},
      {0.380816,0.375159,0.369550,0.363991,0.358501,0.353061,0.347664,0.342349,0.337116,0.331931},
      {0.326810,0.321788,0.316898,0.312051,0.307267,0.302752,0.298231,0.293810,0.289751,0.285784},
      {0.281724,0.278268,0.275267,0.271325,0.269301,0.266271,0.266834,0.263387,0.264002,0.267541},
        
      {0.983513,0.967568,0.952124,0.937158,0.922645,0.908561,0.894883,0.881592,0.868666,0.856088},
      {0.843840,0.831904,0.820265,0.808907,0.797817,0.786980,0.776384,0.766016,0.755865,0.745921},
      {0.736173,0.726611,0.717226,0.708009,0.698953,0.690050,0.681292,0.672672,0.664185,0.655824},
      {0.647584,0.639458,0.631443,0.623532,0.615723,0.608009,0.600389,0.592857,0.585411,0.578046},
      {0.570760,0.563550,0.556413,0.549346,0.542345,0.535406,0.528526,0.521701,0.514926,0.508201},
      {0.501523,0.494891,0.488303,0.481759,0.475257,0.468797,0.462378,0.456000,0.449662,0.443364},
      {0.437106,0.430889,0.424711,0.418574,0.412477,0.406423,0.400410,0.394441,0.388517,0.382636},
      {0.376805,0.371023,0.365294,0.359616,0.353994,0.348441,0.342946,0.337513,0.332155,0.326902},
      {0.321704,0.316604,0.311612,0.306770,0.301954,0.297360,0.292914,0.288574,0.284420,0.280509},
      {0.277032,0.273768,0.270723,0.267774,0.265423,0.263488,0.263733,0.262517,0.260921,0.264691},
        
      {0.983514,0.967570,0.952126,0.937160,0.922646,0.908560,0.894880,0.881584,0.868654,0.856069},
      {0.843813,0.831867,0.820216,0.808846,0.797740,0.786886,0.776271,0.765882,0.755708,0.745738},
      {0.735962,0.726370,0.716953,0.707702,0.698609,0.689666,0.680866,0.672202,0.663667,0.655256},
      {0.646963,0.638781,0.630707,0.622735,0.614861,0.607080,0.599389,0.591784,0.584260,0.576816},
      {0.569447,0.562150,0.554923,0.547762,0.540664,0.533624,0.526639,0.519704,0.512816,0.505974},
      {0.499175,0.492418,0.485701,0.479024,0.472386,0.465786,0.459224,0.452699,0.446211,0.439760},
      {0.433346,0.426969,0.420630,0.414329,0.408066,0.401844,0.395662,0.389523,0.383428,0.377380},
      {0.371379,0.365430,0.359536,0.353700,0.347925,0.342220,0.336589,0.331033,0.325563,0.320198},
      {0.314934,0.309777,0.304761,0.299900,0.295185,0.290643,0.286344,0.282260,0.278380,0.274907},
      {0.271685,0.268766,0.266372,0.264444,0.262996,0.262203,0.261550,0.261707,0.262413,0.261832}



};

    double f2rat = 1.0;
    double ISO;

    // Loop over xtab to find the correct range for Xb
    for (int i = 0; i < 99; i++) {
      if (Xb >= xtab[i] && Xb < xtab[i + 1]) {
	double xlo = xtab[i];
	double xhi = xtab[i + 1];
	double deltax = xhi - Xb;

	// Loop over q2tab to find the correct range for Q2
	for (int j = 0; j < 9; j++) {
	  if (Q2 >= q2tab[j] && Q2 < q2tab[j + 1]) {
	    double q2lo = q2tab[j];
	    double q2hi = q2tab[j + 1];
	    double deltaq2 = q2hi - q2lo;

	    // Linear interpolation between values
	    double A1 = nonp[i][j];
	    double A2 = nonp[i + 1][j];
	    double A3 = nonp[i][j + 1];
	    double A4 = nonp[i + 1][j + 1];

	    double A12 = (A1 * (xhi - Xb) + A2 * (Xb - xlo)) / deltax;
	    double A34 = (A3 * (xhi - Xb) + A4 * (Xb - xlo)) / deltax;

	    double A1234 = (A12 * (q2hi - Q2) + A34 * (Q2 - q2lo)) / deltaq2;

	    f2rat = A1234;

	    // Apply John's scale factor
	    f2rat = (1.0 - 0.013 * log(Q2 / 16.0)) * (1.0 + f2rat) - 1.0;

	    //return f2rat;
	  }
	}
      }
    }
    //return f2rat; // Default return if no match found
    double A =N+Z;
    ISO= (A/2.0)*( ( 1+f2rat )  / ( Z + N*f2rat )    );
    return ISO;
  }



double weight_delta_corr_Dave(double delta ,  TString spectrometer_option) {
  // this is applied as a weight to the MC, rather than a redefinition of the delta variable
  double delta_corr=1.0;
  double deltatmp;
  if (delta >= 8.0) {
    deltatmp = 8.0;
  } else if (delta <= -8.0) {
    deltatmp = -8.0;
  } else {
    deltatmp = delta;
  }


  double h1 = 1.0069;
  double h2 = 0.34018e-02;
  double h3 = -0.71161e-03;
  double h4 = -0.12060e-04;
  double h5 = 0.11322e-04;
  double h6 = -0.78222e-06;

  if ( spectrometer_option=="HMS") { 

    delta_corr =  h1 + h2 * deltatmp + h3 * TMath::Power(deltatmp, 2) 
                + h4 * TMath::Power(deltatmp, 3) 
                + h5 * TMath::Power(deltatmp, 4) 
                + h6 * TMath::Power(deltatmp, 5);
}
  else if ( spectrometer_option=="SHMS") { delta_corr = 1.0;}

  return delta_corr;

}


double weight_delta_corr(double delta ,  TString spectrometer_option) {

double delta_corr;
if ( spectrometer_option=="HMS") { delta_corr =  0.990337*delta - 0.00236077*TMath::Power(delta,2) + 0.000286814*TMath::Power(delta,3) + 2.09878E-6*TMath::Power(delta,4) - 2.4867E-6*TMath::Power(delta,5) + 1.8646E-7*TMath::Power(delta,6) ;   }
else if ( spectrometer_option=="SHMS") { delta_corr = delta;}

return delta_corr;

}



// Function to calculate hmscereff_xem2
double hmscereff_xem2(double delta ,  TString spectrometer_option) {
  // Parameters
  const double a1 = 0.9994;
  const double b1 = 0.99727;
  const double b2 = -0.00213;
  const double c1 = 0.9962;

  double eff;
  if ( spectrometer_option=="HMS") {
    if (delta <= -1.0) {     eff = a1;   }
    else if (delta > -1.0 && delta < 0.5) {     eff = b1 + b2 * delta;  }
    else {     eff = c1;  }
  }
  else if ( spectrometer_option=="SHMS") {  eff =1.0;}
  return eff;
}


std::pair<double, double> calculateTotalEntriesAndError(TH1* hist) {
  double totalEntries = 0;
  double totalErrorSquared = 0;

  // Loop over all bins
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    double binContent = hist->GetBinContent(i);
    totalEntries += binContent;
    totalErrorSquared += binContent; // For Poisson statistics, the variance is equal to the bin content
  }

  double totalError = sqrt(totalErrorSquared);

  std::cout << "Total Entries: " << totalEntries << std::endl;
  std::cout << "Total Error: " << totalError << std::endl;

  return std::make_pair(totalEntries, totalError);
}

void plotHistograms_Dummy(TH1F* X_Data_LT, TH1F* X_Data_LT_Dummy,  TString out_name);
void plotHistograms_Contamination(TH1F* X_Data_ST, TH1F* X_Data_ST_tmp,  TString out_name);
std::vector<ScaleFactors> run_data_LT , run_data_ST , run_data_LT_Dummy , run_data_ST_tmp; // run_data_ST_tmp is for contamination subtraction (Ca48 , B10 B11 targets)   

void cross_Sections_updated_workingProcess (const char* target_in = "Ca12", const char* angle_in = "20", const char* momentum_spec_in = "2p42", TString spectrometer_option="HMS" , TString ytarg_corr_opt ="ON", TString jacob_corr_opt ="ON" , TString delta_corr_opt="ON" , TString CSB_corr_opt="ON" , TString shift_collimator_opt="00") {

  TString target_input        = target_in;
  TString angle_input         = angle_in;
  TString momentum_spec_input = momentum_spec_in;
  TString spectrometer        = spectrometer_option;

  bool ytar_corr_ON , MC_Jacobian_corr_ON , delta_corr_ON , Coulomb_corr_ON , CSB_corr_ON; 
  if (ytarg_corr_opt =="ON"){ytar_corr_ON        = true ; }       else {ytar_corr_ON        = false;}
  if (jacob_corr_opt =="ON"){MC_Jacobian_corr_ON = true ; }       else {MC_Jacobian_corr_ON = false;}
  if (delta_corr_opt =="ON"){delta_corr_ON       = true ; }       else {delta_corr_ON       = false;}
  if (CSB_corr_opt   =="ON"){CSB_corr_ON         = true ; }       else {CSB_corr_ON         = false;}

  TString spectrometer_lowercase;
  if (      spectrometer=="HMS" || spectrometer=="test") {spectrometer_lowercase = "hms";}
  else if ( spectrometer=="SHMS")                        {spectrometer_lowercase = "shms";}

  TString f_st_mc_st , f_lt_mc_st  ,mc_file_lt_st , mc_file_st_st  ; 
  double angle , p_spec , AM_ST , shift_collimator;
  int delta_hi , delta_lo , N, Z;
  double density_ST, length_ST  ;
  if (shift_collimator_opt=="00"){       shift_collimator=0.0 ; }
  else if (shift_collimator_opt=="p1"){  shift_collimator=-1.0; }  // p1 for +1cm shift, it carries a - just bc of the definition later
  else if (shift_collimator_opt=="p2"){  shift_collimator=-2.0; }
  else if (shift_collimator_opt=="p3"){  shift_collimator=-3.0; }
  else if (shift_collimator_opt=="m1"){  shift_collimator= 1.0; }
  else if (shift_collimator_opt=="m2"){  shift_collimator= 2.0; }

  double slope = 2.545 ;
  double vertical_shift = shift_collimator*sqrt(1 + pow(slope, 2));


  cout << "SPECTROMETER: "<<spectrometer<<endl;

    angle=angle_input.Atof();
    angle = angle +  0.0974 ; // 0.0017 radiasn is 0.0974 degrees

    //  angle=20.005;
  if  ( spectrometer=="test" ){
    //WORK ON THIS A BIT MORE
    mc_file_lt_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_d2cryo22_hms.out";
    mc_file_st_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_12carbon22_hms.out";
    density_ST    = 2.00;    length_ST     = 0.287;    delta_hi      = 20;    delta_lo      = -20;
    f_st_mc_st    = "/work/smoran/xem2/sim/test_5p87_carbon.root";    f_lt_mc_st    = "/work/smoran/xem2/sim/test_5p87_deuterium.root";
    p_spec  = 5.004878;    cout<< "RUNNING TEST CASE: " <<endl;
    run_data_ST       = {      {5828 , 142189  , 1 , 0.9969 , 1       , 0.999881 , 2}    };
    run_data_LT       = {      {4780 , 90719.6 , 1 , 0.9972 , 1.00004 , 0.999827 , 2}    };
    run_data_LT_Dummy = {      {4772 , 79655.2 , 1 , 0.9976 , 1.00004 , 0.999886 , 2}    };

}


  mc_file_st_st = Form("/work/smoran/xem2/cross_sections/EXTERNALS_tables/OUTV4_IT3/xem2_emc_rc_%s_22_hms.out" , (const char*)target_input);



  if      (target_input=="He03") {density_ST    = 0.041; length_ST     = 9.9980;  AM_ST         = 3.0;   }
  else if (target_input=="LH02") {density_ST    = 0.07248;length_ST    = 9.9980;  AM_ST         = 1.0;   }
  else if (target_input=="He04") {density_ST    = 0.1293;length_ST     = 9.9980;  AM_ST         = 4.0;   }
  else if (target_input=="Li06") {density_ST    = 0.463; length_ST     = 0.53347730; AM_ST      = 6.0;   } // modified to get the numbers Dave got for the thickness  
  else if (target_input=="Li07") {density_ST    = 0.540; length_ST     = 0.49272835; AM_ST      = 7.0;   } // modified to get the numbers Dave got for the thickness 
  else if (target_input=="Be09") {density_ST    = 1.848; length_ST     = 0.535;   AM_ST         = 9.0;   }
  else if (target_input=="BC10") {density_ST    = 2.325 ;length_ST     = 0.2450;  AM_ST         = 10.0;  }
  else if (target_input=="BC11") {density_ST    = 2.434 ;length_ST     = 0.2600;  AM_ST         = 11.0;  }
  else if (target_input=="Ca12") {density_ST    = 2.00;  length_ST     = 0.287;   AM_ST         = 12.0;  }
  else if (target_input=="Al27") {density_ST    = 2.699; length_ST     = 0.1704;  AM_ST         = 27.0;  }
  else if (target_input=="Ca40") {density_ST    = 1.546; length_ST     = 0.5079;  AM_ST         = 40.0;  }
  else if (target_input=="Ca48") {density_ST    = 1.855; length_ST     = 0.5667;  AM_ST         = 48.0;  }
  else if (target_input=="Ti48") {density_ST    = 4.540; length_ST     = 0.0648;  AM_ST         = 48.0;  }
  else if (target_input=="Fe54") {density_ST    = 7.605; length_ST     = 0.0483;  AM_ST         = 54.0;  }
  else if (target_input=="Ni58") {density_ST    = 8.787; length_ST     = 0.0274;  AM_ST         = 58.0;  }
  else if (target_input=="Cu63") {density_ST    = 8.960; length_ST     = 0.1051;  AM_ST         = 63.0;  }
  else if (target_input=="Ni64") {density_ST    = 9.696; length_ST     = 0.0269;  AM_ST         = 64.0;  }
  else if (target_input=="Ag10") {density_ST    = 10.50; length_ST     = 0.0503;  AM_ST         = 108.0; }
  else if (target_input=="Sn11") {density_ST    = 7.310; length_ST     = 0.0624;  AM_ST         = 119.0; }
  else if (target_input=="Au19") {density_ST    = 19.32; length_ST     = 0.0209;  AM_ST         = 197.0; }
  else if (target_input=="Th23") {density_ST    = 11.72; length_ST     = 0.0349;  AM_ST         = 232.0; }

  delta_hi      =  10 ; //20;  10 is for daves simulations
  delta_lo      = -10 ; //-20; -10 is for daves simulations
  mc_file_lt_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/OUTV4_IT3/xem2_emc_rc_LD2_22_hms.out";
  f_lt_mc_st    = Form("/work/smoran/xem2/sim/V2/%s_%sdeg_deuterium_%s.root", (const char*)spectrometer_lowercase   ,  (const char*)angle_input , (const char*)momentum_spec_input) ;
  f_st_mc_st    = Form("/work/smoran/xem2/sim/V2/%s_%sdeg_%s_%s.root"   , (const char*)spectrometer_lowercase ,  (const char*)angle_input  ,  (const char*)target_input , (const char*)momentum_spec_input) ;

  TString modified_momentum = momentum_spec_input;
  modified_momentum.ReplaceAll("p", ".");
  p_spec = modified_momentum.Atof();

  //Dave's offset applied to the central momentum
  p_spec = p_spec * (1.0 - 0.0049738 - 0.0044562 * exp(-0.5 * pow((p_spec - 4.6682) / 0.51466, 2)));


  //  cout<<__LINE__<<endl;
  fillScaleFactors("run_numbers_info.txt", "LD2"        , angle_input, spectrometer, momentum_spec_input, p_spec,run_data_LT       );
  fillScaleFactors("run_numbers_info.txt", "Dummy"      , angle_input, spectrometer, momentum_spec_input, p_spec,run_data_LT_Dummy );
  fillScaleFactors("run_numbers_info.txt", target_input , angle_input, spectrometer, momentum_spec_input, p_spec,run_data_ST       );

  //special targets, Ca48, BC10, BC11, they need contamination subtraction

  if (target_input=="Ca48")                         {fillScaleFactors("run_numbers_info.txt", "Ca40" , angle_input, spectrometer, momentum_spec_input, p_spec,run_data_ST_tmp       ); }
  if (target_input=="BC10" || target_input=="BC11") {fillScaleFactors("run_numbers_info.txt", "Ca12" , angle_input, spectrometer, momentum_spec_input, p_spec,run_data_ST_tmp       ); }

  cout << "The target input is: " << target_input << " , and the p is = "<<momentum_spec_input <<endl;
  TString out_name = "cross_sections_" + target_input + "_angle_" + angle_input + "_momentum_" + momentum_spec_input + "_" + spectrometer;
  const double theta_central = angle*3.14159265/180.0  + 0.0017;  // the 0.0017 is dave's offset to the kinematics'
  Variables st_vars_mc, lt_vars_mc ;

  TFile *f_st_mc   = new TFile((const char*) f_st_mc_st  );
  TFile *f_lt_mc   = new TFile((const char*) f_lt_mc_st  );
  TTree *t_st_mc ,*t_lt_mc;
  if (spectrometer=="HMS" || spectrometer=="test"){
	  t_st_mc   = (TTree*) f_st_mc  ->Get("h10");
	  t_lt_mc   = (TTree*) f_lt_mc  ->Get("h10");
	}
  else if (spectrometer=="SHMS"){
  cout<< "THE SPECTROMETER IS SHMS!"<<endl;
  
	  t_st_mc   = (TTree*) f_st_mc  ->Get("h1411");
	  t_lt_mc   = (TTree*) f_lt_mc  ->Get("h1411");
	}	
	
  setBranchAddresses_mc(t_st_mc    , st_vars_mc,spectrometer);
  setBranchAddresses_mc(t_lt_mc    , lt_vars_mc,spectrometer);

double low_lim_hist , high_lim_histo;

if (spectrometer=="HMS" || spectrometer=="test")      {
  if (var_dep == "delta") {  low_lim_hist=-8.0  ; high_lim_histo=8.0;} // if you are binning in delta
  if (var_dep == "xb")    {  low_lim_hist=0.1   ; high_lim_histo=1.5;} // if you are binning in xb
}
else if (spectrometer=="SHMS"){low_lim_hist=-10.1 ; high_lim_histo=22.1;}

  TH1F *X_MC_LT = new TH1F("X SIM LT ", "X SIM LT", N_BINS, low_lim_hist, high_lim_histo);

  if (var_dep == "delta") { X_MC_LT->GetXaxis()->SetTitle("delta");}
  else if (var_dep == "xb") { X_MC_LT->GetXaxis()->SetTitle("xb");}

  X_MC_LT->GetYaxis()->SetTitle("Counts");

  TH1F *X_MC_ST         = (TH1F*)X_MC_LT->Clone();
  TH1F *X_Data_ST       = (TH1F*)X_MC_LT->Clone();
  TH1F *X_Data_LT       = (TH1F*)X_MC_LT->Clone();
  TH1F *X_Data_LT_Dummy = (TH1F*)X_MC_LT->Clone();
  TH1F *X_Data_ST_tmp   = (TH1F*)X_MC_LT->Clone();

  // just for testing purposes:
  TH1F *X_MC_ST_Unweighted   = (TH1F*)X_MC_LT->Clone();
  TH1F *X_MC_LT_Unweighted   = (TH1F*)X_MC_LT->Clone();

  X_MC_ST  ->Sumw2();
  X_MC_LT  ->Sumw2();
  X_Data_ST->Sumw2();
  X_Data_LT->Sumw2();
  X_Data_ST_tmp->Sumw2();
  X_Data_LT_Dummy->Sumw2();
  X_MC_ST_Unweighted->Sumw2();
  X_MC_LT_Unweighted->Sumw2();

  vector<double>   V1_lt,V2_lt,V3_lt,V4_lt,V5_lt,V6_lt,V7_lt,V8_lt,V9_lt,V10_lt,V11_lt,V12_lt,V13_lt;
  vector<double>   V1_st,V2_st,V3_st,V4_st,V5_st,V6_st,V7_st,V8_st,V9_st,V10_st,V11_st,V12_st,V13_st;

  ImportRadcor(  V1_lt,V2_lt,V3_lt,V4_lt,V5_lt,V6_lt,V7_lt,V8_lt,V9_lt,V10_lt,V11_lt,V12_lt,V13_lt   ,  (const char*)mc_file_lt_st  );  
  ImportRadcor(  V1_st,V2_st,V3_st,V4_st,V5_st,V6_st,V7_st,V8_st,V9_st,V10_st,V11_st,V12_st,V13_st   ,  (const char*)mc_file_st_st  );

  const int size_st= V1_st.size();
  const int size_lt= V1_lt.size();
  cout<<"size of V1_ST: "<<size_st<< " , and V1_LT: "<< size_lt <<endl;

  TGraph2DErrors* gr2D_SigmaRad_ST  		= createGraph2D(V2_st, V3_st, V9_st, size_st);
  TGraph2DErrors* gr2D_SigmaBorn_ST 		= createGraph2D(V2_st, V3_st, V6_st, size_st);
  TGraph2DErrors* gr2D_SigmaRad_LT  		= createGraph2D(V2_lt, V3_lt, V9_lt, size_lt);
  TGraph2DErrors* gr2D_SigmaBorn_LT 		= createGraph2D(V2_lt, V3_lt, V6_lt, size_lt);
  TGraph2DErrors* gr2D_CoulombCorrection_ST     = createGraph2D(V2_st, V3_st, V13_st,size_st);

  Long64_t nentries_st_mc   = t_st_mc->GetEntries(); 
  Long64_t nentries_lt_mc   = t_lt_mc->GetEntries();
  double mc_scale_ST , mc_scale_LT;

  mc_scale_ST = calculate_mc_scale_factor( nentries_st_mc, p_spec, density_ST, length_ST, AM_ST, sim_charge_ST, delta_hi, delta_lo , "ST" );
  mc_scale_LT = calculate_mc_scale_factor( nentries_lt_mc, p_spec, density_LT, length_LT, AM_LT, sim_charge_LT, delta_hi, delta_lo , "LT" );

  cout<< "MC scale afctors: " << endl;

  cout << "ST: " << mc_scale_ST  <<endl;
  cout << "LT: " << mc_scale_LT <<endl;




  if (trigger_select=="ELREAL"){
  // here Im keeping all the entries with trigger == 2 (ELREAL)
  run_data_ST.erase(std::remove_if(run_data_ST.begin(), run_data_ST.end(), [](const ScaleFactors& data) {
        return data.trigger != 2;
      }), run_data_ST.end());
    
  run_data_LT.erase(std::remove_if(run_data_LT.begin(), run_data_LT.end(), [](const ScaleFactors& data) {
        return data.trigger != 2;
      }), run_data_LT.end());
    
  run_data_LT_Dummy.erase(std::remove_if(run_data_LT_Dummy.begin(), run_data_LT_Dummy.end(), [](const ScaleFactors& data) {
        return data.trigger != 2;
      }), run_data_LT_Dummy.end());
  if (target_input=="Ca48" || target_input=="BC10" || target_input=="BC11") {   
    run_data_ST_tmp.erase(std::remove_if(run_data_ST_tmp.begin(), run_data_ST_tmp.end(), [](const ScaleFactors& data) {
	  return data.trigger != 2;
	}), run_data_ST_tmp.end());

     }

  }    

  ////////////////////
  /////   DATA   /////
  ////////////////////

  double total_SF_ST_tmp = 0.0; double production_scale_ST_tmp ; 

    double total_SF_ST = 0.0;
    for (const auto& data : run_data_ST) { total_SF_ST += (data.charge*data.Eff_Fid*data.Eff_Elec_LiveTime*data.Eff_Comp_LiveTime) / data.PS; }
    double production_scale_ST = 1.0 / total_SF_ST;     cout<< "production_scale_ST : " << production_scale_ST <<endl;
    double total_SF_LT = 0.0;
    for (const auto& data : run_data_LT) { total_SF_LT += (data.charge*data.Eff_Fid*data.Eff_Elec_LiveTime*data.Eff_Comp_LiveTime) / data.PS; }
    double production_scale_LT = 1.0 / total_SF_LT;     cout<< "production_scale_LT : " << production_scale_LT <<endl;
    double total_SF_LT_Dummy = 0.0;
    for (const auto& data : run_data_LT_Dummy) { total_SF_LT_Dummy += (data.charge*data.Eff_Fid*data.Eff_Elec_LiveTime*data.Eff_Comp_LiveTime) / data.PS; }
    double production_scale_LT_Dummy = 1.0 / total_SF_LT_Dummy;     cout<< "production_scale_LT_Dummy : " << production_scale_LT_Dummy <<endl;

    if (target_input=="Ca48" || target_input=="BC10" || target_input=="BC11") {    
      for (const auto& data : run_data_ST_tmp) { total_SF_ST_tmp += (data.charge*data.Eff_Fid*data.Eff_Elec_LiveTime*data.Eff_Comp_LiveTime) / data.PS; }
      production_scale_ST_tmp = 1.0 / total_SF_ST_tmp;     cout<< "production_scale_ST (contamination subtraction) : " << production_scale_ST_tmp <<endl;
    }

    double cuts_delta_min_st , cuts_delta_max_st , cuts_ytar_min_st , cuts_ytar_max_st , cuts_ph_min_st , cuts_ph_max_st , cuts_th_min_st , cuts_th_max_st , cuts_npeSum_st , cuts_etotTrack_st;
    double cuts_delta_min_lt , cuts_delta_max_lt , cuts_ytar_min_lt , cuts_ytar_max_lt , cuts_ph_min_lt , cuts_ph_max_lt , cuts_th_min_lt , cuts_th_max_lt , cuts_npeSum_lt , cuts_etotTrack_lt;
    double cuts_delta_min_st_mc, cuts_delta_max_st_mc, cuts_ytar_min_st_mc, cuts_ytar_max_st_mc, cuts_ph_min_st_mc, cuts_ph_max_st_mc, cuts_th_min_st_mc, cuts_th_max_st_mc, cuts_npeSum_st_mc, cuts_etotTrack_st_mc;
    double cuts_delta_min_lt_mc, cuts_delta_max_lt_mc, cuts_ytar_min_lt_mc, cuts_ytar_max_lt_mc, cuts_ph_min_lt_mc, cuts_ph_max_lt_mc, cuts_th_min_lt_mc, cuts_th_max_lt_mc, cuts_npeSum_lt_mc, cuts_etotTrack_lt_mc;

    std::ofstream outFile_runs_ST(       out_name + "_yield_runs_ST.txt",     std::ios_base::app );
    std::ofstream outFile_runs_LT(       out_name + "_yield_runs_LT.txt",     std::ios_base::app );
    std::ofstream outFile_runs_LT_Dummy( out_name + "_yield_runs_Al_LT.txt",  std::ios_base::app );
    std::ofstream outFile_runs_ST_tmp(   out_name + "_yield_runs_ST_tmp.txt", std::ios_base::app );

    outFile_runs_ST       << "run\ttrigger\tYield\tError\n";
    outFile_runs_LT       << "run\ttrigger\tYield\tError\n";
    outFile_runs_LT_Dummy << "run\ttrigger\tYield\tError\n";

    if (target_input=="Ca48" || target_input=="BC10" || target_input=="BC11") {    outFile_runs_ST_tmp       << "run\ttrigger\tYield\tError\n";     }


    if (spectrometer=="HMS" || spectrometer=="test"){
      cuts_delta_min_st  = -8.0 ;
      cuts_delta_max_st  = 8.0 ;
      cuts_ytar_min_st  = -100;
      cuts_ytar_max_st  = 100;
      cuts_ph_min_st  = -0.032;
      cuts_ph_max_st  = 0.032;
      cuts_th_min_st  = -0.085;
      cuts_th_max_st  = 0.085;
      cuts_npeSum_st  = 2.0;
      cuts_etotTrack_st  = 0.7;

      cuts_delta_min_lt  = -8.0 ;
      cuts_delta_max_lt  = 8.0 ;
      cuts_ytar_min_lt  = -100;
      cuts_ytar_max_lt  = 100;
      cuts_ph_min_lt  = -0.032;
      cuts_ph_max_lt  = 0.032;
      cuts_th_min_lt  = -0.085;
      cuts_th_max_lt  = 0.085;
      cuts_npeSum_lt  = 2.0;
      cuts_etotTrack_lt  = 0.7;

      cuts_delta_min_st_mc = -8.0;
      cuts_delta_max_st_mc = 8.0;
      cuts_ytar_min_st_mc  = -100;   // not used
      cuts_ytar_max_st_mc  = 100;    // not used
      cuts_ph_min_st_mc    = -0.032;
      cuts_ph_max_st_mc    = 0.032;
      cuts_th_min_st_mc    = -0.085;
      cuts_th_max_st_mc    = 0.085;
      cuts_npeSum_st_mc    = 2.0;   // not used
      cuts_etotTrack_st_mc = 0.7;   // not used
        
      cuts_delta_min_lt_mc = -8.0;
      cuts_delta_max_lt_mc = 8.0;
      cuts_ytar_min_lt_mc  = -100;
      cuts_ytar_max_lt_mc  = 100;
      cuts_ph_min_lt_mc    = -0.032;
      cuts_ph_max_lt_mc    = 0.032;
      cuts_th_min_lt_mc    = -0.085;
      cuts_th_max_lt_mc    = 0.085;
      cuts_npeSum_lt_mc    = 2.0;
      cuts_etotTrack_lt_mc = 0.7;

    }
    else if (spectrometer=="SHMS"){
      cuts_delta_min_st  = -10.0 ;
      cuts_delta_max_st  = 22.0 ;
      cuts_ytar_min_st  = -1.5;
      cuts_ytar_max_st  = 1.5;
      cuts_ph_min_st  = -0.04;
      cuts_ph_max_st  = 0.04;
      cuts_th_min_st  = -0.07;
      cuts_th_max_st  = 0.07;
      cuts_npeSum_st  = 2.0;
      cuts_etotTrack_st  = 0.7;


      cuts_delta_min_lt  = -10.0 ;
      cuts_delta_max_lt  = 22.0 ;
      cuts_ytar_min_lt  = -0.75;
      cuts_ytar_max_lt  = 0.75;
      cuts_ph_min_lt  = -0.04;
      cuts_ph_max_lt  = 0.04;
      cuts_th_min_lt  = -0.07;
      cuts_th_max_lt  = 0.07;
      cuts_npeSum_lt  = 2.0;
      cuts_etotTrack_lt  = 0.7;

      cuts_delta_min_st_mc = -10.0;
      cuts_delta_max_st_mc = 22.0;
      cuts_ytar_min_st_mc = -1.5;
      cuts_ytar_max_st_mc = 1.5;
      cuts_ph_min_st_mc = -0.04;
      cuts_ph_max_st_mc = 0.04;
      cuts_th_min_st_mc = -0.07;
      cuts_th_max_st_mc = 0.07;
      cuts_npeSum_st_mc = 2.0;
      cuts_etotTrack_st_mc = 0.7;

      cuts_delta_min_lt_mc = -10.0;
      cuts_delta_max_lt_mc = 22.0;
      cuts_ytar_min_lt_mc = -0.75;
      cuts_ytar_max_lt_mc = 0.75;
      cuts_ph_min_lt_mc = -0.04;
      cuts_ph_max_lt_mc = 0.04;
      cuts_th_min_lt_mc = -0.07;
      cuts_th_max_lt_mc = 0.07;
      cuts_npeSum_lt_mc = 2.0;
      cuts_etotTrack_lt_mc = 0.7;

    }


    for (int i = 0; i < run_data_ST.size(); i++) { cout<< "Looking at run (ST): " << run_data_ST[i].run_number << endl; }

    //////////    ST   ///////////

    for (int i = 0; i < run_data_ST.size(); i++) {
      if (spectrometer=="test"){spectrometer="HMS";}
      //      cout<< "Looking at run : " << run_data_ST[i].run_number << endl;
    TFile *f_st_data = new TFile( Form("/work/smoran/xem2/data/%s/%s_replay_production_%d_-1.root" ,(const char*)spectrometer,(const char*)spectrometer_lowercase ,run_data_ST[i].run_number ));
    TTree *t_st_data = (TTree*) f_st_data->Get("T");
    Variables st_vars_data;

    setBranchAddresses_data(t_st_data, st_vars_data,spectrometer);
    Long64_t nentries_st_data = t_st_data->GetEntries();
    cout<<"nentries_data : " << nentries_st_data <<endl;


    TH1F *X_tmp   = (TH1F*)X_MC_LT->Clone();
    	for (int j = 0; j < nentries_st_data; j++){
    		t_st_data->GetEntry(j);
		
		
if(st_vars_data.ptardp >= cuts_delta_min_st && st_vars_data.ptardp <= cuts_delta_max_st 
   && st_vars_data.decal  > cuts_etotTrack_st
   && st_vars_data.ptarth > cuts_th_min_st    &&   st_vars_data.ptarth < cuts_th_max_st   
   && st_vars_data.ptarph > cuts_ph_min_st    &&   st_vars_data.ptarph < cuts_ph_max_st    
   && st_vars_data.npeSum > cuts_npeSum_st    
   && st_vars_data.ptary  < cuts_ytar_max_st  && st_vars_data.ptary  > cuts_ytar_min_st  /*&& st_vars_data.ptarx < 2.0*/
   && st_vars_data.AvgCurr >5.0 && st_vars_data.CurrentFlag ==1

  ) {
		
      double Eprime ;
      Eprime =  p_spec*(1+0.01*st_vars_data.ptardp);
      double thetarad;
      if (spectrometer=="HMS" || spectrometer=="test")      {thetarad= TMath::ACos((cos(theta_central) + st_vars_data.ptarph*sin(theta_central))/TMath::Sqrt(1. + st_vars_data.ptarph*st_vars_data.ptarph+st_vars_data.ptarth*st_vars_data.ptarth));      }
      else if (spectrometer=="SHMS"){thetarad= TMath::ACos((cos(theta_central) - st_vars_data.ptarph*sin(theta_central))/TMath::Sqrt(1. + st_vars_data.ptarph*st_vars_data.ptarph+st_vars_data.ptarth*st_vars_data.ptarth));      }
      
      double thetadeg = thetarad*(180.0/3.14);                   
      double weight_F2 , var_x_axis  , weight_hms_cer_eff;
      weight_hms_cer_eff = hmscereff_xem2(st_vars_data.ptardp , "HMS"); 

      if (var_dep=="delta")  { var_x_axis = st_vars_data.ptardp; }
      else if (var_dep=="xb"){ var_x_axis = deltaToX( st_vars_data.ptardp , p_spec , Einitial , angle);  }
      if (F2_option) { weight_F2 = weight_calculateF2( thetadeg ,  Eprime ,  Einitial)               ;}      else {  weight_F2 =1.0;}

      if (shift_collimator_opt =="00"){ X_tmp->Fill( var_x_axis*weight_hms_cer_eff ); }

      else {  if (abs(166.37*st_vars_data.ptarth) <= (11.645 -shift_collimator) && abs(st_vars_data.ptary + 166.37*st_vars_data.ptarph )<= (4.575 -shift_collimator) &&
		  (166.37*st_vars_data.ptarth) < (2.545 * (st_vars_data.ptary + 166.37*st_vars_data.ptarph) + (17.4625 - vertical_shift)) &&
		  (166.37*st_vars_data.ptarth) > ((2.545)*(st_vars_data.ptary + 166.37*st_vars_data.ptarph)  - (17.4625 - vertical_shift)) &&
		  (166.37*st_vars_data.ptarth) < ((-2.545)*(st_vars_data.ptary + 166.37*st_vars_data.ptarph) + (17.4625 - vertical_shift)) &&
		  (166.37*st_vars_data.ptarth) > ((-2.545)*(st_vars_data.ptary + 166.37*st_vars_data.ptarph) - (17.4625 - vertical_shift))     ) 
	        {  
	              X_tmp->Fill( var_x_axis  , weight_F2   );
		}       
            }

      		}
    	} // end cycle through ST entries

	auto result = calculateTotalEntriesAndError(X_tmp);
	double scaling_factor_run =  run_data_ST[i].PS / (run_data_ST[i].charge * run_data_ST[i].Eff_Fid * run_data_ST[i].Eff_Elec_LiveTime * run_data_ST[i].Eff_Comp_LiveTime );
	outFile_runs_ST << run_data_ST[i].run_number<< "\t"<< run_data_ST[i].trigger << "\t"<< result.first*scaling_factor_run << "\t" <<result.second*scaling_factor_run <<"\n" ;
    X_Data_ST->Add(X_tmp);
    } //end cycle through different ST runs

    //    cout<< "ST Yields: " << endl;
    // printHistogramContents(  X_Data_ST );

    outFile_runs_ST.close();


    //////////    LT   ///////////

      for (int i = 0; i < run_data_LT.size(); i++) {
	if (spectrometer=="test"){spectrometer="HMS";}    
    TFile *f_lt_data = new TFile( Form("/work/smoran/xem2/data/%s/%s_replay_production_%d_-1.root" ,(const char*)spectrometer,(const char*)spectrometer_lowercase , run_data_LT[i].run_number ));
    TTree *t_lt_data = (TTree*) f_lt_data->Get("T");
    Variables lt_vars_data;
    
    setBranchAddresses_data(t_lt_data, lt_vars_data,spectrometer);
    Long64_t nentries_lt_data = t_lt_data->GetEntries();
    TH1F *X_tmp   = (TH1F*)X_MC_LT->Clone();
    	for (int j = 0; j < nentries_lt_data; j++){
    		t_lt_data->GetEntry(j);
	    	if( lt_vars_data.ptardp >= cuts_delta_min_lt && lt_vars_data.ptardp <= cuts_delta_max_lt &&
		    lt_vars_data.decal  > cuts_etotTrack_lt && lt_vars_data.ptarth > cuts_th_min_lt    &&
		    lt_vars_data.ptarth < cuts_th_max_lt    && lt_vars_data.ptarph > cuts_ph_min_lt    &&
		    lt_vars_data.ptarph < cuts_ph_max_lt    && lt_vars_data.npeSum > cuts_npeSum_lt    &&
		    lt_vars_data.ptary  < cuts_ytar_max_lt  && lt_vars_data.ptary  > cuts_ytar_min_lt  && /* lt_vars_data.ptarx < 2.0 &&*/
		    lt_vars_data.AvgCurr >5.0  && lt_vars_data.CurrentFlag ==1
){
		    double Eprime =  p_spec*(1+0.01*lt_vars_data.ptardp);
		    double thetarad;
      		    if (spectrometer=="HMS" || spectrometer=="test")      {thetarad= TMath::ACos((cos(theta_central) + lt_vars_data.ptarph*sin(theta_central))/TMath::Sqrt(1. + lt_vars_data.ptarph*lt_vars_data.ptarph+lt_vars_data.ptarth*lt_vars_data.ptarth));      }
                    else if (spectrometer=="SHMS"){thetarad= TMath::ACos((cos(theta_central) - lt_vars_data.ptarph*sin(theta_central))/TMath::Sqrt(1. + lt_vars_data.ptarph*lt_vars_data.ptarph+lt_vars_data.ptarth*lt_vars_data.ptarth));      }
		    double thetadeg = thetarad*(180.0/3.14); 
		    double weight_F2  , var_x_axis , weight_hms_cer_eff;
		    weight_hms_cer_eff = hmscereff_xem2(lt_vars_data.ptardp , "HMS");


		    if (var_dep=="delta")  { var_x_axis = lt_vars_data.ptardp;  }
		    else if (var_dep=="xb"){ var_x_axis = deltaToX( lt_vars_data.ptardp , p_spec , Einitial , angle);  }
		    if (F2_option) { weight_F2 = weight_calculateF2( thetadeg ,  Eprime ,  Einitial)            ;}      else { weight_F2 =1.0;}

		if (shift_collimator_opt =="00"){ X_tmp->Fill( var_x_axis*weight_hms_cer_eff ); }
		else {  if (abs(166.37*lt_vars_data.ptarth) <= (11.645 -shift_collimator) && abs(lt_vars_data.ptary + 166.37*lt_vars_data.ptarph )<= (4.575 -shift_collimator) &&
			    (166.37*lt_vars_data.ptarth) < (2.545 * (lt_vars_data.ptary + 166.37*lt_vars_data.ptarph) + (17.4625 - vertical_shift)) &&
			    (166.37*lt_vars_data.ptarth) > ((2.545)*(lt_vars_data.ptary + 166.37*lt_vars_data.ptarph)  - (17.4625 - vertical_shift)) &&
			    (166.37*lt_vars_data.ptarth) < ((-2.545)*(lt_vars_data.ptary + 166.37*lt_vars_data.ptarph) + (17.4625 - vertical_shift)) &&
			    (166.37*lt_vars_data.ptarth) > ((-2.545)*(lt_vars_data.ptary + 166.37*lt_vars_data.ptarph) - (17.4625 - vertical_shift))     )
		    {
                      X_tmp->Fill( var_x_axis  , weight_F2   );
		    }
		}

      		}
    	}
        auto result = calculateTotalEntriesAndError(X_tmp);
        double scaling_factor_run =  run_data_LT[i].PS / (run_data_LT[i].charge * run_data_LT[i].Eff_Fid * run_data_LT[i].Eff_Elec_LiveTime * run_data_LT[i].Eff_Comp_LiveTime );
        outFile_runs_LT << run_data_LT[i].run_number<< "\t"<<run_data_LT[i].trigger<< "\t"<< result.first*scaling_factor_run << "\t" <<result.second*scaling_factor_run <<"\n" ;

    X_Data_LT->Add(X_tmp);
     }

      outFile_runs_LT.close();

      ////////   dummy   ////////

      for (int i = 0; i < run_data_LT_Dummy.size(); i++) {
	if (spectrometer=="test"){spectrometer="HMS";}
	TFile *f_lt_data_Dummy = new TFile( Form("/work/smoran/xem2/data/%s/%s_replay_production_%d_-1.root" ,(const char*)spectrometer,(const char*)spectrometer_lowercase, run_data_LT_Dummy[i].run_number ));
	TTree *t_lt_data_Dummy = (TTree*) f_lt_data_Dummy->Get("T");
	Variables lt_vars_data_Dummy;

	setBranchAddresses_data(t_lt_data_Dummy, lt_vars_data_Dummy,spectrometer);
	Long64_t nentries_lt_data_Dummy = t_lt_data_Dummy->GetEntries();
	TH1F *X_tmp   = (TH1F*)X_MC_LT->Clone();
        for (int j = 0; j < nentries_lt_data_Dummy; j++){
	  t_lt_data_Dummy->GetEntry(j);
	  if(
	     lt_vars_data_Dummy.ptardp >= cuts_delta_min_lt && lt_vars_data_Dummy.ptardp <= cuts_delta_max_lt 
	     && lt_vars_data_Dummy.decal  > cuts_etotTrack_lt 
	     && lt_vars_data_Dummy.ptarth > cuts_th_min_lt    &&  lt_vars_data_Dummy.ptarth < cuts_th_max_lt    
	     && lt_vars_data_Dummy.ptarph > cuts_ph_min_lt    &&  lt_vars_data_Dummy.ptarph < cuts_ph_max_lt   
	     && lt_vars_data_Dummy.npeSum > cuts_npeSum_lt    
	     && lt_vars_data_Dummy.ptary  < cuts_ytar_max_lt  && lt_vars_data_Dummy.ptary  > cuts_ytar_min_lt /*&& lt_vars_data_Dummy.ptarx <2.0 */
	     && lt_vars_data_Dummy.AvgCurr >5.0  && lt_vars_data_Dummy.CurrentFlag ==1
 ){
	     
	     double Eprime =  p_spec*(1+0.01*lt_vars_data_Dummy.ptardp);
	     double thetarad;
      	     if (spectrometer=="HMS"||spectrometer=="test")      {thetarad= TMath::ACos((cos(theta_central) + lt_vars_data_Dummy.ptarph*sin(theta_central))/TMath::Sqrt(1. + lt_vars_data_Dummy.ptarph*lt_vars_data_Dummy.ptarph+lt_vars_data_Dummy.ptarth*lt_vars_data_Dummy.ptarth));      }
             else if (spectrometer=="SHMS"){thetarad= TMath::ACos((cos(theta_central) - lt_vars_data_Dummy.ptarph*sin(theta_central))/TMath::Sqrt(1. + lt_vars_data_Dummy.ptarph*lt_vars_data_Dummy.ptarph+lt_vars_data_Dummy.ptarth*lt_vars_data_Dummy.ptarth));      }
	     double thetadeg = thetarad*(180.0/3.14); 
	     double weight_F2 , var_x_axis , weight_hms_cer_eff;

	     weight_hms_cer_eff = hmscereff_xem2(lt_vars_data_Dummy.ptardp , "HMS");

	     if (var_dep=="delta")  { var_x_axis = lt_vars_data_Dummy.ptardp;  }
	     else if (var_dep=="xb"){ var_x_axis = deltaToX( lt_vars_data_Dummy.ptardp , p_spec , Einitial , angle);  }
	     if (F2_option) {  weight_F2 = weight_calculateF2( thetadeg ,  Eprime ,  Einitial)               ;}	     else { weight_F2 =1.0;}
	     	     

	     if (shift_collimator_opt =="00"){ X_tmp->Fill( var_x_axis*weight_hms_cer_eff  ); }
	     else {  if (abs(166.37*lt_vars_data_Dummy.ptarth) <= (11.645 -shift_collimator) && abs(lt_vars_data_Dummy.ptary + 166.37*lt_vars_data_Dummy.ptarph )<= (4.575 -shift_collimator) &&
			 (166.37*lt_vars_data_Dummy.ptarth) < (2.545 * (lt_vars_data_Dummy.ptary + 166.37*lt_vars_data_Dummy.ptarph) + (17.4625 - vertical_shift)) &&
			 (166.37*lt_vars_data_Dummy.ptarth) > ((2.545)*(lt_vars_data_Dummy.ptary + 166.37*lt_vars_data_Dummy.ptarph)  - (17.4625 - vertical_shift)) &&
			 (166.37*lt_vars_data_Dummy.ptarth) < ((-2.545)*(lt_vars_data_Dummy.ptary + 166.37*lt_vars_data_Dummy.ptarph) + (17.4625 - vertical_shift)) &&
			 (166.37*lt_vars_data_Dummy.ptarth) > ((-2.545)*(lt_vars_data_Dummy.ptary + 166.37*lt_vars_data_Dummy.ptarph) - (17.4625 - vertical_shift))     )
		 {
		   X_tmp->Fill( var_x_axis  , weight_F2   );
		 }
	     }


	  }
        }

        auto result = calculateTotalEntriesAndError(X_tmp);
        double scaling_factor_run =  run_data_LT_Dummy[i].PS / (run_data_LT_Dummy[i].charge * run_data_LT_Dummy[i].Eff_Fid * run_data_LT_Dummy[i].Eff_Elec_LiveTime * run_data_LT_Dummy[i].Eff_Comp_LiveTime );
        outFile_runs_LT_Dummy << run_data_LT_Dummy[i].run_number<<"\t"<<run_data_LT_Dummy[i].trigger<< "\t"<< result.first*scaling_factor_run << "\t" <<result.second*scaling_factor_run <<"\n" ;

	X_Data_LT_Dummy->Add(X_tmp);
      } //end cycle through dummy 

      outFile_runs_LT_Dummy.close();



      ////////////    contamination subtraction for Ca48 , BC10, and BC11   /////////////////////////



      if (target_input=="Ca48" || target_input=="BC10" || target_input=="BC11") {  


	//////////     tmp     //////////////
                                                                                                                                                                                                                                                                                                                                                                         
	for (int i = 0; i < run_data_ST_tmp.size(); i++) {
	  if (spectrometer=="test"){spectrometer="HMS";}
	  TFile *f_st_data_tmp = new TFile( Form("/work/smoran/xem2/data/%s/%s_replay_production_%d_-1.root" ,(const char*)spectrometer,(const char*)spectrometer_lowercase, run_data_ST_tmp[i].run_number ));
	  TTree *t_st_data_tmp = (TTree*) f_st_data_tmp->Get("T");
	  Variables st_vars_data_tmp;

	  setBranchAddresses_data(t_st_data_tmp, st_vars_data_tmp,spectrometer);
	  Long64_t nentries_st_data_tmp = t_st_data_tmp->GetEntries();
	  TH1F *X_tmp   = (TH1F*)X_MC_ST->Clone();
	  for (int j = 0; j < nentries_st_data_tmp; j++){
	    t_st_data_tmp->GetEntry(j);
	    if(
             st_vars_data_tmp.ptardp >= cuts_delta_min_st && st_vars_data_tmp.ptardp <= cuts_delta_max_st
             && st_vars_data_tmp.decal  > cuts_etotTrack_st
             && st_vars_data_tmp.ptarth > cuts_th_min_st    &&  st_vars_data_tmp.ptarth < cuts_th_max_st
             && st_vars_data_tmp.ptarph > cuts_ph_min_st    &&  st_vars_data_tmp.ptarph < cuts_ph_max_st
             && st_vars_data_tmp.npeSum > cuts_npeSum_st
             && st_vars_data_tmp.ptary  < cuts_ytar_max_st  && st_vars_data_tmp.ptary  > cuts_ytar_min_st /*&& st_vars_data_tmp.ptarx <2.0*/
             && st_vars_data_tmp.AvgCurr >5.0  && st_vars_data_tmp.CurrentFlag ==1
	       ){

	      double Eprime =  p_spec*(1+0.01*st_vars_data_tmp.ptardp);
	      double thetarad;
	      if (spectrometer=="HMS"||spectrometer=="test")      {thetarad= TMath::ACos((cos(theta_central) + st_vars_data_tmp.ptarph*sin(theta_central))/TMath::Sqrt(1. + st_vars_data_tmp.ptarph*st_vars_data_tmp.ptarph+st_vars_data_tmp.ptarth*st_vars_data_tmp.ptarth));      }
	      else if (spectrometer=="SHMS"){thetarad= TMath::ACos((cos(theta_central) - st_vars_data_tmp.ptarph*sin(theta_central))/TMath::Sqrt(1. + st_vars_data_tmp.ptarph*st_vars_data_tmp.ptarph+st_vars_data_tmp.ptarth*st_vars_data_tmp.ptarth));      }
	      double thetadeg = thetarad*(180.0/3.14);
	      double weight_F2 , var_x_axis;
	      if (var_dep=="delta")  { var_x_axis = st_vars_data_tmp.ptardp;  }
	      else if (var_dep=="xb"){ var_x_axis = deltaToX( st_vars_data_tmp.ptardp , p_spec , Einitial , angle);  }
	      if (F2_option) {  weight_F2 = weight_calculateF2( thetadeg ,  Eprime ,  Einitial)               ;}
	      else { weight_F2 =1.0;}

	      X_tmp->Fill(var_x_axis, weight_F2);


	      if (shift_collimator_opt =="00"){ X_tmp->Fill( var_x_axis  ); }
	      else {  if (abs(166.37*st_vars_data_tmp.ptarth) <= (11.645 -shift_collimator) && abs(st_vars_data_tmp.ptary + 166.37*st_vars_data_tmp.ptarph )<= (4.575 -shift_collimator) &&
			  (166.37*st_vars_data_tmp.ptarth) < (2.545 * (st_vars_data_tmp.ptary + 166.37*st_vars_data_tmp.ptarph) + (17.4625 - vertical_shift)) &&
			  (166.37*st_vars_data_tmp.ptarth) > ((2.545)*(st_vars_data_tmp.ptary + 166.37*st_vars_data_tmp.ptarph)  - (17.4625 - vertical_shift)) &&
			  (166.37*st_vars_data_tmp.ptarth) < ((-2.545)*(st_vars_data_tmp.ptary + 166.37*st_vars_data_tmp.ptarph) + (17.4625 - vertical_shift)) &&
			  (166.37*st_vars_data_tmp.ptarth) > ((-2.545)*(st_vars_data_tmp.ptary + 166.37*st_vars_data_tmp.ptarph) - (17.4625 - vertical_shift))     )
		  {
		    X_tmp->Fill( var_x_axis  , weight_F2   );
		  }
	      }

	    }
	  }

	  auto result = calculateTotalEntriesAndError(X_tmp);
	  double scaling_factor_run =  run_data_ST_tmp[i].PS / (run_data_ST_tmp[i].charge * run_data_ST_tmp[i].Eff_Fid * run_data_ST_tmp[i].Eff_Elec_LiveTime * run_data_ST_tmp[i].Eff_Comp_LiveTime );
	  outFile_runs_ST_tmp << run_data_ST_tmp[i].run_number<<"\t"<<run_data_ST_tmp[i].trigger<< "\t"<< result.first*scaling_factor_run << "\t" <<result.second*scaling_factor_run <<"\n" ;

	  X_Data_ST_tmp->Add(X_tmp);
	} //end cycle through tmp                                                                                                                                                                                                                                                                                                                                                           
	outFile_runs_ST_tmp.close();
      
      }


      std::cout<< "Production scale factors for ST, LT , and dummy: "<< production_scale_ST << " , "<< production_scale_LT << " , " << production_scale_LT_Dummy << std::endl;

      X_Data_ST->Scale(production_scale_ST);
      X_Data_LT->Scale(production_scale_LT);
      X_Data_LT_Dummy->Scale(production_scale_LT_Dummy);

      if (target_input=="Ca48" || target_input=="BC10" || target_input=="BC11") {  X_Data_ST_tmp->Scale(production_scale_ST_tmp);    }

      writeHistogramToTextFile( X_Data_ST , p_spec , Einitial, angle , out_name + "_CNY_ST_TABLE.txt" ) ;
      writeHistogramToTextFile( X_Data_LT , p_spec , Einitial, angle , out_name + "_CNY_LT_TABLE.txt" ) ;

      float R , thickness_LD2 , thickness_dummy;
      float LD2_entrance = 0.168;                        //mm   
      float LD2_exit     = 0.2024;                       // mm  
      float LD2_average  = (LD2_entrance + LD2_exit)/2;  //take average of the entrance and exit to get thickness, but it still has the mm units. Needs to match  thickness_dummy units (g/(cm^2)) 
      thickness_LD2      = LD2_average*(0.1)*(2.7);      // 0.1 converts mm to cm. 2.7 g/(cm^3) is the density of Aluminium, converts cm to g/(cm^2), it now matches the units of thickness_dummy units (g/(cm^2) 
      float dummy_upstream   = 0.240;       //  g/(cm^2)                    
      float dummy_downstream = 0.236;       //  g/(cm^2) 
      thickness_dummy = (dummy_upstream + dummy_downstream)/2;   // g/(cm^2), average of dummy upstream and downstream                                           
      R =  (thickness_LD2 )/ thickness_dummy; // 0.996 cryo target contraction factor NOT USED FOR DUMMY
      //X_Data_LT_Dummy->Scale(R); we scale later, just in case we do CSB

      double R_tmp ; 
      if (target_input=="Ca48") {  R_tmp = 0.113 ; X_Data_ST_tmp->Scale(R_tmp); }
      if (target_input=="BC10") {  R_tmp = 0.231 ; X_Data_ST_tmp->Scale(R_tmp); }
      if (target_input=="BC11") {  R_tmp = 0.236 ; X_Data_ST_tmp->Scale(R_tmp); }

      ////////////////////  
      /////   SIM    /////  
      ////////////////////                                                                                                                                                     
      std::vector<double> hsdp_values_ST;
      std::vector<double> weight_values_ST;
      std::vector<double> hsdp_values_LT;
      std::vector<double> weight_values_LT;

      cout<< "SIM LT:" <<endl;
      int weight_zero_LT_count=0;
      int progress_step_LT = nentries_lt_mc / 100.0;
      for (int j=1; j<nentries_lt_mc; j++) {
	t_lt_mc->GetEntry(j);


	if (lt_vars_mc.stop_id==0.0 /*&& lt_vars_mc.hpsyptar < cuts_ph_max_lt_mc    && lt_vars_mc.hpsyptar > cuts_ph_min_lt_mc    &&
	                               lt_vars_mc.hpsxptar < cuts_th_max_lt_mc    && lt_vars_mc.hpsxptar > cuts_th_min_lt_mc     &&
				       lt_vars_mc.hpsytar  < cuts_ytar_max_lt_mc  && lt_vars_mc.hpsytar  > cuts_ytar_min_lt_mc*/ ){

	  double Eprime , Eprimei;
	  Eprime  =  p_spec*(1+0.01*lt_vars_mc.hsdp) ;
	  Eprimei =  p_spec*(1+0.01*lt_vars_mc.hsdpi);
	  X_MC_LT_Unweighted->Fill( lt_vars_mc.hsdp) ;
	  double thetarad , thetaradi;
	  if (spectrometer=="HMS"||spectrometer=="test"){
	   thetarad = TMath::ACos((cos(theta_central) + lt_vars_mc.hpsyptar*sin(theta_central))/TMath::Sqrt( 1. + lt_vars_mc.hpsxptar*lt_vars_mc.hpsxptar  +lt_vars_mc.hpsyptar*lt_vars_mc.hpsyptar  ));
	   thetaradi= TMath::ACos((cos(theta_central) + lt_vars_mc.hpsyptari*sin(theta_central))/TMath::Sqrt(1. + lt_vars_mc.hpsxptari*lt_vars_mc.hpsxptari+lt_vars_mc.hpsyptari*lt_vars_mc.hpsyptari));
	  }
	  
	  else if (spectrometer=="SHMS"){
	   thetarad = TMath::ACos((cos(theta_central) - lt_vars_mc.hpsyptar*sin(theta_central))/TMath::Sqrt(1. + lt_vars_mc.hpsxptar*lt_vars_mc.hpsxptar+lt_vars_mc.hpsyptar*lt_vars_mc.hpsyptar));
	   thetaradi= TMath::ACos((cos(theta_central) - lt_vars_mc.hpsyptari*sin(theta_central))/TMath::Sqrt(1. + lt_vars_mc.hpsxptari*lt_vars_mc.hpsxptari+lt_vars_mc.hpsyptari*lt_vars_mc.hpsyptari));
	  }	  
	  
	  double thetadeg  = thetarad*(180.0/3.14);
	  double thetadegi = thetaradi*(180.0/3.14);
	  thetadeg  = findNearestValue(thetadeg, angle);
	  thetadegi = findNearestValue(thetadegi, angle);
	  double weight = Interpolate1DFromGraph2DErrors(gr2D_SigmaRad_LT, Eprimei, thetadegi );  
	  double weight_F2;
	  if (F2_option) { weight_F2 = weight_calculateF2( thetadeg ,  Eprime ,  Einitial)               ;}
	  else {weight_F2 = 1.0;}
	  
	  double weight_ytarCorr , weight_MC_JacobianCorr , weight_delta;
	  
	  if (ytar_corr_ON){weight_ytarCorr = weight_ytar_corr( lt_vars_mc.hpsytar , spectrometer);}
	  else weight_ytarCorr=1.0;
	  if (MC_Jacobian_corr_ON){weight_MC_JacobianCorr = weight_MC_Jacobian_corr ( lt_vars_mc.hpsxptari   , lt_vars_mc.hpsyptari,  spectrometer  );}
	  else weight_MC_JacobianCorr=1.0;
	  if (delta_corr_ON) {weight_delta = weight_delta_corr_Dave(lt_vars_mc.hsdpi , spectrometer );}
	  else weight_delta = 1.0;
	  double var_x_axis;
	  if      (var_dep=="delta")  { var_x_axis = lt_vars_mc.hsdp;  }
	  else if (var_dep=="xb"   )  { var_x_axis = deltaToX( lt_vars_mc.hsdp , p_spec , Einitial , angle);  }
	  
	  if (  lt_vars_mc.hsdp <= cuts_delta_max_lt_mc && lt_vars_mc.hsdp >= cuts_delta_min_lt_mc     ){	   


	    if (shift_collimator_opt =="00"){ X_MC_LT->Fill( var_x_axis) ; }//,weight*weight_ytarCorr*weight_MC_JacobianCorr*weight_delta ); }
	    else {  if (abs(166.37*lt_vars_mc.hpsxptar) <= (11.645 -shift_collimator) && abs(lt_vars_mc.hpsytar + 166.37*lt_vars_mc.hpsyptar )<= (4.575 -shift_collimator) &&
			(166.37*lt_vars_mc.hpsxptar) < (2.545 * (lt_vars_mc.hpsytar + 166.37*lt_vars_mc.hpsyptar) + (17.4625 - vertical_shift)) &&
			(166.37*lt_vars_mc.hpsxptar) > ((2.545)*(lt_vars_mc.hpsytar + 166.37*lt_vars_mc.hpsyptar)  - (17.4625 - vertical_shift)) &&
			(166.37*lt_vars_mc.hpsxptar) < ((-2.545)*(lt_vars_mc.hpsytar + 166.37*lt_vars_mc.hpsyptar) + (17.4625 - vertical_shift)) &&
			(166.37*lt_vars_mc.hpsxptar) > ((-2.545)*(lt_vars_mc.hpsytar + 166.37*lt_vars_mc.hpsyptar) - (17.4625 - vertical_shift))     )
		{
		  X_MC_LT->Fill( var_x_axis,weight*weight_F2*weight_ytarCorr*weight_MC_JacobianCorr*weight_delta );
		}
	    }

	  } 
	  
	  hsdp_values_LT.push_back(lt_vars_mc.hsdp);
	  weight_values_LT.push_back(weight);
	  if (weight==0) {  cout <<"E prime: " << Eprime << " , angle: "<< angle << " , delta is "<< lt_vars_mc.hsdp << endl; weight_zero_LT_count=weight_zero_LT_count+1;}
	}

	if (j % progress_step_LT == 0) {
	  int progress = (j * 100) / nentries_lt_mc;
	  cout << "\rProgress: " << progress << "% completed , for LT" << flush;
	}

      }
      cout << "There were "<< weight_zero_LT_count << " zeroes in LT." <<endl;



      X_MC_LT->Scale(mc_scale_LT);
      X_MC_LT_Unweighted->Scale(mc_scale_LT);

      writeHistogramToTextFile( X_MC_LT , p_spec , Einitial, angle , out_name + "_LT_MC_xs_unwgted_TABLE.txt" ) ;




      cout<< "SIM ST:" <<endl;
      int weight_zero_ST_count=0;
      int progress_step_ST = nentries_st_mc / 100;
      for (int j=1; j<nentries_st_mc; j++) {
	t_st_mc->GetEntry(j);
	if (   st_vars_mc.stop_id==0.0  /*&& 
				          st_vars_mc.hpsyptar < cuts_ph_max_st_mc    && st_vars_mc.hpsyptar > cuts_ph_min_st_mc     &&
            				  st_vars_mc.hpsxptar < cuts_th_max_st_mc    && st_vars_mc.hpsxptar > cuts_th_min_st_mc     &&
					  st_vars_mc.hpsytar  < cuts_ytar_max_st_mc  && st_vars_mc.hpsytar  > cuts_ytar_min_st_mc  */   ){
	  double Eprime, Eprimei ;
	  Eprime  =  p_spec*(1+0.01*st_vars_mc.hsdp );
	  Eprimei =  p_spec*(1+0.01*st_vars_mc.hsdpi);
	  X_MC_ST_Unweighted->Fill( st_vars_mc.hsdp);
	  double thetarad, thetaradi;
	  if (spectrometer=="HMS" || spectrometer=="test"){
	    thetarad = TMath::ACos((cos(theta_central) + st_vars_mc.hpsyptar*sin(theta_central))/TMath::Sqrt(1. + st_vars_mc.hpsxptar*st_vars_mc.hpsxptar+st_vars_mc.hpsyptar*st_vars_mc.hpsyptar));
	    thetaradi= TMath::ACos((cos(theta_central) + st_vars_mc.hpsyptari*sin(theta_central))/TMath::Sqrt(1. + st_vars_mc.hpsxptari*st_vars_mc.hpsxptari+st_vars_mc.hpsyptari*st_vars_mc.hpsyptari));
}
	else if (spectrometer=="SHMS"){
	 thetarad = TMath::ACos((cos(theta_central) - st_vars_mc.hpsyptar*sin(theta_central))/TMath::Sqrt(1. + st_vars_mc.hpsxptar*st_vars_mc.hpsxptar+st_vars_mc.hpsyptar*st_vars_mc.hpsyptar));
	 thetaradi= TMath::ACos((cos(theta_central) - st_vars_mc.hpsyptari*sin(theta_central))/TMath::Sqrt(1. + st_vars_mc.hpsxptari*st_vars_mc.hpsxptari+st_vars_mc.hpsyptari*st_vars_mc.hpsyptari));
	}
	  double thetadeg  = thetarad*(180.0/3.14);
	  double thetadegi = thetaradi*(180.0/3.14);
	  thetadeg  = findNearestValue(thetadeg, angle);
	  thetadegi = findNearestValue(thetadegi, angle);
	  double weight = Interpolate1DFromGraph2DErrors(gr2D_SigmaRad_ST, Eprimei, thetadegi );  // using theta calculated
	  // applying the CC to the simulation instead of the data
	  //      double weight_CC = Interpolate1DFromGraph2DErrors(gr2D_CoulombCorrection_ST, Eprime, thetadeg );                    
	  double weight_F2;
	  if (F2_option) { weight_F2 = weight_calculateF2( thetadeg ,  Eprime ,  Einitial)               ;}
	  else {weight_F2 = 1.0;}
	  
	  double weight_ytarCorr , weight_MC_JacobianCorr , weight_delta;
	  
	  if (ytar_corr_ON){weight_ytarCorr = weight_ytar_corr( st_vars_mc.hpsytar , spectrometer);}
	  else weight_ytarCorr=1.0;
	  if (MC_Jacobian_corr_ON){weight_MC_JacobianCorr = weight_MC_Jacobian_corr ( st_vars_mc.hpsxptari   , st_vars_mc.hpsyptari,  spectrometer  );}
	  else weight_MC_JacobianCorr=1.0;
	  if (delta_corr_ON) {weight_delta = weight_delta_corr_Dave(st_vars_mc.hsdpi , spectrometer );}
	  else weight_delta = 1.0;
          double var_x_axis;
          if      (var_dep=="delta")  { var_x_axis = st_vars_mc.hsdp;  }
          else if (var_dep=="xb"   )  { var_x_axis = deltaToX( st_vars_mc.hsdp , p_spec , Einitial , angle);  }
	  
	  if ( st_vars_mc.hsdp <= cuts_delta_max_st_mc && st_vars_mc.hsdp >= cuts_delta_min_st_mc  ){




	    //	    X_MC_ST->Fill( var_x_axis , weight*weight_F2*weight_ytarCorr*weight_MC_JacobianCorr*weight_delta );


	    if (shift_collimator_opt =="00"){ X_MC_ST->Fill( var_x_axis); }  //  ,weight*weight_ytarCorr*weight_MC_JacobianCorr*weight_delta ); }
	    else {  if (abs(166.37*st_vars_mc.hpsxptar) <= (11.645 -shift_collimator) && abs(st_vars_mc.hpsytar + 166.37*st_vars_mc.hpsyptar )<= (4.575 -shift_collimator) &&
			(166.37*st_vars_mc.hpsxptar) < (2.545 * (st_vars_mc.hpsytar + 166.37*st_vars_mc.hpsyptar) + (17.4625 - vertical_shift)) &&
			(166.37*st_vars_mc.hpsxptar) > ((2.545)*(st_vars_mc.hpsytar + 166.37*st_vars_mc.hpsyptar)  - (17.4625 - vertical_shift)) &&
			(166.37*st_vars_mc.hpsxptar) < ((-2.545)*(st_vars_mc.hpsytar + 166.37*st_vars_mc.hpsyptar) + (17.4625 - vertical_shift)) &&
			(166.37*st_vars_mc.hpsxptar) > ((-2.545)*(st_vars_mc.hpsytar + 166.37*st_vars_mc.hpsyptar) - (17.4625 - vertical_shift))     )
		{
		  X_MC_ST->Fill( var_x_axis,weight*weight_F2*weight_ytarCorr*weight_MC_JacobianCorr*weight_delta );
		}
	    }


	  }
	  


	  hsdp_values_ST.push_back(st_vars_mc.hsdp);
	  weight_values_ST.push_back(weight);

	  if (weight==0){  cout <<" weight is zero ; E prime: " << Eprime << " , angle: "<< angle << " , delta is "<< st_vars_mc.hsdp << endl; weight_zero_ST_count=weight_zero_ST_count+1;}
	  //if (weight_zero_ST_count>3) { cout<< "exiting the code, because weight_zero_ST_count>3"<<endl;  exit(0);}
	}

	if (j % progress_step_ST == 0) {
	  int progress = (j * 100) / nentries_st_mc;
	  cout << "\rProgress: " << progress << "% completed , for ST" << flush;
	}
      }
      cout << "There were "<< weight_zero_ST_count << " zeroes in ST." <<endl;
      

      X_MC_ST->Scale(mc_scale_ST);
      X_MC_ST_Unweighted->Scale(mc_scale_ST);


      writeHistogramToTextFile( X_MC_ST , p_spec , Einitial, angle , out_name + "_ST_MC_xs_unwgted_TABLE.txt" ) ;
      writeHistogramToTextFile( X_MC_LT , p_spec , Einitial, angle , out_name + "_LT_MC_xs_unwgted_TABLE.txt" ) ;




      TCanvas *cc1;
      PlotWithRatio(cc1, X_Data_ST, X_Data_LT, X_MC_ST, X_MC_LT,X_MC_ST_Unweighted ,X_MC_LT_Unweighted , out_name);
      delete cc1; 
      plotHistograms_Dummy(X_Data_LT, X_Data_LT_Dummy, out_name);
      if (target_input=="Ca48" || target_input=="BC10" || target_input=="BC11") {plotHistograms_Contamination(X_Data_ST, X_Data_ST_tmp, out_name);}
      plottingWeightFactors( hsdp_values_ST, weight_values_ST ,hsdp_values_LT, weight_values_LT ,  out_name);

      cout<< "UNWEIGHTED:::\n\n\n";
      for (int i=0 ; i < 20 ; i++){  cout<<X_MC_ST_Unweighted->GetBinContent(1+i)<<endl;  }

  // DATA / MC , both properly scaled.

      if (target_input=="Ca48" || target_input=="BC10" || target_input=="BC11") {      X_Data_ST->Add(X_Data_ST_tmp, -1); }  



      // APPLYING CSB CORRECTION
      TH1F* h_X_CSB_corr_ST_dp       = new TH1F("h_X_CSB_corr_ST", "h_X_CSB_corr_ST", X_Data_ST->GetNbinsX(), X_Data_ST->GetXaxis()->GetXmin(), X_Data_ST->GetXaxis()->GetXmax());
      TH1F* h_X_CSB_corr_LT_dp       = new TH1F("h_X_CSB_corr_LT", "h_X_CSB_corr_LT", X_Data_LT->GetNbinsX(), X_Data_LT->GetXaxis()->GetXmin(), X_Data_LT->GetXaxis()->GetXmax());
      TH1F* h_X_CSB_corr_LT_Dummy_dp = new TH1F("h_X_CSB_corr_LT_Dummy", "h_X_CSB_corr_LT_Dummy", X_Data_LT_Dummy->GetNbinsX(), X_Data_LT_Dummy->GetXaxis()->GetXmin(), X_Data_LT_Dummy->GetXaxis()->GetXmax());

      for (int i = 1; i <= X_Data_ST->GetNbinsX(); ++i) {
	double binContent        = X_Data_ST->GetBinContent(i);
	double binError          = X_Data_ST->GetBinError(i);
	double binCenter         = X_Data_ST->GetBinCenter(i);
	double Eprime            = p_spec*(1+0.01*binCenter)  ;
	double theta_central_deg = theta_central * 180.0/3.141592 ;
	double weight            = weight_CSB_corr(target_input, spectrometer_option, angle, Eprime) ;
	double updatedValue      = binContent * weight;
	cout << "for ST , CSB corr: weight is " << weight <<endl;
	h_X_CSB_corr_ST_dp->SetBinContent(i, updatedValue);
	h_X_CSB_corr_ST_dp->SetBinError(i, binError * weight);

      }

      for (int i = 1; i <= X_Data_LT->GetNbinsX(); ++i) {
	double binContent        = X_Data_LT->GetBinContent(i);
	double binError          = X_Data_LT->GetBinError(i);
	double binCenter         = X_Data_LT->GetBinCenter(i);
	double Eprime            = p_spec*(1+0.01*binCenter)  ;
	double theta_central_deg = theta_central * 180.0/3.141592 ;
	double weight            = weight_CSB_corr( "LD2", spectrometer_option, angle, Eprime) ;
	double updatedValue      = binContent * weight;
	h_X_CSB_corr_LT_dp->SetBinContent(i, updatedValue);
        h_X_CSB_corr_LT_dp->SetBinError(i, binError * weight);

      }


      for (int i = 1; i <= X_Data_LT_Dummy->GetNbinsX(); ++i) {
	double binContent        = X_Data_LT_Dummy->GetBinContent(i);
	double binError          = X_Data_LT_Dummy->GetBinError(i);
	double binCenter         = X_Data_LT_Dummy->GetBinCenter(i);
	double Eprime            = p_spec*(1+0.01*binCenter)  ;
	double theta_central_deg = theta_central * 180.0/3.141592 ;
	double weight            = weight_CSB_corr( "Dummy", spectrometer_option, angle, Eprime) ;
	double updatedValue      = binContent * weight;
        h_X_CSB_corr_LT_Dummy_dp->SetBinContent(i, updatedValue);
        h_X_CSB_corr_LT_Dummy_dp->SetBinError(i, binError * weight);

      }


	if (CSB_corr_ON) {
	  X_Data_ST->Reset();
	  X_Data_LT->Reset();
	  X_Data_LT_Dummy->Reset();	

	  X_Data_ST->Add(h_X_CSB_corr_ST_dp);
	  X_Data_LT->Add(h_X_CSB_corr_LT_dp);
	  X_Data_LT_Dummy->Add(h_X_CSB_corr_LT_Dummy_dp);  
	}

	writeHistogramToTextFile( X_Data_ST , p_spec , Einitial, angle , out_name + "_CNY_ST_CSB_corr_TABLE.txt" ) ;
	writeHistogramToTextFile( X_Data_LT , p_spec , Einitial, angle , out_name + "_CNY_LT_CSB_corr_TABLE.txt" ) ;


  X_Data_LT_Dummy->Scale(R);
  X_Data_ST->Divide(X_MC_ST);
  X_Data_LT->Add(X_Data_LT_Dummy, -1); // performing dummy subtraction
  X_Data_LT->Divide(X_MC_LT);

  writeHistogramToTextFile( X_MC_ST , p_spec , Einitial, angle , out_name + "_ST_MC_xs_unwgted_TABLE.txt" ) ;
  writeHistogramToTextFile( X_MC_LT , p_spec , Einitial, angle , out_name + "_LT_MC_xs_unwgted_TABLE.txt" ) ;

  writeHistogramToTextFile( X_Data_ST , p_spec , Einitial, angle , out_name + "_ST_xs_unwgted_TABLE.txt" ) ;
  writeHistogramToTextFile( X_Data_LT , p_spec , Einitial, angle , out_name + "_LT_xs_unwgted_TABLE.txt" ) ;

  std::vector<double> hsdp_values_ST_Born;
  std::vector<double> weight_values_ST_Born;
  std::vector<double> hsdp_values_LT_Born;
  std::vector<double> weight_values_LT_Born;


  TH1F* h_X_weighted_ST_dp = new TH1F("h_X_weighted_ST", "h_X_weighted_ST", X_Data_ST->GetNbinsX(), X_Data_ST->GetXaxis()->GetXmin(), X_Data_ST->GetXaxis()->GetXmax());     
  TH1F* h_X_weighted_LT_dp = new TH1F("h_X_weighted_LT", "h_X_weighted_LT", X_Data_LT->GetNbinsX(), X_Data_LT->GetXaxis()->GetXmin(), X_Data_LT->GetXaxis()->GetXmax());     
  // Loop over all bins in the original histogram                                                                                                                            
  for (int i = 1; i <= X_Data_ST->GetNbinsX(); ++i) {                             
    double binContent        = X_Data_ST->GetBinContent(i);
    double binError          = X_Data_ST->GetBinError(i);
    double binCenter         = X_Data_ST->GetBinCenter(i); 
    double Eprime            = p_spec*(1+0.01*binCenter)  ;
    double theta_central_deg = theta_central * 180.0/3.141592 ;                                                                                                   
    double weight            = gr2D_SigmaBorn_ST->Interpolate(Eprime,theta_central_deg); 
    double updatedValue      = binContent * weight;                                                                                                                
    h_X_weighted_ST_dp->SetBinContent(i, updatedValue);                            
    h_X_weighted_ST_dp->SetBinError(i, binError * weight); 

    hsdp_values_ST_Born.push_back(binCenter);
    weight_values_ST_Born.push_back(weight);

                                                                                                         
  }                                                                                                                                                                                     
 
  for (int i = 1; i <= X_Data_LT->GetNbinsX(); ++i) {
    double binContent        = X_Data_LT->GetBinContent(i);      
    double binError          = X_Data_LT->GetBinError(i);                        
    double binCenter         = X_Data_LT->GetBinCenter(i);  
    double Eprime            = p_spec*(1+0.01*binCenter)  ;
    double theta_central_deg = theta_central * 180.0/3.141592 ;    
    double weight            = gr2D_SigmaBorn_LT->Interpolate(Eprime,theta_central_deg); 
    double updatedValue      = binContent * weight;                                      
    h_X_weighted_LT_dp->SetBinContent(i, updatedValue);       
    h_X_weighted_LT_dp->SetBinError(i, binError * weight);    

    hsdp_values_LT_Born.push_back(binCenter);
    weight_values_LT_Born.push_back(weight);

  }

  writeHistogramToTextFile( h_X_weighted_ST_dp , p_spec , Einitial, angle , out_name + "_ST_xs_TABLE.txt" ) ;
  writeHistogramToTextFile( h_X_weighted_LT_dp , p_spec , Einitial, angle , out_name + "_LT_xs_TABLE.txt" ) ;
                                                                                                                                                                    
  h_X_weighted_ST_dp->Divide(h_X_weighted_LT_dp);   
  h_X_weighted_ST_dp->Scale(AM_LT/AM_ST);   

  plottingWeightFactors_Born( hsdp_values_ST_Born, weight_values_ST_Born ,hsdp_values_LT_Born, weight_values_LT_Born ,  out_name);

  for (int i = 1; i <= h_X_weighted_ST_dp->GetNbinsX(); ++i){

    cout << h_X_weighted_ST_dp->GetBinContent(i) << " +- " << h_X_weighted_ST_dp->GetBinError(i) << endl;
}

  TCanvas *c1_dp = new TCanvas("c1_dp", "Cross-Sections as a function of delta", 500, 500);
  h_X_weighted_ST_dp->Draw();
  h_X_weighted_ST_dp->SetLineColor(kRed);
  h_X_weighted_ST_dp->SetMarkerColor(kRed); 
  h_X_weighted_ST_dp->SetMinimum(0);        
  h_X_weighted_ST_dp->SetMaximum(1.7);      
  c1_dp->SetGrid(); 
  printHistogramContents(h_X_weighted_ST_dp);  
  c1_dp->Print(out_name + "_DELTA.pdf");         
  c1_dp->Print(out_name + "_DELTA.root"); 

  writeHistogramToTextFile( h_X_weighted_ST_dp , p_spec , Einitial, angle , out_name + "_TABLE.txt" ) ;


}


int main(int argc, char** argv) {
  if (argc != 10) {
    std::cerr << "Usage: " << argv[0] << " <target_input> <angle_input> <momentum_spec_input> <spectrometer> <ytar_corr> <MCJacobian_corr> <delta_corr> <CSB_corr> <Shift_Collimator>" << std::endl;
    std::cerr << "Example: " << argv[0] << " Ca12 20 2p420 HMS ON ON ON ON 00" << std::endl;
    std::cerr << "Example of a TEST case: " << argv[0] << " Ca12 20 5p870 test ON ON ON ON 00" << std::endl;
    return 1;
  }

  cross_Sections_updated_workingProcess(argv[1], argv[2], argv[3], argv[4] , argv[5], argv[6], argv[7],argv[8] , argv[9] );
  return 0;
}
 
// functions


void ImportRadcor(vector<double> &v1,vector<double> &v2,vector<double> &v3,vector<double> &v4,vector<double> &v5,vector<double> &v6,vector<double> &v7,vector<double> &v8,vector<double> &v9, vector<double> &v10,vector<double> &v11,vector<double> &v12,vector<double> &v13, const char * filename){

  double i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16;
  ifstream infile(filename);

  if(infile.fail()){
    cout << "Cannot open the file: " << filename << endl;
    exit(1);
  }
  else{
    while(!infile.eof()){
      infile >> i1 >> i2 >> i3 >> i4>>i5>>i6>>i7>>i8>>i9>>i10>>i11>>i12>>i13;

      v1.push_back (i1);
      v2.push_back (i2);
      v3.push_back (i3);
      v4.push_back (i4);
      v5.push_back (i5);
      v6.push_back (i6);
      v7.push_back (i7);
      v8.push_back (i8);
      v9.push_back (i9);
      v10.push_back(i10);
      v11.push_back(i11);
      v12.push_back(i12);
      v13.push_back(i13);
    }

    infile.close();
    v1.pop_back();
    v2.pop_back();
    v3.pop_back();
    v4.pop_back();
    v5.pop_back();
    v6.pop_back();
    v7.pop_back();
    v8.pop_back();
    v9.pop_back();
    v10.pop_back();
    v11.pop_back();
    v12.pop_back();
    v13.pop_back();
  }
}


/*
void setBranchAddresses_mc(TTree *tree, Variables &vars, TString spectrometer_option) {
if (spectrometer_option=="HMS" || spectrometer_option=="test"){
  tree->SetBranchAddress("hsdelta", &vars.hsdp);
  tree->SetBranchAddress("hsxpfp",  &vars.hspxp);
  tree->SetBranchAddress("hsypfp",  &vars.hspyp);
  tree->SetBranchAddress("stop_id", &vars.stop_id);
  tree->SetBranchAddress("hsxfp",   &vars.hpsxfp);
  tree->SetBranchAddress("hsyfp",   &vars.hpsyfp);
  tree->SetBranchAddress("hsxptar", &vars.hpsxptar);
  tree->SetBranchAddress("hsyptar", &vars.hpsyptar);
  tree->SetBranchAddress("hsytar",  &vars.hpsytar);
  tree->SetBranchAddress("hsyptari", &vars.hpsyptari);
  tree->SetBranchAddress("hsxptari", &vars.hpsxptari);
  tree->SetBranchAddress("hsdeltai", &vars.hsdpi);


}
  else if (spectrometer_option=="SHMS"){
  cout<< "TAKING SHMS BRANCHES MC"<<endl;
  tree->SetBranchAddress("psdelta", &vars.hsdp);
  tree->SetBranchAddress("psxpfp",  &vars.hspxp);
  tree->SetBranchAddress("psypfp",  &vars.hspyp);
  tree->SetBranchAddress("stop_id", &vars.stop_id);
  tree->SetBranchAddress("psxfp",   &vars.hpsxfp);
  tree->SetBranchAddress("psyfp",   &vars.hpsyfp);
  tree->SetBranchAddress("psxptar", &vars.hpsxptar);
  tree->SetBranchAddress("psyptar", &vars.hpsyptar);
  tree->SetBranchAddress("psytar",  &vars.hpsytar);
}
  
}

*/

void setBranchAddresses_mc(TTree *tree, Variables &vars, TString spectrometer_option) {
  // Turn off all branches
  tree->SetBranchStatus("*", 0);

  if (spectrometer_option=="HMS" || spectrometer_option=="test") {
    // Enable only the branches you need for HMS
    tree->SetBranchStatus("hsdelta", 1);
    tree->SetBranchStatus("hsxpfp", 1);
    tree->SetBranchStatus("hsypfp", 1);
    tree->SetBranchStatus("stop_id", 1);
    tree->SetBranchStatus("hsxfp", 1);
    tree->SetBranchStatus("hsyfp", 1);
    tree->SetBranchStatus("hsxptar", 1);
    tree->SetBranchStatus("hsyptar", 1);
    tree->SetBranchStatus("hsytar", 1);
    tree->SetBranchStatus("hsyptari", 1);
    tree->SetBranchStatus("hsxptari", 1);
    tree->SetBranchStatus("hsdeltai", 1);

    tree->SetBranchAddress("hsdelta", &vars.hsdp);
    tree->SetBranchAddress("hsxpfp",  &vars.hspxp);
    tree->SetBranchAddress("hsypfp",  &vars.hspyp);
    tree->SetBranchAddress("stop_id", &vars.stop_id);
    tree->SetBranchAddress("hsxfp",   &vars.hpsxfp);
    tree->SetBranchAddress("hsyfp",   &vars.hpsyfp);
    tree->SetBranchAddress("hsxptar", &vars.hpsxptar);
    tree->SetBranchAddress("hsyptar", &vars.hpsyptar);
    tree->SetBranchAddress("hsytar",  &vars.hpsytar);
    tree->SetBranchAddress("hsyptari", &vars.hpsyptari);
    tree->SetBranchAddress("hsxptari", &vars.hpsxptari);
    tree->SetBranchAddress("hsdeltai", &vars.hsdpi);
  }
  else if (spectrometer_option=="SHMS") {
    cout<< "TAKING SHMS BRANCHES MC"<<endl;
    // Enable only the branches you need for SHMS
    tree->SetBranchStatus("psdelta", 1);
    tree->SetBranchStatus("psxpfp", 1);
    tree->SetBranchStatus("psypfp", 1);
    tree->SetBranchStatus("stop_id", 1);
    tree->SetBranchStatus("psxfp", 1);
    tree->SetBranchStatus("psyfp", 1);
    tree->SetBranchStatus("psxptar", 1);
    tree->SetBranchStatus("psyptar", 1);
    tree->SetBranchStatus("psytar", 1);

    tree->SetBranchAddress("psdelta", &vars.hsdp);
    tree->SetBranchAddress("psxpfp",  &vars.hspxp);
    tree->SetBranchAddress("psypfp",  &vars.hspyp);
    tree->SetBranchAddress("stop_id", &vars.stop_id);
    tree->SetBranchAddress("psxfp",   &vars.hpsxfp);
    tree->SetBranchAddress("psyfp",   &vars.hpsyfp);
    tree->SetBranchAddress("psxptar", &vars.hpsxptar);
    tree->SetBranchAddress("psyptar", &vars.hpsyptar);
    tree->SetBranchAddress("psytar",  &vars.hpsytar);
  }
}



void setBranchAddresses_data(TTree *tree, Variables &vars, TString spectrometer_option) {
  // Turn off all branches
  tree->SetBranchStatus("*", 0);

  if (spectrometer_option == "HMS" || spectrometer_option == "test") {
    // Enable only the branches you need for HMS
    tree->SetBranchStatus("H.gtr.dp", 1);
    tree->SetBranchStatus("H.gtr.th", 1);
    tree->SetBranchStatus("H.gtr.ph", 1);
    tree->SetBranchStatus("H.gtr.x", 1);
    tree->SetBranchStatus("H.gtr.y", 1);
    tree->SetBranchStatus("H.dc.x_fp", 1);
    tree->SetBranchStatus("H.dc.y_fp", 1);
    tree->SetBranchStatus("H.dc.xp_fp", 1);
    tree->SetBranchStatus("H.dc.yp_fp", 1);
    tree->SetBranchStatus("H.kin.W", 1);
    tree->SetBranchStatus("H.kin.Q2", 1);
    tree->SetBranchStatus("H.kin.nu", 1);
    tree->SetBranchStatus("H.kin.scat_ang_deg", 1);
    tree->SetBranchStatus("H.cal.etottracknorm", 1);
    tree->SetBranchStatus("H.cer.npeSum", 1);
    tree->SetBranchStatus("H.bcm.bcm4a.AvgCurrent", 1);
    tree->SetBranchStatus("H.bcm.CurrentFlag", 1);

    tree->SetBranchAddress("H.gtr.dp",            &vars.ptardp); 
    tree->SetBranchAddress("H.gtr.th",            &vars.ptarth); 
    tree->SetBranchAddress("H.gtr.ph",            &vars.ptarph); 
    tree->SetBranchAddress("H.gtr.x",             &vars.ptarx);  
    tree->SetBranchAddress("H.gtr.y",             &vars.ptary);  
    tree->SetBranchAddress("H.dc.x_fp",           &vars.dpsxfp); 
    tree->SetBranchAddress("H.dc.y_fp",           &vars.dpsyfp); 
    tree->SetBranchAddress("H.dc.xp_fp",          &vars.dpsxptar);
    tree->SetBranchAddress("H.dc.yp_fp",          &vars.dpsyptar);
    tree->SetBranchAddress("H.kin.W",             &vars.dW);     
    tree->SetBranchAddress("H.kin.Q2",            &vars.dQ2);    
    tree->SetBranchAddress("H.kin.nu",            &vars.dnu);    
    tree->SetBranchAddress("H.kin.scat_ang_deg",  &vars.dtheta); 
    tree->SetBranchAddress("H.cal.etottracknorm", &vars.decal);
    tree->SetBranchAddress("H.cer.npeSum",        &vars.npeSum);
    tree->SetBranchAddress("H.bcm.bcm4a.AvgCurrent", &vars.AvgCurr);
    tree->SetBranchAddress("H.bcm.CurrentFlag"     , &vars.CurrentFlag);
  } 
  else if (spectrometer_option == "SHMS") {
    cout<< "TAKING SHMS BRANCHES"<<endl;
    // Enable only the branches you need for SHMS
    tree->SetBranchStatus("P.gtr.dp", 1);
    tree->SetBranchStatus("P.gtr.th", 1);
    tree->SetBranchStatus("P.gtr.ph", 1);
    tree->SetBranchStatus("P.gtr.x", 1);
    tree->SetBranchStatus("P.gtr.y", 1);
    tree->SetBranchStatus("P.dc.x_fp", 1);
    tree->SetBranchStatus("P.dc.y_fp", 1);
    tree->SetBranchStatus("P.dc.xp_fp", 1);
    tree->SetBranchStatus("P.dc.yp_fp", 1);
    tree->SetBranchStatus("P.kin.W", 1);
    tree->SetBranchStatus("P.kin.Q2", 1);
    tree->SetBranchStatus("P.kin.nu", 1);
    tree->SetBranchStatus("P.kin.scat_ang_deg", 1);
    tree->SetBranchStatus("P.cal.etottracknorm", 1);
    tree->SetBranchStatus("P.ngcer.npeSum", 1);

    tree->SetBranchAddress("P.gtr.dp",            &vars.ptardp); 
    tree->SetBranchAddress("P.gtr.th",            &vars.ptarth); 
    tree->SetBranchAddress("P.gtr.ph",            &vars.ptarph); 
    tree->SetBranchAddress("P.gtr.x",             &vars.ptarx);  
    tree->SetBranchAddress("P.gtr.y",             &vars.ptary);  
    tree->SetBranchAddress("P.dc.x_fp",           &vars.dpsxfp); 
    tree->SetBranchAddress("P.dc.y_fp",           &vars.dpsyfp); 
    tree->SetBranchAddress("P.dc.xp_fp",          &vars.dpsxptar);
    tree->SetBranchAddress("P.dc.yp_fp",          &vars.dpsyptar);
    tree->SetBranchAddress("P.kin.W",             &vars.dW);     
    tree->SetBranchAddress("P.kin.Q2",            &vars.dQ2);    
    tree->SetBranchAddress("P.kin.nu",            &vars.dnu);    
    tree->SetBranchAddress("P.kin.scat_ang_deg",  &vars.dtheta); 
    tree->SetBranchAddress("P.cal.etottracknorm", &vars.decal);
    tree->SetBranchAddress("P.ngcer.npeSum",      &vars.npeSum);
  }
}


double calculate_mc_scale_factor(int nentries_mc , double p_spec, double density , double length, double AM , double sim_charge , int delta_hi, int delta_lo , TString targType ) {
  cout<< "NOTE THE p_spec : " << p_spec <<endl;
  double domega           = (phi_hi - phi_lo)*(theta_hi-theta_lo) / 1000. /1000.; // differential solid angle in sr
  double ep_max           = p_spec*(1.+0.01*delta_hi);
  double ep_min           = p_spec*(1.+0.01*delta_lo);
  double deltaEprime      = ep_max-ep_min;
  double target_thickness = density*length;
  if (targType=="LT")     { target_thickness = (target_thickness -0.0168-0.02024)*0.996 ; } // cryotarget correction factor
  double luminosity       = target_thickness*sim_charge/AM*N_A/Q_E*1e-39;
  double factor           = (luminosity * domega * deltaEprime) / nentries_mc ;
  return factor;
}

double calculate_BjorkenX ( double p_spec , double hsdp , double hpsyptar , double hpsxptar , double& Eprime , double& thetadeg, double theta_central, TString spectrometer_option ) {
  Eprime   = p_spec*(1+0.01*hsdp);
  double thetarad; 
 if      ( spectrometer_option=="HMS" || spectrometer_option =="test")  thetarad = TMath::ACos((cos(theta_central) + hpsyptar*sin(theta_central))/TMath::Sqrt(1. + hpsxptar*hpsxptar+hpsyptar*hpsyptar));
 else if ( spectrometer_option=="SHMS") thetarad = TMath::ACos((cos(theta_central) - hpsyptar*sin(theta_central))/TMath::Sqrt(1. + hpsxptar*hpsxptar+hpsyptar*hpsyptar)); 
 
  thetadeg = thetarad*(180.0/3.1415);
  double Q2       = 4.0 * Einitial * Eprime * (TMath::Sin(thetarad / 2.0) * TMath::Sin(thetarad / 2.0));
  double nu       = Einitial - Eprime;
  double W2       = -(Q2) + (Mp * Mp) + ( 2.0 * Mp * nu);
  double x        = (Q2)/2./Mp/(nu);
  return x;
}


TGraph2DErrors* createGraph2D( const std::vector<double>& V2, const std::vector<double>& V3, const std::vector<double>& VData, int size) {
  TGraph2DErrors *graph = new TGraph2DErrors();
  for (Int_t j = 0; j < size; j++) {
    graph->SetPoint(j, V2[j], V3[j], VData[j]);
    graph->SetPointError(j, 0., 0., 0.);
  }
  return graph;
}


void printHistogramContents(TH1F *hist) {

  int numBins = hist->GetNbinsX();
  std::cout << "x_values = [";
  for (int i = 1; i <= numBins; ++i) {
    double x = hist->GetXaxis()->GetBinCenter(i);
    if (i == numBins) {
      std::cout << x;
    } else {
      std::cout << x << ", ";
    }
  }
  std::cout << "]" << std::endl;

  std::cout << "bin_contents = [";
  for (int i = 1; i <= numBins; ++i) {
    double content = hist->GetBinContent(i);
    if (i == numBins) {
      std::cout << content;
    } else {
      std::cout << content << ", ";
    }
  }
  std::cout << "]" << std::endl;

  std::cout << "bin_errors = [";
  for (int i = 1; i <= numBins; ++i) {
    double error = hist->GetBinError(i);
    if (i == numBins) {
      std::cout << error;
    } else {
      std::cout << error << ", ";
    }
  }
  std::cout << "]" << std::endl;
}




void PlotWithRatio(TCanvas *&canvas, TH1F *Data_ST, TH1F *Data_LT, TH1F *MC_ST, TH1F *MC_LT, TH1F *MC_ST_UNWeigthed,TH1F *MC_LT_UNWeigthed , TString out_name) {
  canvas = new TCanvas("cc1", "Cross-Sections", 800, 600);
  canvas->Divide(4,2);
  for (int i = 1; i <= 8; ++i) {
    canvas->cd(i);
    gPad->SetTicks(1,1);
    gPad->SetGrid();
    gStyle->SetGridColor(kGray);
    //gPad->SetLogy();                                                                      
  }

  canvas->cd(2);
  Data_ST->SetTitle("ST");
  Data_ST->SetMarkerSize(2);
  Data_ST->SetStats(0);
  Data_ST->Draw("E");
  MC_ST->SetLineColor(kRed);
  MC_ST->SetMarkerSize(3);
  MC_ST->Draw("same");

  canvas->cd(1);

  MC_ST_UNWeigthed->SetLineColor(kGreen);
  MC_ST_UNWeigthed->SetMarkerSize(3);
  MC_ST_UNWeigthed->SetTitle("MC ST Unweighted");
  MC_ST_UNWeigthed->SetStats(0);
  MC_ST_UNWeigthed->Draw("E");


  canvas->cd(4);
  Data_LT->SetTitle("LD2");
  Data_LT->SetMarkerSize(3);
  Data_LT->SetStats(0);
  Data_LT->Draw("E");
  MC_LT->SetLineColor(kRed);
  MC_LT->SetMarkerSize(3);
  MC_LT->Draw("same");

  canvas->cd(3);
  MC_LT_UNWeigthed->SetLineColor(kGreen);
  MC_LT_UNWeigthed->SetMarkerSize(3);
  MC_LT_UNWeigthed->SetStats(0);
  MC_LT_UNWeigthed->SetTitle("MC LD2 Unweighted");
  MC_LT_UNWeigthed->Draw("E");

  // Ratio plot                                                                                 
  canvas->cd(5);
  TH1F *ratio_ST = (TH1F*)Data_ST->Clone("ratio_ST");
  ratio_ST->Divide(MC_ST);
  ratio_ST->SetTitle("Data ST / MC ST");
  ratio_ST->SetMarkerSize(2);
  ratio_ST->GetYaxis()->SetTitle("Data / MC");
  ratio_ST->GetYaxis()->SetNdivisions(505);
  ratio_ST->GetYaxis()->SetTitleSize(0.06);
  ratio_ST->GetYaxis()->SetTitleOffset(0.55);
  ratio_ST->SetStats(0);
  ratio_ST->Draw("E");

  canvas->cd(6);
  TH1F *ratio_LT = (TH1F*)Data_LT->Clone("ratio_LT");
  ratio_LT->Divide(MC_LT);
  ratio_LT->SetTitle("Data LD2 / MC LD2");
  ratio_LT->SetMarkerSize(2);
  ratio_LT->GetYaxis()->SetTitle("Data / MC");
  ratio_LT->GetYaxis()->SetNdivisions(505);
  ratio_LT->GetYaxis()->SetTitleSize(0.06);
  ratio_LT->GetYaxis()->SetTitleOffset(0.55);
  ratio_LT->SetStats(0);
  ratio_LT->Draw("E");

  TH1F *ratio_CrossSections = (TH1F*)ratio_ST->Clone("ratio_crossSections");
  ratio_CrossSections->Divide(ratio_LT);
  canvas->cd(7);
  ratio_CrossSections->SetTitle("Cross-Section Ratio (ST/LD2)");
  ratio_CrossSections->SetMarkerSize(2);
  ratio_CrossSections->SetStats(0);
  ratio_CrossSections->Draw("E");

  // Legend                                                                     
  TLegend *legend = new TLegend(0.1, 0.1, 0.4, 0.25); // Adjust legend position                                                                   
  legend->AddEntry(Data_ST, "Data");
  legend->AddEntry(MC_ST, "MC");
  legend->AddEntry(MC_ST_UNWeigthed, "MC unwt");
  legend->SetTextSize(0.04);
  canvas->cd(1);
  legend->Draw();

  canvas->SaveAs(out_name + "_Intermediate_Steps.pdf");
}

double deltaToX(double delta , double p_spec , double Einitial , double angle) {
  double Eprime            = p_spec*(1+0.01*delta);
  double angle_radians     = angle * 3.1418/180.0;
  double Q2                = 4.0 * Einitial * Eprime * (TMath::Sin( angle_radians  / 2.0) * TMath::Sin(angle_radians / 2.0));
  double nu                = Einitial - Eprime;
  double X_Bjorken         = (Q2)/2./Mp/(nu);
  return X_Bjorken;
}

double XTodelta(double X_Bjorken, double p_spec, double Einitial, double angle) {
  double angle_radians = angle * 3.1418 / 180.0;
  double sin2 = TMath::Sin(angle_radians / 2.0) * TMath::Sin(angle_radians / 2.0);
  double Eprime = (2 * Mp * Einitial * X_Bjorken) / (X_Bjorken * 2 * Mp + 4 * Einitial * sin2);
  double delta = (Eprime / p_spec - 1) * 100.0;

  return delta;
}
void convertHistogramToGraph(TH1F* h_delta, TGraphErrors*& graph,  double p_spec , double Einitial , double angle , TString out_name) {
  int nBins = h_delta->GetNbinsX();

  std::vector<double> xValues;
  std::vector<double> yValues;
  std::vector<double> exValues;
  std::vector<double> eyValues;

  for (int i = 1; i <= nBins; ++i) {
    double delta   = h_delta->GetXaxis()->GetBinCenter(i);
    double x       = deltaToX(delta, p_spec, Einitial, angle);
    //    cout<< "FOR BIN "<< i <<" DELTA IS: " << delta << ", AND X IS: "<< x<<endl;                                                                                                                                                                                           
    double content = h_delta->GetBinContent(i);
    double error   = h_delta->GetBinError(i);

    xValues.push_back(x);
    yValues.push_back(content);
    exValues.push_back(0.0); // Assume no error in x for simplicity                                                                                                                                                                                                             
    eyValues.push_back(error);
  }


  graph = new TGraphErrors(nBins, xValues.data(), yValues.data(), exValues.data(), eyValues.data());
  graph->SetTitle("CrossSections as a function of Xb;x;Ratio");


  TCanvas* canvas = new TCanvas("canvas", "Histograms and Graphs", 800, 600);
  canvas->Divide(2, 1);
  canvas->cd(1);
  h_delta->Draw();
  canvas->cd(2);
  graph->Draw("AP"); // "AP" means draw with axis and points                                                                              

  canvas->SaveAs(out_name + "_TGraphErrors_Comparison_delta_X.pdf");

}

void printGraphContents(TGraphErrors *graph) {
  int numPoints = graph->GetN();

  std::cout << "x_values = [";
  for (int i = 0; i < numPoints; ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    if (i == numPoints - 1) {
      std::cout << x;
    } else {
      std::cout << x << ", ";
    }
  }
  std::cout << "]" << std::endl;

  std::cout << "y_values = [";
  for (int i = 0; i < numPoints; ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    if (i == numPoints - 1) {
      std::cout << y;
    } else {
      std::cout << y << ", ";
    }
  }
  std::cout << "]" << std::endl;

  std::cout << "x_errors = [";
  for (int i = 0; i < numPoints; ++i) {
    double x, y;
    double ex, ey;
    graph->GetPoint(i, x, y);
    ex = graph->GetErrorX(i);
    if (i == numPoints - 1) {
      std::cout << ex;
    } else {
      std::cout << ex << ", ";
    }
  }
  std::cout << "]" << std::endl;

  std::cout << "y_errors = [";
  for (int i = 0; i < numPoints; ++i) {
    double x, y;
    double ex, ey;
    graph->GetPoint(i, x, y);
    ey = graph->GetErrorY(i);
    if (i == numPoints - 1) {
      std::cout << ey;
    } else {
      std::cout << ey << ", ";
    }
  }
  std::cout << "]" << std::endl;
}

void writeHistogramToTextFile(TH1F* hist, double p_spec , double Einitial , double angle, const char* filename) {
  std::ofstream outFile(filename);
  if (!outFile.is_open()) { std::cerr << "Error opening file: " << filename << std::endl;return;  }
  outFile << "Bin\tXb\tdelta\tQ2\tEp\txi\tratio\tError\n";  
  int nBins = hist->GetNbinsX();
  for (int i = 1; i <= nBins; ++i) {
    int binNumber     = i;
    double delta , xb; 
    if (var_dep=="delta") { delta      = hist->GetXaxis()->GetBinCenter(i);  }
    if (var_dep=="xb")    { xb         = hist->GetXaxis()->GetBinCenter(i);  }

    double binContent = hist->GetBinContent(i);
    double binError   = hist->GetBinError(i);

    if (var_dep=="delta") {  xb         = deltaToX(delta, p_spec, Einitial, angle);} 
    if (var_dep=="xb")    {  delta      = XTodelta( xb ,  p_spec, Einitial, angle);}

    double Eprime     = p_spec*(1+0.01*delta);
    double Q2         = 2.0*Mp*xb*(Einitial-Eprime);
    double nu         = Einitial-Eprime;
    double xi         = (2*xb) / ( 1+TMath::Sqrt(1+( Q2/(nu*nu))));

    outFile << binNumber << "\t" << xb <<"\t"<< delta << "\t" << Q2<< "\t"<< Eprime<< "\t"<< xi <<"\t" << binContent << "\t" << binError << "\n";
  }



  outFile.close();
  std::cout << "Histogram data written to " << filename << std::endl;
}

void scaleGraph(TGraph* graph, double scaleFactor) {
  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    y /= scaleFactor;  // Divide by scaleFactor (in this case, 6)                                                                                                                                                                                                                                                            
    graph->SetPoint(i, x, y);
  }
}

void plottingWeightFactors( std::vector<double> hsdp_values_ST , std::vector<double> weight_values_ST, std::vector<double> hsdp_values_LT , std::vector<double> weight_values_LT, TString out_name){
  int n_ST = hsdp_values_ST.size();
  TGraph *graph_ST = new TGraph(n_ST, &hsdp_values_ST[0], &weight_values_ST[0]);
  graph_ST->SetTitle("MC Weight vs delta;delta;MC Weight");
  graph_ST->SetMarkerStyle(20);
  graph_ST->SetMarkerColor(kBlue);

  int n_LT = hsdp_values_LT.size();
  TGraph *graph_LT = new TGraph(n_LT, &hsdp_values_LT[0], &weight_values_LT[0]);
  graph_LT->SetTitle("MC Weight vs delta;delta;MC Weight");
  graph_LT->SetMarkerStyle(20);
  graph_LT->SetMarkerColor(kRed);

  double y_max = std::max(*std::max_element(weight_values_ST.begin(), weight_values_ST.end()),
                          *std::max_element(weight_values_LT.begin(), weight_values_LT.end()));

  graph_ST->GetXaxis()->SetLimits(-8.5, 8.5);
  graph_ST->GetYaxis()->SetRangeUser(0.0, y_max*1.1);
  scaleGraph(graph_ST, 6.0); // for carbon is 6.0, make it more general later!!
  TCanvas *c2 = new TCanvas("c2", "Weight vs hsdp", 800, 600);
  graph_ST->Draw("AP");
  graph_LT->Draw("Psame");

  TLegend *legend = new TLegend(0.8, 0.8, 0.95, 0.95);
  legend->AddEntry(graph_ST, "ST", "p");
  legend->AddEntry(graph_LT, "LD2", "p");
  legend->Draw();

  c2->Update();
  c2->SaveAs(out_name + "_weights_MC.png");

}

void plottingWeightFactors_Born( std::vector<double> hsdp_values_ST , std::vector<double> weight_values_ST, std::vector<double> hsdp_values_LT , std::vector<double> weight_values_LT, TString out_name){
  int n_ST = hsdp_values_ST.size();
  TGraph *graph_ST = new TGraph(n_ST, &hsdp_values_ST[0], &weight_values_ST[0]);
  graph_ST->SetTitle(" Sigma Born vs delta;delta;Born Cross-section");
  graph_ST->SetMarkerStyle(20);
  graph_ST->SetMarkerColor(kBlue);

  int n_LT = hsdp_values_LT.size();
  TGraph *graph_LT = new TGraph(n_LT, &hsdp_values_LT[0], &weight_values_LT[0]);
  graph_LT->SetTitle(" Sigma Born vs delta;delta;Born Cross-section");
  graph_LT->SetMarkerStyle(20);
  graph_LT->SetMarkerColor(kRed);

  double y_max = std::max(*std::max_element(weight_values_ST.begin(), weight_values_ST.end()),
                          *std::max_element(weight_values_LT.begin(), weight_values_LT.end()));

  graph_ST->GetXaxis()->SetLimits(-8.5, 8.5);
  graph_ST->GetYaxis()->SetRangeUser(0.0, y_max*1.1);

  scaleGraph(graph_ST, 6.0); // for carbon is 6.0, make it more general later!! 
  TCanvas *c_tmp = new TCanvas("c22", "Weight vs delta", 800, 600);
  graph_ST->Draw("AP");
  graph_LT->Draw("Psame");

  TLegend *legend = new TLegend(0.8, 0.8, 0.95, 0.95);
  legend->AddEntry(graph_ST, "ST", "p");
  legend->AddEntry(graph_LT, "LD2", "p");
  legend->Draw();

  c_tmp->Update();
  c_tmp->SaveAs(out_name + "_weights_Born.png");

}

void mergePDFs( TString out_name ) {

  const char* pdfFiles[] = {
    out_name + "_Intermediate_Steps.pdf",
    out_name + "_DELTA.pdf"
  };
  int nFiles = sizeof(pdfFiles) / sizeof(pdfFiles[0]);

  const char* mergedOutput = out_name + "merged_output.pdf";
  std::string command = "pdftk ";
  for (int i = 0; i < nFiles; ++i) {
    command += pdfFiles[i];
    command += " ";
  }
  command += "cat output ";
  command += mergedOutput;
  int status = gSystem->Exec(command.c_str());
  if (status == 0) {
    std::cout << "Successfully merged PDFs into " << mergedOutput << std::endl;
  } else {
    std::cerr << "Failed to merge PDFs with pdftk" << std::endl;
  }

}

double linearInterpolation(double x0, double y0, double x1, double y1, double x) {
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}


bool areEqual(double a, double b, double tolerance) {
  return std::fabs(a - b) < tolerance;
}

double Interpolate1DFromGraph2DErrors(TGraph2DErrors* graph, double Eprime, double angle) {

  std::vector<double> x_values, y_values, z_values;
  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y, z;
    graph->GetPoint(i, x, y, z);                                                                           
    if (areEqual(y, angle, 1e-8)) { // Only consider points where y is equal to the angle  
      x_values.push_back(x);
      y_values.push_back(y);
      z_values.push_back(z);
    }
  }

  if (x_values.size() < 2) {
    cout<< "the angle is : "<<angle <<endl;
    std::cout << "Error: Not enough points for interpolation" << std::endl;
    return NAN;
  }

  std::vector<std::tuple<double, double>> points;
  for (size_t i = 0; i < x_values.size(); ++i) {
    points.emplace_back(x_values[i], z_values[i]);
  }
  std::sort(points.begin(), points.end());

  size_t index = 0;
  while (index < points.size() && std::get<0>(points[index]) < Eprime) {
    index++;
  }

  if (index == 0 || index == points.size()) {
    std::cout << "Error: Interpolation point is outside the range of the data" << std::endl;
    return NAN;
  }

  double x0 = std::get<0>(points[index - 1]);
  double x1 = std::get<0>(points[index]);
  double y0 = std::get<1>(points[index - 1]);
  double y1 = std::get<1>(points[index]);

  double interpolated_value = linearInterpolation(x0, y0, x1, y1, Eprime);

  return interpolated_value;
}






















double findNearestValue(double input, double angle) {
    std::vector<double> validValues;
    double minValue, maxValue;
    const double tolerance = 0.5; // Tolerance for angle matching

    // Check for angles within tolerance
    if (fabs(angle - 20.0) <= tolerance) {
        minValue = 14.0;
        maxValue = 28.0;
    } else if (fabs(angle - 26.0) <= tolerance) {
        minValue = 19.0;
        maxValue = 34.0;
    } else if (fabs(angle - 35.0) <= tolerance) {
        minValue = 29.0;
        maxValue = 43.0;
    } else {
    std::cout << "Angle " << angle << " is not within tolerance of 20, 26, or 35 degrees." << std::endl;
        return -1; // Return an error code or handle as needed
    }

    // Populate the valid values vector
    for (double val = minValue; val <= maxValue; val += 0.2) {
        validValues.push_back(val);
    }

    // Check if the input is out of range
    if (input < minValue) {
    std::cout << "Input value " << input << " is too small. Minimum allowed value is " << minValue << "." << std::endl;
        return minValue;
    }
    if (input > maxValue) {
    std::cout << "Input value " << input << " is too large. Maximum allowed value is " << maxValue << "." << std::endl;
        return maxValue;
    }

    // Find the nearest value in the validValues vector
    double nearest = validValues[0];
    double minDiff = fabs(input - nearest);
    for (double val : validValues) {
        double diff = fabs(input - val);
        if (diff < minDiff) {
            minDiff = diff;
            nearest = val;
        }
    }

    return nearest;
}































/*





double findNearestValue(double input, double angle) {
  std::vector<double> validValues;
  double minValue , maxValue;
  if      (angle==20.0){minValue=14.0; maxValue=28.0 ;     }
  else if (angle==26.0){minValue=19.0; maxValue=34.0 ;     }
  else if (angle==35.0){minValue=29.0; maxValue=43.0 ;     }
  for (double val = minValue; val <= maxValue; val += 0.2) {
    validValues.push_back(val);
  }

  if (input < minValue) {
    std::cout << "Input value " << input << " is too small. Minimum allowed value is " << minValue << "." << std::endl;
    return minValue;
  }
  if (input > maxValue) {
    std::cout << "Input value " << input << " is too large. Maximum allowed value is " << maxValue << "." << std::endl;
    return maxValue;
  }

  double nearestValue = validValues[0];
  double minDifference = std::abs(input - validValues[0]);

  for (const auto& value : validValues) {
    double difference = std::abs(input - value);
    if (difference < minDifference) {
      minDifference = difference;
      nearestValue = value;
    }
  }

  return nearestValue;
}

*/
void plotHistograms_Dummy(TH1F* X_Data_LT, TH1F* X_Data_LT_Dummy,  TString out_name) {
  TCanvas *c = new TCanvas("c", "Histograms", 800, 600);
  c->SetLogy();                                    
  TH1F *X_Data_LT_before = (TH1F*)X_Data_LT->Clone("X_Data_LT_before");
  TH1F *X_Data_LT_after = (TH1F*)X_Data_LT->Clone("X_Data_LT_after");

  X_Data_LT_after->Add(X_Data_LT_Dummy, -1); // Subtracting dummy                                                                                                     
  X_Data_LT_before->SetLineColor(kBlue);
  X_Data_LT_Dummy->SetLineColor(kRed);
  X_Data_LT_after->SetLineColor(kGreen);

  double min_val = 1e-3;
  X_Data_LT_before->SetMinimum(min_val);
  X_Data_LT_before->Draw("hist");
  X_Data_LT_Dummy->Draw("hist same");
  X_Data_LT_after->Draw("hist same");

  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(X_Data_LT_before, "LD2", "l");
  legend->AddEntry(X_Data_LT_Dummy, "dummy", "l");
  legend->AddEntry(X_Data_LT_after, "LD2-dummy", "l");
  legend->Draw();

  // Save the canvas                                                                                      
  c->SaveAs(out_name + "dummy_comparison.png");
}


void plotHistograms_Contamination(TH1F* X_Data_ST, TH1F* X_Data_ST_tmp,  TString out_name) {
  TCanvas *c = new TCanvas("c", "Histograms", 800, 600);
  c->SetLogy();
  TH1F *X_Data_ST_before = (TH1F*)X_Data_ST->Clone("X_Data_ST_before");
  TH1F *X_Data_ST_after = (TH1F*)X_Data_ST->Clone("X_Data_ST_after");

  X_Data_ST_after->Add(X_Data_ST_tmp, -1); // Subtracting contamination
  X_Data_ST_before->SetLineColor(kBlue);
  X_Data_ST_tmp->SetLineColor(kRed);
  X_Data_ST_after->SetLineColor(kGreen);

  double min_val = 1e-3;
  X_Data_ST_before->SetMinimum(min_val);
  X_Data_ST_before->Draw("hist");
  X_Data_ST_tmp->Draw("hist same");
  X_Data_ST_after->Draw("hist same");

  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(X_Data_ST_before, "LD2", "l");
  legend->AddEntry(X_Data_ST_tmp, "contamination", "l");
  legend->AddEntry(X_Data_ST_after, "LD2-contamination", "l");
  legend->Draw();

  // Save the canvas  
  c->SaveAs(out_name + "contamination_comparison.png");
}

double weight_calculateF2(double theta_deg , double Eprime , double Eb){

double nu = Eb - Eprime;
double angle_radians     = theta_deg * 3.1418/180.0;
double Q2                = 4.0 * Einitial * Eprime * (TMath::Sin( angle_radians  / 2.0) * TMath::Sin(angle_radians / 2.0));
double R = 0.32/Q2 ; 
double factor1 = nu /(1+ 2*TMath::Tan(angle_radians  / 2.0)*TMath::Tan(angle_radians  / 2.0) * (  (1+ (nu*nu)/(Q2)   )/(1+R) )   );
double alpha = 1/137.0;
double mottxs = ( alpha*alpha * TMath::Cos(angle_radians  / 2.0) * TMath::Cos(angle_radians  / 2.0)    ) / (4*Eb*Eb* TMath::Power(TMath::Sin(angle_radians / 2.0), 4) ) ; 

double weight  = factor1/mottxs;
return weight;

}




void fillScaleFactors(const TString& table_file,
                      const TString& target_input,
                      const TString& angle_input,
                      const TString& spectrometer,
                      const TString& momentum_spec_input,
                      double p_spec,
                      std::vector<ScaleFactors>& run_data_vector) {

  std::ifstream file(table_file.Data());

  //  cout<< "checking file : " << table_file.Data() <<endl; // run_numbers_info.txt

  std::string line;

  if (file.is_open()) {
    // Clear existing data                                                                                                                                               
    run_data_vector.clear();

    while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#') {
        continue; // Skip empty lines or comments                                                                                                                        
      }

      std::stringstream ss(line);
      std::string target, angle, spectr, momentum;
      double p_spec_in_file;
      ScaleFactors sf ;
      //      cout << "p_spec_in_file " << p_spec_in_file << " ,   p_spec: " << p_spec<<endl;
      ss >> target >> angle >> spectr >> momentum >> p_spec_in_file
         >> sf.run_number >> sf.charge >> sf.PS >> sf.Eff_Fid
         >> sf.Eff_Elec_LiveTime >> sf.Eff_Comp_LiveTime >> sf.trigger;

      if (TString(target) == target_input && TString(angle) == angle_input &&
          TString(spectr) == spectrometer && TString(momentum) == momentum_spec_input &&
          fabs(p_spec_in_file -p_spec) < 0.016) {
	run_data_vector.push_back(sf);
      }
    }

    file.close();
  } else {
    std::cerr << "Unable to open file: " << table_file << std::endl;
  }
}










