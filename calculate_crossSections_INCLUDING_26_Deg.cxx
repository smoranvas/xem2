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
const float  N_A            = 6.02*1e+23;         // number of particles per one gram
const float  Q_E            = 1.60*1e-19;         //1C = 6.24e18 electrons
const int    phi_hi         =  100;
const int    phi_lo         = -100;
const int    theta_hi       =  100;
const int    theta_lo       = -100;
const double Einitial       = 10.544;
const double Mp             = 0.938;
const float  density_LT     = 0.16743;       //Has g/cm3 units 
const double length_LT      = 10.0;          //Has cm units  
const double sim_charge_ST  = 1.;            //Has uC units 
const double sim_charge_LT  = 1.;            //Has uC units 
const double AM_LT          = 2.0;
const int    N_BINS         = 80;
const double TOLERANCE      = 1e-8;
const bool   F2_option      = false;     //it should be false always, calculating F2 event by event is not the right approach
const char* var_dep         = "delta";      // delta or xb
const char* trigger_select  = "ELREAL";  // ELREAL or anything else to include 3/4 triggers (not reason to include those in physics analysis)

struct Variables {
  // Variables for SIMULATION
  float hsdp;
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
  float hsdpi;
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
};

struct ScaleFactors {
  int    run_number;
  double charge;
  int    PS;
  double Eff_Fid ;
  double Eff_Elec_LiveTime;
  double Eff_Comp_LiveTime;
  double trigger ; 
};
void fillScaleFactors(const std::string& table_file, 
                      const std::string& target_input, 
                      const std::string& angle_input, 
                      const std::string& spectrometer, 
                      const std::string& momentum_spec_input, 
                      double p_spec,
                      std::vector<ScaleFactors>& run_data_LT, 
                      std::vector<ScaleFactors>& run_data_ST, 
                      std::vector<ScaleFactors>& run_data_LT_Dummy) ;

void ImportRadcor(vector<double> &v1,vector<double> &v2,vector<double> &v3,vector<double> &v4,vector<double> &v5,vector<double> &v6,vector<double> &v7,vector<double> &v8,vector<double> &v9, vector<double> &v10,vector<double> &v11,vector<double> &v12,vector<double> &v13, const char * filename);
void setBranchAddresses_mc(TTree *tree, Variables &vars, TString spectrometer_option); 
void setBranchAddresses_data(TTree *tree, Variables &vars, TString spectrometer_option) ;
double calculate_mc_scale_factor(int nentries_mc , double p_spec, double density , double length, double AM , double sim_charge , int delta_hi, int delta_lo ) ;
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
double final_weight = 1.0/weight;
return final_weight;

}


double weight_delta_corr(double delta ,  TString spectrometer_option){

double delta_corr;
if ( spectrometer_option=="HMS") { delta_corr =  0.990337*delta - 0.00236077*TMath::Power(delta,2) + 0.000286814*TMath::Power(delta,3) + 2.09878E-6*TMath::Power(delta,4) - 2.4867E-6*TMath::Power(delta,5) + 1.8646E-7*TMath::Power(delta,6) ;   }
else if ( spectrometer_option=="SHMS") { delta_corr = delta;}

return delta_corr;

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

  // Print to console
  std::cout << "Total Entries: " << totalEntries << std::endl;
  std::cout << "Total Error: " << totalError << std::endl;

  // Return total entries and total error as a pair
  return std::make_pair(totalEntries, totalError);
}

void plotHistograms_Dummy(TH1F* X_Data_LT, TH1F* X_Data_LT_Dummy,  TString out_name);

std::vector<ScaleFactors> run_data_LT , run_data_ST , run_data_LT_Dummy;   

void cross_Sections_updated_workingProcess (const char* target_in = "C12", const char* angle_in = "20", const char* momentum_spec_in = "2p42", TString spectrometer_option="HMS" , TString ytarg_corr_opt ="ON", TString jacob_corr_opt ="ON" , TString delta_corr_opt="ON" ) {

  TString target_input        = target_in;
  TString angle_input         = angle_in;
  TString momentum_spec_input = momentum_spec_in;
  TString spectrometer        = spectrometer_option;

  bool ytar_corr_ON , MC_Jacobian_corr_ON , delta_corr_ON , Coulomb_corr_ON ; 
  if (ytarg_corr_opt =="ON"){ytar_corr_ON        = true ; }       else {ytar_corr_ON        = false;}
  if (jacob_corr_opt =="ON"){MC_Jacobian_corr_ON = true ; }       else {MC_Jacobian_corr_ON = false;}
  if (delta_corr_opt =="ON"){delta_corr_ON       = true ; }       else {delta_corr_ON       = false;}

  TString spectrometer_lowercase;
  if ( spectrometer=="HMS" || spectrometer=="test") {spectrometer_lowercase = "hms";}
  else if ( spectrometer=="SHMS") {spectrometer_lowercase = "shms";}

  TString f_st_mc_st , f_lt_mc_st  ,mc_file_lt_st , mc_file_st_st  ; 
  double angle , p_spec , AM_ST;
  int delta_hi , delta_lo;
  double density_ST, length_ST  ;

  if  ( target_input=="C12")  {AM_ST         = 12.0;}
  if  ( angle_input=="20"  )  {angle         = 20.0;}
  else if (angle_input=="26") {angle         = 26.0;}
  else if (angle_input=="35") {angle         = 35.0;}

  if  ( target_input=="C12" && angle_input=="20" && spectrometer=="test"){

    mc_file_lt_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_d2cryo22_hms.out";
    mc_file_st_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_12carbon22_hms.out";
    density_ST    = 2.00;
    length_ST     = 0.287;
    delta_hi      = 20;
    delta_lo      = -20;
    f_st_mc_st    = "/work/smoran/xem2/sim/test_5p87_carbon.root";
    f_lt_mc_st    = "/work/smoran/xem2/sim/test_5p87_deuterium.root";

    p_spec  = 5.878;

    run_data_ST = {
      {5828 , 142189 , 1 , 0.9969 , 1 , 0.999881 , 2}
    };
    run_data_LT = {
      {4780 , 90719.6 , 1 , 0.9972 , 1.00004 , 0.999827 , 2}
    };
    run_data_LT_Dummy = {
      {4772 , 79655.2 , 1 , 0.9976 , 1.00004 , 0.999886 , 2}
    };

}

  if (  target_input=="C12" && angle_input=="35" && spectrometer=="HMS" ){

    mc_file_lt_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_d2cryo22_hms.out";
    mc_file_st_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_12carbon22_hms.out";
    density_ST    = 2.00;
    length_ST     = 0.287;
    delta_hi      = 20;
    delta_lo      = -20;
    f_st_mc_st    = Form("/work/smoran/xem2/sim/hms_35deg_carbon_%s.root"   , (const char*)momentum_spec_input) ;    
    f_lt_mc_st    = Form("/work/smoran/xem2/sim/hms_35deg_deuterium_%s.root", (const char*)momentum_spec_input) ;


if      ( momentum_spec_input=="4p08"   ){
	p_spec  = 4.08; 
	run_data_ST = {

	  {5100 , 80743 , 1 , 1 , 1 , 1 , 2},
	  {5134 , 129744 , 1 , 1 , 1 , 0.99999 , 2},
	  {5135 , 24931.3 , 1 , 1 , 1 , 0.99994 , 2},
	  {5136 , 46274.5 , 1 , 1 , 1 , 1 , 2},
	  {5137 , 158974 , 1 , 0.9884 , 1 , 1 , 2},
	  {5138 , 150458 , 1 , 1 , 1 , 0.999969 , 2},
	  {5139 , 171916 , 1 , 1 , 1 , 0.99999 , 2},
	  {5140 , 30341.5 , 1 , 1 , 1 , 1 , 2},
	  {5141 , 128178 , 1 , 1 , 1 , 0.999985 , 2},
	  {5142 , 99343.3 , 1 , 1 , 1 , 0.999965 , 2},
	  {5143 , 64562.4 , 1 , 1 , 1 , 0.999972 , 2},
	  {5145 , 299866 , 1 , 1 , 1 , 0.999995 , 2},
	  {5146 , 71350.6 , 1 , 1 , 1 , 0.999803 , 1},
	  {5147 , 24951.6 , 1 , 1 , 1 , 0.999791 , 1},
	  {5223 , 76972 , 1 , 1 , 1 , 0.99979 , 1},
	  {5224 , 221318 , 1 , 1 , 1 , 0.99999 , 2},
	  {5225 , 205837 , 1 , 1 , 1 , 0.99998 , 2},
	  {5226 , 107437 , 1 , 0.9804 , 1 , 1 , 2},
	  {6083 , 279844 , 1 , 1 , 1 , 1 , 2},
	  {6087 , 270723 , 1 , 0.9922 , 1 , 0.999993 , 2},
	  {6104 , 318198 , 1 , 1 , 1 , 1 , 2},
	  {6107 , 279925 , 1 , 0.9928 , 1 , 0.999985 , 2},
	  {6108 , 19404 , 1 , 1 , 1 , 1 , 2},
	  {6235 , 44945 , 1 , 1 , 1 , 0.999816 , 1},
	  {6236 , 141034 , 1 , 1 , 1 , 0.999982 , 2},
	  {6237 , 144619 , 1 , 0.975 , 1 , 0.999989 , 2},
	  {6238 , 84197.1 , 1 , 1 , 1 , 1 , 2},
	  {6239 , 143683 , 1 , 1 , 1 , 0.999991 , 2},
	  {6240 , 194026 , 1 , 1 , 1 , 0.999992 , 2},
	  {6258 , 181030 , 1 , 1 , 1 , 1 , 2},
	  {6259 , 141770 , 1 , 1 , 1 , 1 , 2},
	  {6260 , 125072 , 1 , 1 , 1 , 1 , 2},
	  {6261 , 73597.8 , 1 , 0.9714 , 1 , 0.999982 , 2},
	  {6262 , 56632.6 , 1 , 1 , 1 , 1 , 2}//,
	  //{6753 , 96829.3 , 1 , 1 , 1 , 0.999989 , 2}, // no file in folder
	  //{6754 , 66194.8 , 1 , 1 , 1 , 1 , 2},// no file in folder 
	  //{6688 , 38033.3 , 1 , 1 , 1 , 1 , 2},// no file in folder 
	  //{6826 , 66609.2 , 1 , 1 , 1 , 0.999986 , 2},// no file in folder 
	  //{6827 , 46462.2 , 1 , 1 , 1 , 1 , 2},// no file in folder 
	  //{6828 , 107415 , 1 , 1 , 1 , 0.999969 , 2}// no file in folder 
 };

	run_data_LT = {
	  {5098 , 125116 , 1 , 1 , 1 , 0.999978 , 2},
	  {5099 , 48850.1 , 1 , 1 , 1 , 1 , 2},
	  {5101 , 82435.3 , 1 , 1 , 1 , 0.999983 , 2},
	  {5102 , 129860 , 1 , 1 , 1 , 1 , 2},
	  {5103 , 137318 , 1 , 1 , 1 , 0.999968 , 2},
	  {5104 , 70640.7 , 1 , 1 , 1 , 1 , 2},
	  {5116 , 177216 , 1 , 1 , 1 , 0.999986 , 2},
	  {5117 , 283446 , 1 , 1 , 0.999993 , 0.999982 , 2},
	  {5118 , 147968 , 1 , 0.9833 , 1 , 0.999979 , 2},
	  {5119 , 130084 , 1 , 1 , 1 , 0.99999 , 2},
	  {5120 , 166292 , 1 , 1 , 0.999989 , 0.999959 , 2},
	  {5121 , 52510.6 , 1 , 1 , 1 , 1 , 2},
	  {5122 , 178954 , 1 , 1 , 1 , 0.99999 , 2},
	  {5123 , 127080 , 1 , 0.9825 , 1 , 1 , 2},
	  {5124 , 119308 , 1 , 1 , 1 , 0.999989 , 2},
	  {5125 , 160472 , 1 , 1 , 1 , 0.99999 , 2},
	  {5126 , 157809 , 1 , 1 , 1 , 0.99999 , 2},
	  {5131 , 172194 , 1 , 0.9865 , 0.999989 , 0.999969 , 2},
	  {5132 , 203323 , 1 , 1 , 1 , 0.999991 , 2},
	  {5197 , 162196 , 1 , 0.983 , 1 , 0.99999 , 2},
	  {5198 , 183049 , 1 , 1 , 1 , 0.99998 , 2},
	  {5199 , 175037 , 1 , 1 , 1 , 0.99999 , 2},
	  {5201 , 186047 , 1 , 1 , 0.999989 , 1 , 2},
	  {5202 , 164189 , 1 , 1 , 1 , 0.99999 , 2},
	  {5204 , 122758 , 1 , 1 , 1 , 1 , 2},
	  {5205 , 45309.2 , 1 , 1 , 1 , 0.99994 , 2},
	  {5206 , 72024.1 , 1 , 1 , 1 , 0.999984 , 2},
	  {5207 , 200510 , 1 , 1 , 1 , 0.999981 , 2},
	  {5208 , 196115 , 1 , 1 , 1 , 0.999981 , 2},
	  {5210 , 130771 , 1 , 1 , 1 , 1 , 2},
	  {5619 , 91211.7 , 1 , 0.963 , 1 , 0.999989 , 2},
	  {5620 , 112084 , 1 , 1 , 1 , 1 , 2},
	  {5621 , 102187 , 1 , 1 , 1 , 0.99999 , 2},
	  {5622 , 106492 , 1 , 1 , 1 , 1 , 2},
	  {5623 , 98277.5 , 1 , 1 , 1 , 1 , 2},
	  {5624 , 106928 , 1 , 0.9767 , 1 , 0.99998 , 2},
	  {5627 , 22240.7 , 1 , 1 , 1 , 1 , 2},
	  {5628 , 119267 , 1 , 1 , 1 , 0.999991 , 2},
	  {5629 , 105798 , 1 , 0.9762 , 1 , 1 , 2},
	  {5630 , 96216.9 , 1 , 1 , 1 , 0.999989 , 2},
	  {5631 , 85932.1 , 1 , 1 , 1 , 1 , 2},
	  {5632 , 103865 , 1 , 1 , 1 , 0.99999 , 2},
	  {5633 , 100989 , 1 , 1 , 1 , 0.999979 , 2},
	  {5636 , 95224.3 , 1 , 1 , 1 , 0.999991 , 2},
	  {5637 , 101567 , 1 , 1 , 1 , 0.999979 , 2},
	  {5638 , 42796.9 , 1 , 1 , 1 , 1 , 2},
	  {5673 , 198290 , 1 , 1 , 1 , 0.999958 , 2},
	  {5674 , 152771 , 1 , 1 , 1 , 0.99999 , 2},
	  {5675 , 218404 , 1 , 0.9881 , 1 , 0.999991 , 2},
	  {5676 , 206822 , 1 , 1 , 0.999989 , 0.99999 , 2},
	  {5677 , 172225 , 1 , 1 , 1 , 0.999956 , 2},
	  {5682 , 117152 , 1 , 1 , 1 , 1 , 2},
	  {5683 , 166859 , 1 , 1 , 1 , 0.999986 , 2},
	  {5684 , 231833 , 1 , 0.9888 , 1 , 1 , 2},
	  {5685 , 197851 , 1 , 0.9892 , 1 , 0.999964 , 2},
	  {5686 , 185338 , 1 , 0.9753 , 1 , 0.999987 , 2},
	  {5687 , 186875 , 1 , 1 , 1 , 0.999988 , 2},
	  {5690 , 219607 , 1 , 1 , 0.999989 , 0.99999 , 2},
	  {5691 , 184218 , 1 , 1 , 1 , 1 , 2},
	  {5692 , 226927 , 1 , 1 , 1 , 1 , 2},
	  {5693 , 258336 , 1 , 1 , 1 , 0.999972 , 2},
	  {5694 , 272122 , 1 , 1 , 1 , 0.999982 , 2},
	  {5696 , 239511 , 1 , 1 , 1 , 0.999971 , 2},
	  {5697 , 220704 , 1 , 1 , 1 , 1 , 2},
	  {5698 , 240087 , 1 , 1 , 1 , 0.999981 , 2},
	  {5699 , 240335 , 1 , 0.9878 , 1 , 0.999981 , 2},
	  {5700 , 194660 , 1 , 0.9787 , 1 , 0.999957 , 2},
	  {5703 , 213747 , 1 , 1 , 1 , 0.999979 , 2},
	  {5704 , 206072 , 1 , 1 , 1 , 0.99999 , 2},
	  {5705 , 230701 , 1 , 1 , 1 , 0.999989 , 2},
	  {5706 , 219102 , 1 , 1 , 1 , 0.999969 , 2},
	  {5707 , 220302 , 1 , 0.9889 , 1 , 0.999989 , 2},
	  {5708 , 117134 , 1 , 1 , 1 , 1 , 2},
	  {5711 , 209596 , 1 , 1 , 1 , 0.999989 , 2},
	  {5712 , 225028 , 1 , 1 , 1 , 1 , 2},
	  {5713 , 98571.9 , 1 , 1 , 1 , 1 , 2},
	  {5714 , 248222 , 1 , 0.9886 , 1 , 0.999981 , 2},
	  {5715 , 238301 , 1 , 1 , 1 , 1 , 2},
	  {5716 , 143142 , 1 , 0.98 , 1 , 0.999983 , 2}
};



	run_data_LT_Dummy = {

	  {5108 , 131242 , 1 , 1 , 1 , 1 , 2},
	  {5109 , 127375 , 1 , 1 , 1 , 1 , 2},
	  {5110 , 66760.1 , 1 , 1 , 1 , 1 , 2},
	  {5212 , 135193 , 1 , 1 , 1 , 0.99999 , 2},
	  {5213 , 147443 , 1 , 1 , 1 , 1 , 2},
	  {5214 , 131386 , 1 , 0.9844 , 1 , 0.999968 , 2},
	  {5625 , 111175 , 1 , 1 , 1 , 0.999977 , 2},
	  {5626 , 73428 , 1 , 1 , 1 , 0.999984 , 2},
	  {5634 , 134425 , 1 , 1 , 1 , 1 , 2},
	  {5635 , 107329 , 1 , 1 , 1 , 0.999987 , 2},
	  {5678 , 116644 , 1 , 1 , 1 , 0.999989 , 2},
	  {5679 , 73706.2 , 1 , 1 , 1 , 1 , 2},
	  {5680 , 33036.2 , 1 , 1 , 1 , 1 , 2},
	  {5688 , 122739 , 1 , 1 , 1 , 0.999977 , 2},
	  {5695 , 155562 , 1 , 1 , 1 , 1 , 2},
	  {5701 , 150135 , 1 , 1 , 1 , 0.99999 , 2},
	  {5702 , 34596.1 , 1 , 1 , 1 , 1 , 2},
	  {5709 , 127350 , 1 , 1 , 1 , 1 , 2},
	  {5710 , 61851.9 , 1 , 1 , 1 , 0.999955 , 2},
	  {5717 , 130208 , 1 , 0.9524 , 1 , 0.99999 , 2},
	  {5718 , 90618.5 , 1 , 1 , 1 , 1 , 2},
	  {5689 , 116881 , 1 , 1 , 1 , 1 , 2},
	  {6079 , 132556 , 1 , 1 , 1 , 1 , 2},
	  {6091 , 133230 , 1 , 1 , 1 , 1 , 2},
	  {6098 , 20325.3 , 1 , 1 , 1 , 1 , 2},
	  {6099 , 110825 , 1 , 1 , 1 , 1 , 2},
	  {6112 , 137394 , 1 , 1 , 1 , 0.999981 , 2},
	  {6246 , 84344.4 , 1 , 1 , 1 , 0.999988 , 2},
	  {6247 , 78411.7 , 1 , 1 , 1 , 1 , 2},
	  {6274 , 119447 , 1 , 1 , 1 , 0.99999 , 2},
	  {6275 , 146804 , 1 , 1 , 1 , 1 , 2},
	  {6616 , 109364 , 1 , 1 , 1 , 1 , 2},
	  {6617 , 123922 , 1 , 0.9524 , 1 , 1 , 2}/*,
	  {6662 , 39595.2 , 1 , 1 , 1 , 1 , 2},
	  {6663 , 99858.8 , 1 , 1 , 1 , 0.999967 , 2},
	  {6664 , 56546.3 , 1 , 0.96 , 1 , 0.999983 , 2},
	  {6699 , 71155.7 , 1 , 1 , 1 , 1 , 2},
	  {6700 , 16365.3 , 1 , 1 , 1 , 1 , 2},
	  {6707 , 87830.3 , 1 , 1 , 1 , 1 , 2},
	  {6714 , 95571.2 , 1 , 1 , 1 , 1 , 2},
	  {6726 , 74352.6 , 1 , 1 , 1 , 1 , 2},
	  {6740 , 31856.1 , 1 , 1 , 1 , 0.99999 , 2},
	  {6741 , 71104.5 , 1 , 1 , 1 , 1 , 2},
	  {6751 , 72032 , 1 , 1 , 1 , 1 , 2},
	  {6752 , 51775.6 , 1 , 1 , 1 , 0.999988 , 2},
	  {6804 , 94079.5 , 1 , 1 , 1.00001 , 1 , 2},
	  {6805 , 124828 , 1 , 1 , 1 , 1 , 2},
	  {6806 , 130606 , 1 , 0.983 , 1 , 0.999981 , 2},
	  {6807 , 108144 , 1 , 1 , 1 , 1 , 2},
	  {6808 , 103484 , 1 , 1 , 1 , 0.99999 , 2},
	  {6809 , 114936 , 1 , 1 , 0.999989 , 1 , 2},
	  {6867 , 97629.7 , 1 , 1 , 1 , 1 , 2},
	  {6868 , 110843 , 1 , 1 , 1 , 1 , 2},
	  {6869 , 92372.7 , 1 , 1 , 1 , 0.999989 , 2},
	  {6870 , 169991 , 1 , 1 , 1 , 1 , 2},
	  {6871 , 96457.8 , 1 , 1 , 1 , 1 , 2},
	  {6872 , 61665.6 , 1 , 1 , 1 , 1 , 2}*/
	  //for some reason I removed those, they are not in the data folder, take a closer look to those.
};
	}

else if ( momentum_spec_input=="3p57"   ){
	p_spec  = 3.57; 
	run_data_ST = {
	  {5228 , 172720 , 1 , 0.9963 , 1 , 0.99999 , 2},
	  {5229 , 163873 , 1 , 0.9978 , 1 , 1 , 2},
	  {5230 , 189468 , 1 , 0.9921 , 1 , 0.999981 , 2},
	  {5231 , 197821 , 1 , 0.9968 , 1 , 0.999981 , 2},
	  {5232 , 201402 , 1 , 0.9965 , 1 , 0.999991 , 2},
	  {5233 , 145434 , 1 , 0.9983 , 1 , 1 , 2},
	  {5234 , 122300 , 1 , 0.999 , 1 , 0.999984 , 2},
	  {5237 , 76622.9 , 1 , 0.9983 , 1.00002 , 0.999806 , 1},
	  {5291 , 69677.2 , 1 , 0.9983 , 1.00003 , 0.999701 , 1},
	  {5292 , 189082 , 1 , 0.9947 , 1 , 1 , 2},
	  {5293 , 187119 , 1 , 0.9987 , 1 , 1 , 2},
	  {5294 , 157591 , 1 , 0.9992 , 1 , 1 , 2},
	  {5295 , 42222.5 , 1 , 0.997 , 1 , 1 , 2},
	  {5296 , 232611 , 1 , 0.9958 , 1 , 0.999982 , 2},
	  {5297 , 257640 , 1 , 0.9977 , 1 , 1 , 2},
	  {6123 , 274734 , 1 , 0.9973 , 0.999992 , 0.999993 , 2},
	  {6126 , 246235 , 1 , 0.9973 , 1 , 1 , 2},
	  {6291 , 153109 , 1 , 0.9976 , 1 , 0.999991 , 2},
	  {6292 , 177021 , 1 , 0.9979 , 1 , 1 , 2},
	  {6293 , 182105 , 1 , 0.9987 , 0.99999 , 1 , 2},
	  {6295 , 1324.98 , 1 , 0.9231 , 1 , 1 , 2},
	  {6296 , 149312 , 1 , 0.9992 , 1 , 0.99999 , 2},
	  {6297 , 146948 , 1 , 0.9951 , 1 , 0.999975 , 2},
	  {6298 , 78951.2 , 1 , 0.9984 , 1 , 0.999828 , 1},
	  {6299 , 9823.74 , 1 , 1 , 1 , 1 , 2}
	  
	};
	run_data_LT = {
	  {5253 , 176995 , 1 , 0.9949 , 0.999988 , 0.999971 , 2},
	  {5254 , 93492 , 1 , 0.9969 , 1 , 0.999962 , 2},
	  {5255 , 120554 , 1 , 0.9971 , 1 , 1 , 2},
	  {5256 , 147099 , 1 , 0.9975 , 1 , 0.999989 , 2},
	  {5257 , 105876 , 1 , 0.996 , 0.999988 , 0.99999 , 2},
	  {5258 , 169951 , 1 , 0.9938 , 1 , 1 , 2},
	  {5259 , 25611.6 , 1 , 1 , 1 , 0.999927 , 2},
	  {5260 , 176887 , 1 , 0.9951 , 0.999988 , 0.999952 , 2},
	  {5261 , 222550 , 1 , 0.9943 , 1 , 0.999985 , 2},
	  {5262 , 49713.9 , 1 , 0.997 , 1 , 1 , 2}
	
	};
	run_data_LT_Dummy = {
	  {5251 , 132872 , 1 , 0.9959 , 1 , 0.99997 , 2},
	  {5268 , 93557.3 , 1 , 0.9981 , 1 , 1 , 2},
	  {5269 , 10814 , 1 , 1 , 1 , 1 , 2},
	  {5270 , 52427.4 , 1 , 1 , 1 , 1 , 2},
	  {6118 , 138566 , 1 , 0.9949 , 1 , 0.99999 , 2},
	  {6129 , 139432 , 1 , 0.9987 , 1 , 0.99999 , 2},
	  {6305 , 72374 , 1 , 0.9975 , 1 , 0.999989 , 2},
	  {6306 , 338.072 , 1 , 1 , 1 , 0.999978 , 2},
	  {6307 , 124440 , 1 , 0.9972 , 1 , 1 , 2},
	  {6308 , 127360 , 1 , 0.9957 , 1 , 1 , 2}
};
	}
	
else if ( momentum_spec_input=="3p09"   ){
	p_spec  = 3.09; 
	run_data_ST = {
	  {5298 , 7414.26 , 1 , 0.991 , 1 , 1 , 2},
	  {5299 , 144350 , 1 , 0.9963 , 1 , 0.999974 , 2},
	  {5300 , 152868 , 1 , 0.9966 , 1 , 0.999993 , 2},
	  {5301 , 192367 , 1 , 0.997 , 1 , 0.999982 , 2},
	  {5302 , 115652 , 1 , 0.9961 , 1 , 0.999979 , 2},
	  {5304 , 170660 , 1 , 0.9977 , 1.00001 , 0.999987 , 2},
	  {5305 , 157897 , 1 , 0.9964 , 1 , 0.999993 , 2},
	  {5306 , 115957 , 1 , 0.9964 , 1 , 0.99998 , 2},
	  {5307 , 105963 , 1 , 0.9973 , 1 , 0.999992 , 2},
	  {5308 , 197940 , 1 , 0.9975 , 1.00002 , 0.999989 , 2},
	  {5309 , 169863 , 1 , 0.9973 , 1.00001 , 0.999993 , 2},
	  {5310 , 102381 , 1 , 0.9975 , 0.999981 , 0.999989 , 2},
	  {5311 , 72293.1 , 1 , 0.9978 , 0.999899 , 0.999768 , 1},
	  {6365 , 56096.1 , 1 , 0.998 , 1 , 0.999849 , 1},
	  {6366 , 126862 , 1 , 0.9979 , 1 , 0.999982 , 2},
	  {6367 , 115415 , 1 , 0.9987 , 1 , 0.999994 , 2},
	  {6370 , 170180 , 1 , 0.9974 , 1 , 0.999989 , 2},
	  {6371 , 103591 , 1 , 0.9971 , 1 , 0.999987 , 2},
	  {6372 , 101073 , 1 , 0.9968 , 1 , 0.999984 , 2},
	  {6373 , 202199 , 1 , 0.998 , 0.999994 , 0.999992 , 2},
	  {6374 , 116021 , 1 , 0.9966 , 1 , 0.999993 , 2},
	  {6375 , 114948 , 1 , 0.9978 , 0.999989 , 0.999985 , 2},
	  {6376 , 132361 , 1 , 0.9976 , 0.99999 , 0.99998 , 2},
	  //{6380 , 138021 , 1 , 0.9975 , 1 , 0.999987 , 2}, // not in tape library
	  {6384 , 102812 , 1 , 0.9977 , 0.999988 , 0.999975 , 2},
	  {6385 , 150468 , 1 , 0.9962 , 1.00001 , 0.999982 , 2}



};
	run_data_LT = {

	  {5347 , 87789.1 , 1 , 0.997 , 1 , 0.999993 , 2},
	  {5349 , 137196 , 1 , 0.9971 , 1.00001 , 0.999967 , 2},
	  {5350 , 158338 , 1 , 0.997 , 1 , 0.999987 , 2},
	  {5351 , 135162 , 1 , 0.9973 , 0.99999 , 0.999959 , 2},
	  {5352 , 121306 , 1 , 0.9976 , 1 , 0.999983 , 2},
	  {5354 , 124254 , 1 , 0.9965 , 1 , 0.999973 , 2},
	  {5355 , 74997.9 , 1 , 0.9978 , 1 , 0.999955 , 2}	
};
	run_data_LT_Dummy = {
	  {5338 , 145704 , 1 , 0.998 , 1.00001 , 0.999984 , 2},
	  {5339 , 107444 , 1 , 0.9973 , 1.00001 , 0.99999 , 2},
	  {6345 , 117862 , 1 , 0.9967 , 1 , 1 , 2},
	  {6346 , 101462 , 1 , 0.9956 , 1 , 0.999976 , 2},
	  {6347 , 78534.8 , 1 , 0.9962 , 1 , 0.999985 , 2},
	  //{6348 , 78922.7 , 1 , 0.9961 , 1 , 1 , 2}, //not in tape library
	  {6545 , 80636.4 , 1 , 0.9967 , 1 , 1 , 2}



};
	}	
else if ( momentum_spec_input=="2p72"   ){
	p_spec  = 2.72; 
	run_data_ST = {
	  {5393 , 43478.1 , 1 , 0.9975 , 0.999962 , 0.999811 , 1},
	  {5394 , 126517 , 1 , 0.9977 , 1 , 0.999959 , 2},
	  {5395 , 149159 , 1 , 0.9979 , 0.999989 , 0.999958 , 2},
	  {5396 , 106622 , 1 , 0.9979 , 1.00001 , 0.999967 , 2}/*,
	  {5398 , 137733 , 1 , 0.9975 , 1 , 0.999956 , 2},
	  {5399 , 131493 , 1 , 0.9979 , 1 , 0.999958 , 2},
	  {5404 , 41185.8 , 1 , 0.9978 , 0.999989 , 0.999972 , 2},
	  {5405 , 64310.9 , 1 , 0.9979 , 1 , 0.999948 , 2},
	  {5408 , 37303.7 , 1 , 0.9961 , 1 , 0.999986 , 2},
	  {6386 , 22788.8 , 1 , 0.996 , 1 , 1 , 2},
	  {6387 , 130215 , 1 , 0.9972 , 1.00001 , 0.999972 , 2},
	  {6388 , 116782 , 1 , 0.997 , 1 , 0.999959 , 2},
	  {6389 , 100555 , 1 , 0.9975 , 1 , 0.999973 , 2},
	  {6390 , 120819 , 1 , 0.9971 , 1.00001 , 0.999971 , 2},
	  {6391 , 93823.1 , 1 , 0.9974 , 1 , 0.999971 , 2},
	  //{6393 , 137939 , 1 , 0.9972 , 1 , 0.999963 , 2},//not in tape library
	  //{6394 , 88791.2 , 1 , 0.997 , 1.00001 , 0.999971 , 2}, // not in tape library
	  {6395 , 125147 , 1 , 0.9965 , 1.00002 , 0.999981 , 2},
	  {6396 , 122199 , 1 , 0.9977 , 1.00001 , 0.999975 , 2},
	  //{6397 , 46543.2 , 1 , 0.9961 , 1 , 0.99995 , 2},//not in tape library
	  //{6398 , 50379.5 , 1 , 0.9963 , 0.999974 , 0.999817 , 1}//not in tape library*/
	  

};

	run_data_LT = {
	  {5356 , 106531 , 1 , 0.9972 , 1 , 0.999955 , 2},
	  {5357 , 131343 , 1 , 0.9972 , 0.99999 , 0.999934 , 2}/*,
	  {5358 , 38086.5 , 1 , 0.9977 , 1.00004 , 0.999932 , 2},
	  {5360 , 132089 , 1 , 0.9971 , 1.00001 , 0.99993 , 2},
	  {5361 , 35269.7 , 1 , 0.9971 , 1 , 0.999914 , 2},
	  {5363 , 12246.8 , 1 , 0.998 , 1 , 0.999991 , 2},
	  {5364 , 129717 , 1 , 0.9977 , 1 , 0.999938 , 2},
	  {5365 , 52677.6 , 1 , 0.9972 , 1 , 0.999919 , 2}*/
};
	run_data_LT_Dummy = {
	  {5369 , 104836 , 1 , 0.998 , 1 , 0.999942 , 2},
	  {5370 , 106162 , 1 , 0.9981 , 1.00001 , 0.999959 , 2}/*,
	  {6413 , 35025.5 , 1 , 0.9986 , 1 , 0.999988 , 2},
	  {6414 , 105743 , 1 , 0.998 , 1 , 0.999978 , 2},
	  {6415 , 79137.2 , 1 , 0.9973 , 1 , 0.999987 , 2},
	  {6416 , 103945 , 1 , 0.9972 , 0.999989 , 0.999959 , 2},
	  {6417 , 70449.6 , 1 , 0.9968 , 1 , 0.999994 , 2}*/
};
	}	
	
else if ( momentum_spec_input=="2p40"   ){
	p_spec  = 2.40; 
	run_data_ST = {
	  /*
//electrons
	  {5409 , 164352 , 1 , 0.9972 , 1.00001 , 0.999911 , 2},
	  {5410 , 145578 , 1 , 0.9976 , 1.00003 , 0.999925 , 2},
	  {5411 , 67022.9 , 1 , 0.9976 , 0.999975 , 0.999911 , 2},
	  {5412 , 54422.3 , 1 , 0.9976 , 1.00014 , 0.99967 , 1},
	  {6461 , 69610.2 , 1 , 0.9971 , 1.00001 , 0.999969 , 2},
	  {6463 , 125817 , 1 , 0.9972 , 1 , 0.999926 , 2},
	  {6464 , 108511 , 1 , 0.9974 , 1.00001 , 0.999958 , 2},
	  {6465 , 62894.6 , 1 , 0.9972 , 0.999954 , 0.999726 , 1}
	  */
	  {5028 , 19333 , 1 , 0.9875 , 1 , 1 , 3},                              // +2.40 C12
	  {6149 , 39070.6 , 1 , 1 , 0.999967 , 1 , 3},            // +2.40 C12
};
	run_data_LT = {
	  /*
//electrons
	  {5426 , 129662 , 1 , 0.9965 , 1.00001 , 0.999856 , 2},
	  {5427 , 82057 , 1 , 0.9965 , 1.00003 , 0.999858 , 2}
	  */
	  {5024 , 30061.1 , 1 , 0.9946 , 0.999912 , 1 , 3}            // +2.40 LD2
};
	run_data_LT_Dummy = {
	  {5422 , 90651.4 , 1 , 0.9966 , 1.00001 , 0.99994 , 2},
	  {5423 , 43042.7 , 1 , 0.9968 , 1 , 0.999944 , 2},
	  {6466 , 61827.2 , 1 , 0.997 , 0.999981 , 0.999933 ,2 },
	  {6467 , 77569 , 1 , 0.9968 , 1 , 0.999961 , 2}
};
	}	
	
else if ( momentum_spec_input=="2p11"   ){
	p_spec  = 2.11; 
	run_data_ST = {/*
//electrons
	  {5436 , 4143.43 , 2 , 0.9987 , 1.0094 , 0.999899 , 1},
	  {5437 , 34516.1 , 1 , 0.9968 , 0.99996 , 0.998725 , 1},
	  {5439 , 207715 , 1 , 0.997 , 1.00007 , 0.999864 , 2},
	  {6477 , 108650 , 1 , 0.9971 , 0.99998 , 0.999887 , 2},
	  {6478 , 42883.4 , 1 , 0.9977 , 1.00003 , 0.99986 , 2},
	  {6479 , 21989.5 , 1 , 0.9979 , 1.00026 , 0.999392 , 1}
		       */

	  {5032 , 22093.9 , 1 , 0.993 , 0.999967 , 0.999973 , 3},   // +2.11 C12
	  {5033 , 12005.9 , 1 , 1 , 1 , 0.99991 , 3},           // +2.11 C12
	  {5034 , 18670.1 , 1 , 0.9884 , 1 , 0.999963 , 3},           // +2.11 C12
	  {6156 , 67147.7 , 1 , 0.9923 , 0.999979 , 1 , 3}            // +2.11 C12
};
	run_data_LT = {
	  /*
//electrons
	  {5428 , 93390.8 , 1 , 0.9962 , 1.00007 , 0.999628 , 2}
	  */
	  {5045 , 38978.4 , 1 , 0.9972 , 1 , 0.999964 , 3}               // +2.11 LD2
};
	run_data_LT_Dummy = {
	  {5430 , 62124.8 , 1 , 0.9966 , 0.999982 , 0.999925 , 2},
	  {6481 , 56987.7 , 1 , 0.9972 , 1.00002 , 0.999904 , 2}

};
	}
	
else if ( momentum_spec_input=="1p85"   ){
	p_spec  = 1.85; 
	run_data_ST = {
	  /*
//electrons
	  {5440 , 114781 , 1 , 0.9966 , 1.00004 , 0.999701 , 2},
	  {5441 , 30951.7 , 2 , 0.9968 , 1.00732 , 0.999723 , 1},
	  {6487 , 12618.7 , 2 , 0.9973 , 0.997547 , 0.998848 , 1}
	  */
	  {5056 , 35350.3 , 1 , 0.9886 , 0.999959 , 1 , 3},           // +1.86 C12
	  {6173 , 30522.6 , 1 , 0.9956 , 1 , 1 , 3}            // +1.85 C12
};
	run_data_LT = {
	  /*
//electrons
	  {5449 , 50290.8 , 1 , 0.9965 , 0.999799 , 0.998578 , 2}
	  */
	  {5046 , 41920.7 , 1 , 0.9921 , 0.999964 , 0.999988 , 3}  // +1.86 LD2
};
	run_data_LT_Dummy = {
	  {5446 , 32965.1 , 1 , 0.997 , 0.999951 , 0.999767 , 2},
	  {6495 , 30917.1 , 1 , 0.9973 , 0.99996 , 0.999833 , 2}

};
	}	
else if ( momentum_spec_input=="1p63"   ){
	p_spec  = 1.63; 
	run_data_ST = {
	  /*
electrons
	  {5460 , 17600.2 , 5 , 0.9963 , 1.00311 , 0.998581 , 1},
	  {5461 , 68832.2 , 1 , 0.9968 , 0.999615 , 0.999147 , 2},
	  {6504 , 15811.7 , 1 , 0.9978 , 0.527468 , 0.470504 , 1}, 
	  {6505 , 66480.1 , 1 , 0.9976 , 1.0001 , 0.999489 , 2}
	  */
	  {5058 , 26703.3 , 1 , 0.9973 , 0.999902 , 1 , 3},           // +1.63 C12
	  {6196 , 20348.8 , 1 , 0.998 , 0.99998 , 1 , 3}            // +1.63  C12
};
	run_data_LT = {
	  /*
//electrons
	  {5450 , 45971.4 , 2 , 0.9958 , 1 , 0.999794 , 2},
	  {5452 , 18197 , 17 , 0.9958 , 1.03201 , 0.998074 , 1}
	  */
	  {5072 , 22602.1 , 1 , 0.9942 , 0.999973 , 0.999946 , 3},    // +1.63 LD2
	  {5073 , 51831.6 , 1 , 0.9958 , 0.999895 , 0.9999 , 3}        // +1.63 LD2
};
	run_data_LT_Dummy = {
	  {5455 , 17338.3 , 1 , 0.9968 , 1.00008 , 0.999652 , 2},
	  {5456 , 12208 , 2 , 0.996 , 1.00891 , 0.995621 , 1},
	  {6498 , 17266.3 , 1 , 0.9966 , 0.999951 , 0.999742 , 2}
};
	}	
	
else if ( momentum_spec_input=="1p44"   ){
	p_spec  = 1.44; 
	run_data_ST = {
          /*                                                                                                                                                                                                                                  
//electrons   
	  {5463 , 44335.9 , 1 , 0.9964 , 1.00018 , 0.999272 , 2},
	  {5465 , 20484.9 , 9 , 0.9957 , 1.00751 , 0.998634 , 1},
	  {6506 , 55768.9 , 1 , 0.9974 , 1.00033 , 0.999283 , 2},
	  {6507 , 19563.6 , 1 , 0.9977 , 0.346465 , 0.257273, 1} 
	  */
	  {5085 , 84673.2 , 1 , 0.996 , 0.999934 , 0.999919 , 3},    // +1.44 C12
	  {6213 , 41322.9 , 1 , 0.9965 , 0.99997 , 0.999989 , 3}     // +1.44 C12
};
	run_data_LT = {
          /*                                                                                                                                                                                                                                 
//electrons  
	  {5481 , 53019.8 , 3 , 0.9974 , 0.992916 , 0.998713, 2},
	  {5482 , 26267.9 , 33 , 0.997 , 1.03838 , 0.997874, 1}
	  */
	  {5074 , 87653 , 1 , 0.9956 , 0.999859 , 0.999884 , 3}       // +1.44 LD2
};
	run_data_LT_Dummy = {
	  {5477 , 9411.92 , 1 , 0.9973 , 1.00025 , 0.999582 , 2},
	  {5478 , 18208 , 5 , 0.9968 , 0.998397 , 0.998655 , 1},
	  {6509 , 20579 , 1 , 0.9974 , 1.00032 , 0.999618 , 2}
};
	}		
	
else if ( momentum_spec_input=="1p26"   ){
	p_spec  = 1.26; 
	run_data_ST = {
	  /*
//electrons
	  {5496 , 20156.6 , 17 , 0.9978 , 1.01497 , 0.997844 , 1},
	  {5497 , 52393.5 , 2 , 0.9976 , 1.00273 , 0.999863 , 2},
	  {6522 , 17517 , 17 , 0.9966 , 1.01453 , 0.99748, 1},
	  {6523 , 6162.97 , 1 , 0.9974 , 0.986569 , 0.985437 , 2},
	  {6524 , 39881.2 , 1 , 0.9971 , 0.999893 , 0.99882 , 2}
	  */

	  {5088 , 47313 , 1 , 0.9961 , 0.999874 , 0.999937 , 3},     // +1.26 C12
	  {6222 , 53152.7 , 1 , 0.9965 , 0.999931 , 0.999965 , 3}  // +1.26 C12
};
	run_data_LT = {
	  /*
//electrons
	  {5483 , 43465.1 , 5 , 0.9966 , 1.00719 , 0.998821 , 2},
	  {5484 , 23253.7 , 65 , 0.9954 , 1.07585 , 0.99979 , 1}
	  */
	  {5097 , 76009.8 , 2 , 0.9955 , 0.999481 , 0.999998 , 3}     // +1.26 LD2
};
	run_data_LT_Dummy = {
	  {5487 , 17181.7 , 1 , 0.9971 , 1.00009 , 0.9994 , 2},
	  {5488 , 20669.1 , 5 , 0.9971 , 0.991909 , 0.996106 , 1},
	  {6516 , 34295.7 , 9 , 0.9978 , 1.01433 , 0.999575 , 1},
	  {6517 , 25745.4 , 1 , 0.9975 , 1.00004 , 0.999464 , 2}
};
	}		
			
}
  if (  target_input=="C12" && angle_input=="26" && spectrometer=="SHMS" ){

    mc_file_lt_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_d2cryo22_hms.out";
    mc_file_st_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_12carbon22_hms.out";
    angle         = 26.0;
    AM_ST         = 12.0;
    density_ST    = 2.00;
    length_ST     = 0.287;
    delta_hi      = 30;
    delta_lo      = -20;
    f_st_mc_st    = Form("/work/smoran/xem2/sim/shms_26deg_carbon_%s.root"   , (const char*)momentum_spec_input) ;    
    f_lt_mc_st    = Form("/work/smoran/xem2/sim/shms_26deg_deuterium_%s.root", (const char*)momentum_spec_input) ;
    
    
// just for debuging    
if      ( momentum_spec_input=="1p95"   ){
	p_spec  = 1.95; 
	run_data_ST        = { 

//{17321 , 4843.7 , 1 , 0.9991 , 0.999736 , 0.996668}, //17321 , beam off for more than 20 min
{17322 , 30034.6 , 1 , 0.999 , 0.999396 , 0.998085} , 
{17424 , 45160.7 , 1 , 0.999 , 0.999433 , 0.998615}/*,
{17425 , 6242.37 , 1 , 0.9985 , 0.748087 , 0.173945}, // problems with this run, ASK!
{17426 , 2909.72 , 9 , 0.999 , 1.0075 , 0.96406}, // !Says EDTM low
{17427 , 7229.99 , 9 , 0.9992 , 1.01827 , 0.99018}*/// !Says EDTM incorrect
};  
	run_data_LT       = {  
//{17313 , 16969.1 , 65 , 0.999 , 1.11689 , 0.996942},
{17324 , 39436.9 , 3 , 0.9987 , 1.01177 , 0.99895}, 
{17432 , 28602.9 , 2 , 0.9986 , 1.00613 , 0.994859} };
   	run_data_LT_Dummy = {	
{17315 , 12938.9 , 9 , 0.9997 , 1.0178 , 0.99294},
{17323 , 15525.7 , 1 , 0.9991 , 1.00016 , 0.999532},
{17431 , 8562.53 , 2 , 0.9991 , 1.00151 , 0.999955}};
	}
/*else if      ( momentum_spec_input=="2p21"   ){
	p_spec  = 2.21; 
	run_data_ST        = {   };
	run_data_LT       = {    };
   	run_data_LT_Dummy = {	};
	}	
*/	
else if      ( momentum_spec_input=="2p52"   ){
	p_spec  = 2.52; 
	run_data_ST        = {  
{17422 , 5660.47 , 5 , 0.9993 , 0.97468 , 0.989297},
{17423 , 31107.6 , 1 , 0.9991 , 1.00066 , 0.999583} };
	run_data_LT       = { {17416 , 36204.5 , 2 , 0.9988 , 0.99779 , 0.998911}   };
   	run_data_LT_Dummy = {	{17419 , 8975.97 , 1 , 0.9993 , 1 , 0.999747} };
	}	
else if      ( momentum_spec_input=="2p86"   ){
	p_spec  = 2.86; 
	run_data_ST        = {  
{17407 , 36517.8 , 1 , 0.9992 , 1.00015 , 0.999702},
{17408 , 44856 , 1 , 0.9992 , 1.00231 , 1.00173} };
	run_data_LT       = {  {17415 , 36499.5 , 1 , 0.9989 , 1.00141 , 0.999179}  };
   	run_data_LT_Dummy = {{17413 , 18964.7 , 1 , 0.9993 , 1.00016 , 0.999812}	};
	}	
else if      ( momentum_spec_input=="3p25"   ){
	p_spec  = 3.25; 
	run_data_ST        = {  {17405 , 18935.2 , 2 , 0.9991 , 0.982341 , 0.999492},
{17406 , 91995.7 , 1 , 0.9992 , 1.00043 , 0.999832} };
	run_data_LT       = {  {17392 , 28911.8 , 1 , 0.9988 , 1.00087 , 0.999498}  };
   	run_data_LT_Dummy = {{17393 , 21737.6 , 1 , 0.999 , 1.0001 , 0.999936}	};
	}	
else if      ( momentum_spec_input=="3p69"   ){
	p_spec  = 3.69; 
	run_data_ST        = { {17382 , 185196 , 1 , 0.9991 , 1.00004 , 0.999896},
{17383 , 15537.6 , 1 , 0.9992 , 0.996674 , 0.992848},
{17390 , 30112.3 , 1 , 0.9991 , 1.00007 , 0.999909}  };
	run_data_LT       = {  {17387 , 24071.4 , 1 , 0.9989 , 1.00009 , 0.999741},
{17388 , 56362.4 , 1 , 0.9988 , 1.0001 , 0.999715}  };
   	run_data_LT_Dummy = {	{17386 , 52784.3 , 1 , 0.9993 , 1.00003 , 0.999881}};
	}	
else if      ( momentum_spec_input=="4p19"   ){
	p_spec  = 4.19; 
	run_data_ST        = { {17377 , 19376.2 , 1 , 0.9986 , 0.999256 , 0.996743},
{17378 , 197994 , 1 , 0.9991 , 1.00001 , 0.999903},
{17379 , 153962 , 1 , 0.9991 , 1.00002 , 0.999939},
{17380 , 198573 , 1 , 0.9992 , 1 , 0.999914} , {17381 , 175347 , 1 , 0.9991 , 1.00003 , 0.999924}};
	run_data_LT       = {  {17369 , 151316 , 1 , 0.9988 , 1.00014 , 0.999815},
{17370 , 112798 , 1 , 0.9988 , 1.00008 , 0.99982}  };
   	run_data_LT_Dummy = {{17371 , 235135 , 1 , 0.999 , 0.999988 , 0.999932}	};
	}					
	
else if      ( momentum_spec_input=="4p76"   ){
	p_spec  = 4.76; 
	run_data_ST        = { {17350 , 180980 , 1 , 0.9991 , 1.00003 , 0.99991},
{17351 , 155247 , 1 , 0.9991 , 0.999989 , 0.999934},
{17352 , 199213 , 1 , 0.9992 , 1.00001 , 0.999901},
{17353 , 182706 , 1 , 0.999 , 1.00004 , 0.999902},
{17354 , 24279.5 , 1 , 0.9996 , 1 , 0.997655}  };
	run_data_LT       = { {17364 , 175285 , 1 , 0.9988 , 0.999924 , 0.999862}   };
   	run_data_LT_Dummy = {{17360 , 23312.2 , 1 , 0.9993 , 1 , 0.999976},
{17362 , 107859 , 1 , 0.9986 , 1 , 0.999939},
{17363 , 21646.1 , 1 , 0.999 , 0.999933 , 0.999941}	};
	}	
	
else if      ( momentum_spec_input=="5p42"   ){
	p_spec  = 5.42; 
	run_data_ST        = { 
{17345 , 17206.6 , 1 , 0.9981 , 1.00073 , 0.999357},
{17346 , 172812 , 1 , 0.9993 , 0.999989 , 0.999932},
{17347 , 146809 , 1 , 0.9994 , 0.999979 , 0.999897},
{17349 , 62581.6 , 1 , 0.9985 , 0.99993 , 0.999925}  };
	run_data_LT       = { 
{17331 , 67572.9 , 1 , 0.9972 , 0.99996 , 0.999849},
{17332 , 87001.2 , 1 , 0.9994 , 1.00046 , 1.00081},
{17334 , 142651 , 1 , 0.9986 , 0.999968 , 0.999825},
{17335 , 11011.3 , 1 , 0.9987 , 0.999938 , 0.999899},
{17336 , 121833 , 1 , 0.9984 , 0.99994 , 0.999896},
{17337 , 44426.8 , 1 , 0.9976 , 0.699226 , 0.608817},
{17338 , 109881 , 1 , 0.9989 , 0.999896 , 0.999817},
{17339 , 793.706 , 1 , 1 , 1 , 1}   };
   	run_data_LT_Dummy = {	
{17340 , 99245.4 , 1 , 0.9985 , 0.999989 , 0.99996},
{17341 , 122878 , 1 , 0.9991 , 1 , 0.99993}};
	}	
	 
}


  if (  target_input=="C12" && angle_input=="26" && spectrometer=="HMS" ){

    mc_file_lt_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_d2cryo22_hms.out";
    mc_file_st_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_12carbon22_hms.out";
    density_ST    = 2.00;
    length_ST     = 0.287;
    delta_hi      = 20;
    delta_lo      = -20;
    f_st_mc_st    = Form("/work/smoran/xem2/sim/hms_26deg_carbon_%s.root"   , (const char*)momentum_spec_input) ;    
    f_lt_mc_st    = Form("/work/smoran/xem2/sim/hms_26deg_deuterium_%s.root", (const char*)momentum_spec_input) ;
    
if      ( momentum_spec_input=="1p95"   ){
	p_spec  = 1.95; 
	run_data_ST = {
	  {4592 , 32580.8 , 1 , 0.9965 , 0.999865 , 0.999833 , 3} // +1.95 C12

	  /*
//electrons
	  {4424 , 15799 , 1 , 0.9963 , 0.999451 , 0.997928 , 2},
	  {4426 , 43625.5 , 9 , 0.9964 , 1.00859 , 0.999534 , 1},
	  {4429 , 63653.6 , 1 , 0.9965 , 0.999858 , 0.998771 , 2},
	  {4432 , 26073 , 1 , 0.9967 , 0.99895 , 0.997507 , 2},
	  {4433 , 19886.8 , 1 , 0.9965 , 0.999628 , 0.998031 , 2},
	  {4437 , 14651.7 , 1 , 0.9963 , 0.999103 , 0.997681 , 2},
	  {4438 , 14539.3 , 5 , 0.9961 , 1.00299 , 0.996627 , 1},
	  {4439 , 1898.42 , 9 , 0.9968 , 1.0014 , 0.999898 , 1},
	  {4440 , 4656.45 , 9 , 0.9974 , 0.999637 , 0.999469 , 1},
	  {4441 , 20383.3 , 1 , 0.9966 , 0.999064 , 0.998097 , 2},
	  {4442 , 32780.9 , 1 , 0.9966 , 0.998799 , 0.997801 , 2},
	  {4443 , 21186 , 1 , 0.9965 , 1.00022 , 0.999101 , 2},
	  {4444 , 23948.5 , 1 , 0.9966 , 1.00028 , 0.999081 , 2},
	  {4445 , 23306.6 , 1 , 0.9963 , 0.999239 , 0.998203 , 2},
	  {4446 , 19300.2 , 1 , 0.9966 , 1.00027 , 0.999062 , 2},
	  {4447 , 15249.1 , 1 , 0.9963 , 1.0001 , 0.999058 , 2},
	  {4448 , 17352.6 , 1 , 0.9967 , 0.998788 , 0.997938 , 2},
	  {4449 , 16305.1 , 1 , 0.9966 , 0.999361 , 0.997912 , 2},
	  {4450 , 967.946 , 1 , 0.9963 , 0.982229 , 0.978693 , 2},
	  {4465 , 5975.45 , 2 , 0.9976 , 0.998538 , 0.99952 , 2},
	  {4466 , 29939.5 , 2 , 0.9967 , 1.00313 , 0.998648 , 2}
	  */
};
	run_data_LT = {

	  {4597 , 34774.8 , 1 , 0.9958 , 0.998857 , 0.999398 , 3} //  // +1.95 LD2
	  /*
//electrons
	  //{4421 , 2190.83 , 1 , 0.9945 , 0.130979 , 0.120518}, //too short, dont use it 
	  //{4422 , 3754.71 , 17 , 0.9972 , 1.07268 , 0.989397},//also too short
	  {4423 , 16760.7 , 5 , 0.9961 , 0.983526 , 0.999871 , 2},
	  {4428 , 38308.5 , 65 , 0.9958 , 1.02481 , 0.99996 , 1},
	  {4431 , 34362.9 , 5 , 0.9961 , 1.01808 , 0.999745 , 2},
	  {4436 , 13552.5 , 2 , 0.9961 , 1.00361 , 0.994886 , 2},
	  {4452 , 18795.8 , 3 , 0.9965 , 1.0008 , 0.998513 , 2},
	  {4453 , 15719 , 3 , 0.9964 , 0.990034 , 0.997631 , 2},
	  {4454 , 17757.3 , 3 , 0.9966 , 0.998158 , 0.997951 , 2},
	  {4457 , 13634.2 , 33 , 0.9965 , 1.03259 , 0.999796 , 1},
	  {4468 , 37975.4 , 9 , 0.9962 , 1.00091 , 0.999902 , 2}
	  */
};
	run_data_LT_Dummy = {
	  {4460 , 15859.1 , 5 , 0.9965 , 0.998509 , 0.996698 , 1},
	  {4467 , 18418.3 , 2 , 0.9965 , 0.997626 , 0.999896 , 2}
	
	};
	}
	
else if      ( momentum_spec_input=="2p21"   ){
	p_spec  = 2.21; 
	run_data_ST = {
	  {4585 , 37936.6 , 1 , 0.9975 , 0.999956 , 0.999855 , 3}
	  /*
//electrons
	  {4554 , 42701.1 , 2 , 0.9972 , 0.994233 , 0.999858 , 2},
	  {4555 , 25865.1 , 9 , 0.997 , 1.00255 , 0.999113 , 1}

	  */
};
	run_data_LT = {{4560 , 39675.7 , 5 , 0.9969 , 1.00085 , 0.998958 , 2}};
	run_data_LT_Dummy = {{4558 , 10216.3 , 1 , 0.9973 , 0.999866 , 0.999409 , 2}};
	}	
	
else if      ( momentum_spec_input=="2p52"   ){
	p_spec  = 2.52; 
	run_data_ST = {
	  {4572 , 49173.1 , 1 , 0.9959 , 0.999865 , 0.999977 , 3} // +2.52 C12
	  /*
//electrons
	  {4552 , 16200 , 5 , 0.997 , 1.00293 , 0.99779 , 1},
	  {4553 , 30643.7 , 1 , 0.9972 , 0.999335 , 0.998741 , 2}
	  */
};
	run_data_LT = {
	  /*
//electrons
{4547 , 33211.8 , 2 , 0.9972 , 0.997036 , 0.997302 , 2}
	  */
	  {4577 , 23139.8 , 1 , 0.9963 , 1 , 0.999762 , 3}               // +2.52 LD2
};
	run_data_LT_Dummy = {{4549 , 21197.5 , 1 , 0.9973 , 1.00006 , 0.999638 , 2}};
	}	
	
else if      ( momentum_spec_input=="2p86"   ){
	p_spec  = 2.86; 
	run_data_ST = {
	  /*
//electrons
	  {4539 , 23316.9 , 1 , 0.9975 , 0.999855 , 0.999652 , 2},
	  {4540 , 44543.6 , 1 , 0.9975 , 1.00021 , 0.999615 , 2},
	  {4541 , 14501.3 , 1 , 0.9971 , 0.994248 , 0.992931 , 1}
	  */
	  {4568 , 35918.7 , 1 , 0.9883 , 0.999869 , 0.99995 , 3}    // +2.86 C12
};
	run_data_LT = {
	  
//electrons
{4546 , 34184.6 , 1 , 0.9972 , 0.999751 , 0.998351 , 2}
	  


};
	run_data_LT_Dummy = {{4544 , 20590.7 , 1 , 0.9972 , 1.00046 , 0.999785 , 2}};
	}	
else if      ( momentum_spec_input=="3p25"   ){
	p_spec  = 3.25; 
	run_data_ST = {
	  {4537 , 22344.2 , 1 , 0.997 , 0.999893 , 0.999353 , 1},
	  {4538 , 90519 , 1 , 0.9973 , 1 , 0.999819 , 2}
};
	run_data_LT = {{4532 , 28473.9 , 1 , 0.9965 , 1.00024 , 0.999499 , 2}};
	run_data_LT_Dummy = {{4533 , 21628.6 , 1 , 0.9967 , 0.999952 , 0.999869 , 2}};
	}	
else if      ( momentum_spec_input=="3p69"   ){
	p_spec  = 3.69; 
	run_data_ST = {
	  {4522 , 192910 , 1 , 0.9967 , 1.00004 , 0.999894 , 2},
	  {4523 , 15964.7 , 1 , 0.9968 , 1 , 0.999638 , 1},
	  {4530 , 30292.6 , 1 , 0.9968 , 1 , 0.99991 , 2}
};
	run_data_LT = {
	  {4527 , 24280.3 , 1 , 0.996 , 1.00006 , 0.999796 , 2},
	  {4528 , 56853.4 , 1 , 0.9966 , 1.00003 , 0.999762 , 2}
};
	run_data_LT_Dummy = {{4526 , 52869.4 , 1 , 0.9962 , 1.00003 , 0.999908 , 2}};
	}	
else if      ( momentum_spec_input=="4p19"   ){
	p_spec  = 4.19; 
	run_data_ST = {

	  {4517 , 18681.5 , 1 , 0.9966 , 1 , 0.999651 , 1},
	  {4518 , 204153 , 1 , 0.9965 , 1 , 0.999928 , 2},
	  {4519 , 153666 , 1 , 0.9967 , 0.999989 , 0.999915 , 2},
	  {4520 , 196799 , 1 , 0.9966 , 1 , 0.999931 , 2},
	  {4521 , 171757 , 1 , 0.9968 , 1.00001 , 0.99992 , 2}
};
	run_data_LT = {
	  {4509 , 150529 , 1 , 0.9964 , 0.999989 , 0.999885 , 2},
	  {4510 , 114449 , 1 , 0.9962 , 1 , 0.999884 , 2}
};
	run_data_LT_Dummy = {{4511 , 234539 , 1 , 0.9967 , 1 , 0.999953 , 2}};
	}
else if      ( momentum_spec_input=="4p76"   ){
	p_spec  = 4.767; 
	run_data_ST = {
	  {4490 , 181429 , 1 , 0.9961 , 0.999989 , 0.999939 , 2},       
	  {4491 , 152909 , 1 , 0.9964 , 1.00001 , 0.999918 , 2},           
	  {4492 , 199090 , 1 , 0.9959 , 1 , 0.999949 , 2},              
	  {4493 , 178817 , 1 , 0.9953 , 1 , 0.999983 , 2},                  
	  {4494 , 24043.8 , 1 , 0.9968 , 1 , 0.999646 , 1}
};
	run_data_LT = {
	  {4504 , 170331 , 1 , 0.9953 , 1.00002 , 0.999954 , 2}
};
	run_data_LT_Dummy = {
	  {4502 , 108739 , 1 , 0.9959 , 1 , 0.999977 , 2},
	  {4503 , 25395.8 , 1 , 0.9966 , 1 , 0.999967 , 2}};
	}	
else if      ( momentum_spec_input=="5p42"   ){
	p_spec  = 5.42; 
	run_data_ST = {

	  {4485 , 17259 , 1 , 0.9811 , 1 , 0.999733 , 1},
	  {4486 , 172707 , 1 , 0.997 , 1.00001 , 0.999982 , 2},
	  {4487 , 146671 , 1 , 0.9943 , 1 , 0.999937 , 2},
	  {4489 , 64575.1 , 1 , 0.996 , 1 , 1 , 2}
};
	run_data_LT = {

	  {4474 , 67250.3 , 1 , 0.9937 , 1 , 0.999978 , 2},
	  {4475 , 86733.8 , 1 , 0.9981 , 1 , 0.99999 , 2},
	  {4476 , 159725 , 1 , 0.9937 , 1 , 1 , 2},
	  {4477 , 122309 , 1 , 0.9971 , 0.99999 , 0.999984 , 2},
	  {4478 , 157804 , 1 , 0.9961 , 1 , 0.999952 , 2}
};
	run_data_LT_Dummy = {
	  {4480 , 99111.1 , 1 , 0.9938 , 1 , 1 , 2},
	  {4481 , 122826 , 1 , 0.9967 , 1 , 1 , 2}};
	}	
		
}

  if (  target_input=="C12" && angle_input=="20" && spectrometer=="HMS"){

    mc_file_lt_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_d2cryo22_hms.out";
    mc_file_st_st = "/work/smoran/xem2/cross_sections/EXTERNALS_tables/xem2_emc_rc_12carbon22_hms.out";
    density_ST    = 2.00;
    length_ST     = 0.287;
    delta_hi      = 20;
    delta_lo      = -20;
    f_st_mc_st    = Form("/work/smoran/xem2/sim/hms_20deg_carbon_%s.root"   , (const char*)momentum_spec_input) ;    
    f_lt_mc_st    = Form("/work/smoran/xem2/sim/hms_20deg_deuterium_%s.root", (const char*)momentum_spec_input) ;
    //f_st_mc_st    = Form("/work/smoran/xem2/sim/test_offset_no_offset/test_offset_hms_20deg_carbon_%s.root"   , (const char*)momentum_spec_input) ;
    //f_lt_mc_st    = Form("/work/smoran/xem2/sim/test_offset_no_offset/test_offset_hms_20deg_deuterium_%s.root", (const char*)momentum_spec_input) ;

    if      ( momentum_spec_input=="6p60"   ){  
      p_spec  = 6.60;  
      run_data_ST = {

	{4689 , 151092 , 1 , 0.9972 , 1 , 0.999918 , 2},
	{4690 , 138433 , 1 , 0.9978 , 0.999992 , 0.999904 , 2},
	{4691 , 100914 , 1 , 0.9972 , 1 , 0.999812 , 2},
	{4692 , 5478.4 , 1 , 0.9971 , 1 , 0.999853 , 2},
	{4693 , 43179.7 , 1 , 0.9983 , 1 , 0.999863 , 2},
	{4694 , 134082 , 1 , 0.9978 , 0.99999 , 0.99986 , 2},
	{4695 , 105612 , 1 , 0.997 , 1 , 0.999915 , 2},
	{4696 , 70375.4 , 1 , 0.9974 , 1 , 0.999956 , 2},
	{4697 , 36921.1 , 1 , 0.9973 , 1 , 0.9999 , 2},
	{5809 , 135454 , 1 , 0.9967 , 1 , 0.999883 , 2},
	{5810 , 92077.7 , 1 , 0.9968 , 1 , 0.99992 , 2},
      };
      run_data_LT = {

	{4662 , 96044.6 , 1 , 0.9959 , 0.999989 , 0.999792 , 2},
	{4663 , 111752 , 1 , 0.9967 , 1 , 0.999804 , 2},
	{4664 , 120185 , 1 , 0.9971 , 0.999991 , 0.999795 , 2},
	{4665 , 93223.1 , 1 , 0.9973 , 0.999989 , 0.999836 , 2},
	{4666 , 40240.3 , 1 , 0.9971 , 1.00003 , 0.99972 , 2}
      };
      run_data_LT_Dummy = {
	{4668 , 114965 , 1 , 0.9983 , 1 , 0.99987 , 2},// commented out before, dont know why
	{4669 , 42155.6 , 1 , 0.9978 , 1 , 0.999876 , 2},  // !50K failed beecause of runinfo daemoon scrpit crash. Fixed after this
	{5795 , 60289.7 , 1 , 0.9953 , 1 , 0.999871 , 2}

};
}
    else if ( momentum_spec_input=="5p36"  ){  
      p_spec  = 5.36;  
      run_data_ST = {
	{4803 , 20343.3 , 1 , 0.9972 , 1.00019 , 0.99898 , 1},
	{4805 , 119391 , 1 , 0.9973 , 1.00005 , 0.999784 , 2},
	{4806 , 126830 , 1 , 0.9972 , 1 , 0.999776 , 2},
	{4808 , 135348 , 1 , 0.9972 , 1.00001 , 0.999798 , 2},
	{4809 , 114576 , 1 , 0.9972 , 1 , 0.999813 , 2},
	{4810 , 77990.9 , 1 , 0.9971 , 1.00007 , 0.999777 , 2},
      };

      run_data_LT = {
	{4786 , 23881.8 , 1 , 0.9974 , 1.00003 , 0.999701 , 2},
	{4787 , 25873.9 , 1 , 0.9973 , 1.00011 , 0.999693 , 2}
      };

      run_data_LT_Dummy = {
	{4789 , 38361.6 , 1 , 0.9974 , 0.999962 , 0.999857 , 2},
	{5845 , 74630.7 , 1 , 0.9967 , 1 , 0.999865 , 2}
};
}
    /*
//boiling runs
    else if ( momentum_spec_input=="4p00"   ){ 
      p_spec  = 4.00;  
      run_data_ST = {
	{5982 , 9132.81 , 1 , 0.9971 , 0.998838 , 0.998165 , 2},
	{5983 , 4221.02 , 1 , 0.9973 , 1.00025 , 0.999735 , 2},
	{5999 , 6510.65 , 1 , 0.9969 , 1.00006 , 0.999864 , 2},
	{6000 , 11177.9 , 1 , 0.9967 , 0.999936 , 0.99975 , 2},
	{6001 , 16737.9 , 1 , 0.9967 , 1.00016 , 0.999646 , 2},
	{6002 , 19465.9 , 1 , 0.9968 , 1.00016 , 0.999613 , 2},
	{6003 , 19707.9 , 1 , 0.9966 , 1.00006 , 0.999583 , 2},
	{6004 , 28932.3 , 1 , 0.9966 , 0.999952 , 0.99949 , 2},
	{6005 , 16130.9 , 1 , 0.9966 , 1.00033 , 0.999561 , 2},
	{6006 , 27529 , 1 , 0.9966 , 0.999607 , 0.998801 , 2},
	{5553 , 48724.3 , 1 , 0.9971 , 0.999036 , 0.99783 , 2},
	{5554 , 38357.3 , 1 , 0.9969 , 0.999426 , 0.99852 , 2},
	{5555 , 45727.8 , 1 , 0.9971 , 1.00003 , 0.999435 , 2},
	{5556 , 29824.4 , 1 , 0.997 , 1.0002 , 0.999507 , 2},
	{5557 , 22542.7 , 1 , 0.9971 , 1.00016 , 0.999591 , 2},
	{5558 , 14784.8 , 1 , 0.9973 , 0.999947 , 0.999715 , 2},
	{5560 , 7419.64 , 1 , 0.9974 , 0.999894 , 0.99985 , 2}
      };

      run_data_LT = {
	{5547 , 45441.8 , 2 , 0.9965 , 0.989496 , 0.997113 , 2},                    
	{5548 , 44512.4 , 2 , 0.9967 , 1.00075 , 0.997071 , 2},                              
	{5562 , 37314.4 , 2 , 0.9969 , 1.00226 , 0.998673 , 2},                   
	{5563 , 28621.8 , 2 , 0.9968 , 0.990964 , 0.999106 , 2},             
	{5565 , 23053.2 , 1 , 0.9969 , 0.998787 , 0.997443 , 2},        
	{5566 , 12940 , 1 , 0.9971 , 0.999167 , 0.998283 , 2},               
	{5567 , 7521.04 , 1 , 0.9972 , 1.00018 , 0.999627 , 2}
      };
      run_data_LT_Dummy = {
	{5576 , 21364.6 , 1 , 0.9971 , 0.999925 , 0.999625 , 2},
	{5974 , 31430.1 , 1 , 0.9971 , 1.00005 , 0.999638 , 2},
	{5992 , 25248.3 , 1 , 0.9965 , 1.00012 , 0.999646 , 2}
};


}*/
    else if ( momentum_spec_input=="5p87" ){  
      p_spec  = 5.878; 

      run_data_ST = {
	{5828 , 142189 , 1 , 0.9969 , 1 , 0.999881 , 2},
	{4755 , 113564 , 1 , 0.9972 , 0.999978 , 0.999843 , 2},
	{4757 , 154667 , 1 , 0.9975 , 1 , 0.999842 , 2},
	{4758 , 173195 , 1 , 0.9975 , 1 , 0.999859 , 2},
	{4759 , 127265 , 1 , 0.9973 , 1 , 0.999896 , 2},
	{4760 , 119810 , 1 , 0.9972 , 1.00003 , 0.99987 , 2},
	{4761 , 16675.4 , 1 , 0.9974 , 1 , 0.99901 , 1}
      };
      run_data_LT = {
	{4780 , 90719.6 , 1 , 0.9972 , 1.00004 , 0.999827 , 2},
	{4781 , 85806 , 1 , 0.9973 , 1.00001 , 0.999827 , 2},
	{4782 , 72493.4 , 1 , 0.9974 , 1.00002 , 0.999837 , 2},
	{4783 , 101834 , 1 , 0.9972 , 0.999956 , 0.999807 , 2},
	{4784 , 100494 , 1 , 0.9973 , 0.99999 , 0.999828 , 2},
	{4785 , 59009.2 , 1 , 0.9973 , 0.999983 , 0.999799 , 2}
      };
      run_data_LT_Dummy = {
	{4772 , 79655.2 , 1 , 0.9976 , 1.00004 , 0.999886 , 2},
	{4773 , 59520.9 , 1 , 0.9972 , 1.00004 , 0.999906 , 2},
	{4774 , 63688.8 , 1 , 0.9978 , 1.00002 , 0.999898 , 2},
	{4775 , 109281 , 1 , 0.9974 , 1.00001 , 0.999907 , 2},
	{4776 , 26074.2 , 1 , 0.9968 , 1 , 0.9999 , 2},
	{5836 , 68999.5 , 1 , 0.9963 , 1.00002 , 0.999887 , 2},
	{5964 , 31305.3 , 1 , 0.9965 , 1 , 0.99988 , 2}
};
}
    else if ( momentum_spec_input=="3p81" ){
      p_spec  = 3.81;
      run_data_ST = { 
	{4901 , 20563.4 , 1 , 0.9974 , 1 , 0.99911 , 2},
	{4902 , 28245.9 , 1 , 0.9975 , 0.999457 , 0.998705 , 2},
	{4903 , 18713.4 , 1 , 0.9975 , 0.995828 , 0.993058 , 1},
	{4904 , 5531.2 , 1 , 0.9978 , 0.999777 , 0.998021 , 1},
	{5894 , 30610 , 1 , 0.997 , 0.999045 , 0.997967 , 2},
	{6555 , 16882.1 , 1 , 0.9966 , 1.00014 , 0.999493 , 2}

      }; 
      run_data_LT = {  
	{4915 , 15917.3 , 1 , 0.9975 , 0.996461 , 0.994945,2}
      };
      run_data_LT_Dummy = {
	{4913 , 21040.3 , 1 , 0.9975 , 0.999876 , 0.999555 , 2},
	{5901 , 24606.6 , 1 , 0.9969 , 1.00017 , 0.999562 , 2}
};
    }

    else if ( momentum_spec_input=="4p78" ){
      p_spec  = 4.78 ;
      run_data_ST = { 
	{4846 , 97585 , 1 , 0.9973 , 0.999967 , 0.999571 , 2},
	{4847 , 82828 , 1 , 0.9972 , 1.0001 , 0.999717 , 2},
	{4848 , 70820.6 , 1 , 0.9972 , 1 , 0.999659 , 2},
	{4849 , 29273.3 , 1 , 0.9973 , 1.00009 , 0.999098 , 1},
	{5866 , 66385.2 , 1 , 0.9967 , 0.999939 , 0.99965 , 2},
	{5867 , 8107.99 , 1 , 0.9966 , 0.999782 , 0.999625 , 2}
      }; 
 
      run_data_LT = { 
	{4862 , 49207.7 , 1 , 0.9971 , 0.998802 , 0.997358,2},//!I assume this is the run that tripped the HV with 60 uA
	{4863 , 62358.6 , 1 , 0.9971 , 0.999231 , 0.99775,2}
      };

      run_data_LT_Dummy = {   
	{4859 , 19003.5 , 1 , 0.9976 , 0.99986 , 0.999789 , 2},
	{5873 , 33776 , 1 , 0.9963 , 1.00007 , 0.999775 , 2}
 };
   }

    else if ( momentum_spec_input=="4p27" ){
      p_spec  = 4.27 ;
      run_data_ST = { 
	{4875 , 89736.5 , 1 , 0.9973 , 0.999821 , 0.999327 , 2},
	{4876 , 18632.8 , 1 , 0.9974 , 0.99988 , 0.998932 , 1},
	{5883 , 24390.1 , 1 , 0.9968 , 1.00018 , 0.999213 , 2},
	{6554 , 29201.6 , 1 , 0.9966 , 1.00011 , 0.999463 , 2}

      }; 
      run_data_LT = { 
	{4864 , 23808.1 , 1 , 0.9973 , 0.999261 , 0.998408 , 2}

      };
      run_data_LT_Dummy = {
	{5876 , 32947.7 , 1 , 0.9965 , 1.00004 , 0.999665 , 2},
	{4866 , 14411.9 , 1 , 0.9976 , 1.00007 , 0.999687 , 2}
      };
    }
    else if ( momentum_spec_input=="3p40" ){
      p_spec  = 3.40;
      run_data_ST = { 
	{4600 , 33244.5 , 1 , 0.9885 , 1 , 0.999852 , 3},               // +3.40 C12
	{6008 , 33806.7 , 1 , 0.9903 , 0.999724 , 0.999895 , 3}         // +3.40 C12
	/*
//electrons
	{4925 , 23957.9 , 1 , 0.9973 , 0.998078 , 0.996092 , 2},
	{5910 , 35241.7 , 2 , 0.9969 , 0.992951 , 0.999385 , 2},
	{6569 , 8596.91 , 1 , 0.9972 , 0.999666 , 0.997689 , 2}
	*/
      }; 
 
      run_data_LT = { 
	{4599 , 65364.4 , 1 , 0.9919 , 0.999968 , 0.999724 , 3}  // +3.40 LD2
	/*
//electrons
	{4918 , 14734.9 , 2 , 0.9974 , 1.00253 , 0.997325,2}
	*/
      };

      run_data_LT_Dummy = {
	{5904 , 46425.1 , 1 , 0.9969 , 0.999795 , 0.999046 , 2}
	/*
//electrons
	{5904 , 46425.1 , 1 , 0.9969 , 0.999795 , 0.999046 , 2},
	{4920 , 20699.6 , 1 , 0.9974 , 0.9998 , 0.998774 , 2}
	*/
};

    }
    else if ( momentum_spec_input=="3p04" ){
      p_spec  =3.04 ;
      run_data_ST = { 
	/*
//electrons
	{4947 , 14395.6 , 5 , 0.9976 , 0.984388 , 0.995691 , 1},
	{4949 , 46227.9 , 2 , 0.9975 , 0.997463 , 0.999032 , 2},
	{5919 , 31012.7 , 2 , 0.9969 , 1.00294 , 0.997954 , 2},
	{5920 , 42401.9 , 2 , 0.9969 , 1.00407 , 0.998983 , 2}
	*/
	{4609 , 29983.1 , 1 , 0.9931 , 1 , 0.999849 , 3},         // +3.04 C12
	{5528 , 47554.6 , 1 , 0.9916 , 0.999913 , 0.999831 , 3},  // +3.04 C12
	{6027 , 74913.1 , 1 , 0.9909 , 0.999971 , 0.999846 , 3}   // +3.04 C12

      }; 
  
    run_data_LT = { 
      {4610 , 32254.4 , 1 , 0.9934 , 0.999792 , 0.999494 , 3},  // +3.04 LD2
      {5520 , 20584.6 , 1 , 0.9935 , 0.999946 , 0.99979 , 3}    // + 3.04 LD2
      /*
//electrons
      {4956 , 2760.81 , 3 , 0.9976 , 0.999951 , 0.988337 , 2},
      {4957 , 18887.9 , 3 , 0.9975 , 1.00042 , 0.997448 , 2}
      */
    };

    run_data_LT_Dummy = {
      {4954 , 12604.5 , 1 , 0.9974 , 0.997924 , 0.99535 , 2}
      /*
//electrons
      {4954 , 12604.5 , 1 , 0.9974 , 0.997924 , 0.99535 , 2},
      {5925 , 32672.8 , 1 , 0.9969 , 0.999013 , 0.997668 , 2}
      */
};

    }
    else if ( momentum_spec_input=="2p71" ){
      p_spec  = 2.71 ;
      run_data_ST = { 
	/*
//electrons
	{4969 , 27356.2 , 2 , 0.9977 , 0.998359 , 0.997904 , 2},
	{4970 , 13385.7 , 9 , 0.9975 , 1.00852 , 0.997008 , 1},
	{5935 , 39050.4 , 3 , 0.997 , 1.00903 , 0.99846 , 2},
	{6590 , 38932.2 , 2 , 0.9971 , 0.998937 , 0.998498 , 2}

	*/
	{4613 , 37492.5 , 1 , 0.9942 , 0.999952 , 0.999777 , 3},  // +2.71 C12
	{6045 , 45261 , 1 , 0.9934 , 0.999906 , 0.999738 , 3}     // +2.71 C12
      }; 
 
      run_data_LT = { 
	{4611 , 31247.6 , 1 , 0.9946 , 0.999645 , 0.999112 , 3}  // +2.71 LD2
	/*
//electrons
	{4958 , 18426.9 , 5 , 0.9974 , 1.00212 , 0.998242 , 2},
	{4962 , 17490.4 , 5 , 0.9973 , 1.0078 , 0.997956 ,2 }
	*/
      };
      run_data_LT_Dummy = {
	{4960 , 46473.3 , 2 , 0.9976 , 1.00287 , 0.999392 , 2}
	/*
//electrons
	{4960 , 46473.3 , 2 , 0.9976 , 1.00287 , 0.999392 , 2},
	{5929 , 27066.1 , 2 , 0.9969 , 0.995755 , 0.999132 , 2}
	*/

};

    }
    else if ( momentum_spec_input=="2p42" ){
      p_spec  =2.42 ;

      run_data_ST = { 
	{4641 , 41855 , 1 , 0.9952 , 0.999771 , 0.99964 , 3},       // +2.421 C12
	//{4642 , 34541.7 , 65 , 0.9949 , 1.0512 , 0.998875 , 1},    // +2.421 C12
	{6062 , 38613 , 1 , 0.9946 , 0.999729 , 0.99963 , 3}//,       // +2.42 C12
	//	{6063 , 41657.6 , 65 , 0.9957 , 1.05531 , 0.998857 , 1}  // +2.42 C12



	/*
//electrons
	{5004 , 11347.6 , 17 , 0.9974 , 0.997477 , 0.996712 , 1},
	{5005 , 15970.7 , 5 , 0.9975 , 0.990934 , 0.998845 , 2},
	{5006 , 43256.9 , 5 , 0.9974 , 1.00113 , 0.999396 , 2},
	{5949 , 34700.3 , 5 , 0.9969 , 0.993913 , 0.998839 , 2},
	{5950 , 17699.3 , 17 , 0.9965 , 1.01534 , 0.995359 , 1},
	{6591 , 16945.3 , 2 , 0.9971 , 1.00072 , 0.993626 , 2},
	{6592 , 18543.6 , 17 , 0.997 , 1.01001 , 0.998853 , 1}*/
      }; 
 
      run_data_LT = { 

	//{4634 , 32985.9 , 257 , 0.9901 , 1.14568 , 0.999808 , 1},// +2.421 LD2 
	{4636 , 36522.2 , 1 , 0.9944 , 0.999038 , 0.997509 , 3}  // +2.421 LD2


	/*
//electrons
	{5019 , 8647.4 , 33 , 0.997 , 1.05237 , 0.995873, 1},
	{5020 , 13782.4 , 9 , 0.997 , 1.0115 , 0.998482 , 2}

	*/
      };

      run_data_LT_Dummy = {

	{4622 , 20526.7 , 1 , 0.9952 , 0.999868 , 0.999847 , 3},  // +2.421 dummy
	//{4623 , 16950.7 , 65 , 1 , 1.02725 , 0.999902 , 1},       // +2.421 dummy
	{6071 , 24624.7 , 1 , 0.9939 , 0.999456 , 0.999825 , 3}//,  // +2.42 dummy 
	//{6073 , 23257.8 , 33 , 0.9943 , 1.02359 , 0.99879 , 1}    // +2.42 dummy

	/*
//electrons
	{5958 , 19936.3 , 2 , 0.997 , 0.997628 , 0.997547 , 2},
	{5959 , 12915.6 , 9 , 0.9966 , 1.00056 , 0.992086 , 1},
	{5015 , 9896.94 , 9 , 0.9976 , 1.00794 , 0.993851 , 1},
	{5016 , 14157.8 , 3 , 0.9975 , 1.0125 , 0.99876 , 2}
	*/
};


    }
}

  TString out_name = "cross_sections_" + target_input + "_angle_" + angle_input + "_momentum_" + momentum_spec_input + "_" + spectrometer;
 
  const double theta_central = angle*3.1415/180.0;

  Variables st_vars_mc, lt_vars_mc ;

  TFile *f_st_mc   = new TFile((const char*) f_st_mc_st  );
  TFile *f_lt_mc   = new TFile((const char*) f_lt_mc_st  );
  TTree *t_st_mc ,*t_lt_mc;
  if (spectrometer=="HMS" || spectrometer=="test"){
	  t_st_mc   = (TTree*) f_st_mc  ->Get("h1");
	  t_lt_mc   = (TTree*) f_lt_mc  ->Get("h1");
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
  if (var_dep == "delta") {  low_lim_hist=-8.1  ; high_lim_histo=8.1;} // if you are binning in delta
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
  // just for testing purposes:
  TH1F *X_MC_ST_Unweighted   = (TH1F*)X_MC_LT->Clone();
  TH1F *X_MC_LT_Unweighted   = (TH1F*)X_MC_LT->Clone();

  X_MC_ST  ->Sumw2();
  X_MC_LT  ->Sumw2();
  X_Data_ST->Sumw2();
  X_Data_LT->Sumw2();
  X_Data_LT_Dummy->Sumw2();
  X_MC_ST_Unweighted->Sumw2();
  X_MC_LT_Unweighted->Sumw2();
  cout<<__LINE__<<endl;
  vector<double>   V1_lt,V2_lt,V3_lt,V4_lt,V5_lt,V6_lt,V7_lt,V8_lt,V9_lt,V10_lt,V11_lt,V12_lt,V13_lt;
  vector<double>   V1_st,V2_st,V3_st,V4_st,V5_st,V6_st,V7_st,V8_st,V9_st,V10_st,V11_st,V12_st,V13_st;

  ImportRadcor(  V1_lt,V2_lt,V3_lt,V4_lt,V5_lt,V6_lt,V7_lt,V8_lt,V9_lt,V10_lt,V11_lt,V12_lt,V13_lt   ,  (const char*)mc_file_lt_st  );  
  ImportRadcor(  V1_st,V2_st,V3_st,V4_st,V5_st,V6_st,V7_st,V8_st,V9_st,V10_st,V11_st,V12_st,V13_st   ,  (const char*)mc_file_st_st  );

  const int size_st= V1_st.size();
  const int size_lt= V1_lt.size();
  cout<<"size of V1_ST: "<<size_st<< " , and V1_LT: "<< size_lt <<endl;
  cout<<__LINE__<<endl;
  TGraph2DErrors* gr2D_SigmaRad_ST  		= createGraph2D(V2_st, V3_st, V9_st, size_st);
  TGraph2DErrors* gr2D_SigmaBorn_ST 		= createGraph2D(V2_st, V3_st, V6_st, size_st);
  TGraph2DErrors* gr2D_SigmaRad_LT  		= createGraph2D(V2_lt, V3_lt, V9_lt, size_lt);
  TGraph2DErrors* gr2D_SigmaBorn_LT 		= createGraph2D(V2_lt, V3_lt, V6_lt, size_lt);
  TGraph2DErrors* gr2D_CoulombCorrection_ST     = createGraph2D(V2_st, V3_st, V13_st,size_st);
  cout<<__LINE__<<endl;
  Long64_t nentries_st_mc   = t_st_mc->GetEntries(); 
  cout<<__LINE__<<endl;
  Long64_t nentries_lt_mc   = t_lt_mc->GetEntries();
  cout<<__LINE__<<endl;
  double mc_scale_ST , mc_scale_LT;
  cout<<__LINE__<<endl;
  mc_scale_ST = calculate_mc_scale_factor( nentries_st_mc, p_spec, density_ST, length_ST, AM_ST, sim_charge_ST, delta_hi, delta_lo );
  mc_scale_LT = calculate_mc_scale_factor( nentries_lt_mc, p_spec, density_LT, length_LT, AM_LT, sim_charge_LT, delta_hi, delta_lo );
  cout<<__LINE__<<endl;
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
  }    


  ////////////////////
  /////   DATA   /////
  ////////////////////
  cout<<__LINE__<<endl;
    double total_SF_ST = 0.0;
    for (const auto& data : run_data_ST) { total_SF_ST += (data.charge*data.Eff_Fid*data.Eff_Elec_LiveTime*data.Eff_Comp_LiveTime) / data.PS; }
    double production_scale_ST = 1.0 / total_SF_ST;     cout<< "production_scale_ST : " << production_scale_ST <<endl;
    double total_SF_LT = 0.0;
    for (const auto& data : run_data_LT) { total_SF_LT += (data.charge*data.Eff_Fid*data.Eff_Elec_LiveTime*data.Eff_Comp_LiveTime) / data.PS; }
    double production_scale_LT = 1.0 / total_SF_LT;     cout<< "production_scale_LT : " << production_scale_LT <<endl;
    double total_SF_LT_Dummy = 0.0;
    for (const auto& data : run_data_LT_Dummy) { total_SF_LT_Dummy += (data.charge*data.Eff_Fid*data.Eff_Elec_LiveTime*data.Eff_Comp_LiveTime) / data.PS; }
    double production_scale_LT_Dummy = 1.0 / total_SF_LT_Dummy;     cout<< "production_scale_LT_Dummy : " << production_scale_LT_Dummy <<endl;
    cout<<__LINE__<<endl;
    double cuts_delta_min_st , cuts_delta_max_st , cuts_ytar_min_st , cuts_ytar_max_st , cuts_ph_min_st , cuts_ph_max_st , cuts_th_min_st , cuts_th_max_st , cuts_npeSum_st , cuts_etotTrack_st;
    double cuts_delta_min_lt , cuts_delta_max_lt , cuts_ytar_min_lt , cuts_ytar_max_lt , cuts_ph_min_lt , cuts_ph_max_lt , cuts_th_min_lt , cuts_th_max_lt , cuts_npeSum_lt , cuts_etotTrack_lt;
    double cuts_delta_min_st_mc, cuts_delta_max_st_mc, cuts_ytar_min_st_mc, cuts_ytar_max_st_mc, cuts_ph_min_st_mc, cuts_ph_max_st_mc, cuts_th_min_st_mc, cuts_th_max_st_mc, cuts_npeSum_st_mc, cuts_etotTrack_st_mc;
    double cuts_delta_min_lt_mc, cuts_delta_max_lt_mc, cuts_ytar_min_lt_mc, cuts_ytar_max_lt_mc, cuts_ph_min_lt_mc, cuts_ph_max_lt_mc, cuts_th_min_lt_mc, cuts_th_max_lt_mc, cuts_npeSum_lt_mc, cuts_etotTrack_lt_mc;

    std::ofstream outFile_runs_ST( out_name + "_yield_runs_ST.txt", std::ios_base::app );
    std::ofstream outFile_runs_LT( out_name + "_yield_runs_LT.txt", std::ios_base::app );
    std::ofstream outFile_runs_LT_Dummy( out_name + "_yield_runs_Al_LT.txt", std::ios_base::app );
    outFile_runs_ST       << "run\ttrigger\tYield\tError\n";
    outFile_runs_LT       << "run\ttrigger\tYield\tError\n";
    outFile_runs_LT_Dummy << "run\ttrigger\tYield\tError\n";

    for (int i = 0; i < run_data_ST.size(); i++) {
      if (spectrometer=="test"){spectrometer="HMS";}
    TFile *f_st_data = new TFile( Form("/work/smoran/xem2/data/%s/%s_replay_production_%d_-1.root" ,(const char*)spectrometer,(const char*)spectrometer_lowercase ,run_data_ST[i].run_number ));
    //cout<< Form("/cache/hallc/xem2/analysis/ONLINE/REPLAYS/HMS/PRODUCTION/hms_replay_production_%d_-1.root" , run_data_ST[i].run_number ) <<endl;
    TTree *t_st_data = (TTree*) f_st_data->Get("T");
    Variables st_vars_data;

    setBranchAddresses_data(t_st_data, st_vars_data,spectrometer);
    Long64_t nentries_st_data = t_st_data->GetEntries();
    cout<<"nentries_data : " << nentries_st_data <<endl;


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

    TH1F *X_tmp   = (TH1F*)X_MC_LT->Clone();
    	for (int j = 0; j < nentries_st_data; j++){
    		t_st_data->GetEntry(j);
		
		
if(st_vars_data.ptardp > cuts_delta_min_st && st_vars_data.ptardp < cuts_delta_max_st &&
   st_vars_data.decal  > cuts_etotTrack_st && st_vars_data.ptarth > cuts_th_min_st    &&
   st_vars_data.ptarth < cuts_th_max_st    && st_vars_data.ptarph > cuts_ph_min_st    &&
   st_vars_data.ptarph < cuts_ph_max_st    && st_vars_data.npeSum > cuts_npeSum_st    &&
   st_vars_data.ptary  < cuts_ytar_max_st  && st_vars_data.ptary  > cuts_ytar_min_st  && st_vars_data.ptarx < 2.0

  ){
		
      double Eprime ;
      Eprime =  p_spec*(1+0.01*st_vars_data.ptardp);
      double thetarad;
      if (spectrometer=="HMS" || spectrometer=="test")      {thetarad= TMath::ACos((cos(theta_central) + st_vars_data.ptarph*sin(theta_central))/TMath::Sqrt(1. + st_vars_data.ptarph*st_vars_data.ptarph+st_vars_data.ptarth*st_vars_data.ptarth));      }
      else if (spectrometer=="SHMS"){thetarad= TMath::ACos((cos(theta_central) - st_vars_data.ptarph*sin(theta_central))/TMath::Sqrt(1. + st_vars_data.ptarph*st_vars_data.ptarph+st_vars_data.ptarth*st_vars_data.ptarth));      }
      
      double thetadeg = thetarad*(180.0/3.14);                   
      double weight_F2 , var_x_axis;
      if (var_dep=="delta")  { var_x_axis = st_vars_data.ptardp; }
      else if (var_dep=="xb"){ var_x_axis = deltaToX( st_vars_data.ptardp , p_spec , Einitial , angle);  }
      //cout<< "var_x_axis : "<< var_x_axis <<endl;
      //cout<< "st_vars_data.ptardp : "<< st_vars_data.ptardp <<" , pspec: "<<p_spec << " Einitial: "<< Einitial<< " , angle: "<< angle << " , Xb: "<< var_x_axis<< endl;
      if (F2_option) { weight_F2 = weight_calculateF2( thetadeg ,  Eprime ,  Einitial)               ;}      else {  weight_F2 =1.0;}
		X_tmp->Fill( var_x_axis  , weight_F2   );

      		}
    	}
	auto result = calculateTotalEntriesAndError(X_tmp);
	double scaling_factor_run =  run_data_ST[i].PS / (run_data_ST[i].charge * run_data_ST[i].Eff_Fid * run_data_ST[i].Eff_Elec_LiveTime * run_data_ST[i].Eff_Comp_LiveTime );
	outFile_runs_ST << run_data_ST[i].run_number<< "\t"<< run_data_ST[i].trigger << "\t"<< result.first*scaling_factor_run << "\t" <<result.second*scaling_factor_run <<"\n" ;
    X_Data_ST->Add(X_tmp);
  }
    outFile_runs_ST.close();
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
	    	if( lt_vars_data.ptardp > cuts_delta_min_lt && lt_vars_data.ptardp < cuts_delta_max_lt &&
		    lt_vars_data.decal  > cuts_etotTrack_lt && lt_vars_data.ptarth > cuts_th_min_lt    &&
		    lt_vars_data.ptarth < cuts_th_max_lt    && lt_vars_data.ptarph > cuts_ph_min_lt    &&
		    lt_vars_data.ptarph < cuts_ph_max_lt    && lt_vars_data.npeSum > cuts_npeSum_lt    &&
		    lt_vars_data.ptary  < cuts_ytar_max_lt  && lt_vars_data.ptary  > cuts_ytar_min_lt && lt_vars_data.ptarx < 2.0){
		    double Eprime =  p_spec*(1+0.01*lt_vars_data.ptardp);
		    double thetarad;
      		    if (spectrometer=="HMS" || spectrometer=="test")      {thetarad= TMath::ACos((cos(theta_central) + lt_vars_data.ptarph*sin(theta_central))/TMath::Sqrt(1. + lt_vars_data.ptarph*lt_vars_data.ptarph+lt_vars_data.ptarth*lt_vars_data.ptarth));      }
                    else if (spectrometer=="SHMS"){thetarad= TMath::ACos((cos(theta_central) - lt_vars_data.ptarph*sin(theta_central))/TMath::Sqrt(1. + lt_vars_data.ptarph*lt_vars_data.ptarph+lt_vars_data.ptarth*lt_vars_data.ptarth));      }
		    double thetadeg = thetarad*(180.0/3.14); 
		    double weight_F2  , var_x_axis;
		    if (var_dep=="delta")  { var_x_axis = lt_vars_data.ptardp;  }
		    else if (var_dep=="xb"){ var_x_axis = deltaToX( lt_vars_data.ptardp , p_spec , Einitial , angle);  }
		    if (F2_option) { weight_F2 = weight_calculateF2( thetadeg ,  Eprime ,  Einitial)               ;}
		    else { weight_F2 =1.0;}
	      	X_tmp->Fill( var_x_axis , weight_F2 );
      		}
    	}
        auto result = calculateTotalEntriesAndError(X_tmp);
        double scaling_factor_run =  run_data_LT[i].PS / (run_data_LT[i].charge * run_data_LT[i].Eff_Fid * run_data_LT[i].Eff_Elec_LiveTime * run_data_LT[i].Eff_Comp_LiveTime );
        outFile_runs_LT << run_data_LT[i].run_number<< "\t"<<run_data_LT[i].trigger<< "\t"<< result.first*scaling_factor_run << "\t" <<result.second*scaling_factor_run <<"\n" ;

    X_Data_LT->Add(X_tmp);
     }

      outFile_runs_LT.close();

      //dummy
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
	     lt_vars_data_Dummy.ptardp > cuts_delta_min_lt && lt_vars_data_Dummy.ptardp < cuts_delta_max_lt &&
	     lt_vars_data_Dummy.decal  > cuts_etotTrack_lt && lt_vars_data_Dummy.ptarth > cuts_th_min_lt    &&
	     lt_vars_data_Dummy.ptarth < cuts_th_max_lt    && lt_vars_data_Dummy.ptarph > cuts_ph_min_lt    &&
	     lt_vars_data_Dummy.ptarph < cuts_ph_max_lt    && lt_vars_data_Dummy.npeSum > cuts_npeSum_lt    &&
	     lt_vars_data_Dummy.ptary  < cuts_ytar_max_lt  && lt_vars_data_Dummy.ptary  > cuts_ytar_min_lt && lt_vars_data_Dummy.ptarx <2.0 ){
	     
	     double Eprime =  p_spec*(1+0.01*lt_vars_data_Dummy.ptardp);
	     double thetarad;
      	     if (spectrometer=="HMS"||spectrometer=="test")      {thetarad= TMath::ACos((cos(theta_central) + lt_vars_data_Dummy.ptarph*sin(theta_central))/TMath::Sqrt(1. + lt_vars_data_Dummy.ptarph*lt_vars_data_Dummy.ptarph+lt_vars_data_Dummy.ptarth*lt_vars_data_Dummy.ptarth));      }
             else if (spectrometer=="SHMS"){thetarad= TMath::ACos((cos(theta_central) - lt_vars_data_Dummy.ptarph*sin(theta_central))/TMath::Sqrt(1. + lt_vars_data_Dummy.ptarph*lt_vars_data_Dummy.ptarph+lt_vars_data_Dummy.ptarth*lt_vars_data_Dummy.ptarth));      }
	     double thetadeg = thetarad*(180.0/3.14); 
	     double weight_F2 , var_x_axis;
	     if (var_dep=="delta")  { var_x_axis = lt_vars_data_Dummy.ptardp;  }
	     else if (var_dep=="xb"){ var_x_axis = deltaToX( lt_vars_data_Dummy.ptardp , p_spec , Einitial , angle);  }
	     if (F2_option) {  weight_F2 = weight_calculateF2( thetadeg ,  Eprime ,  Einitial)               ;}
	     else { weight_F2 =1.0;}
	     	     
	     X_tmp->Fill(var_x_axis, weight_F2);

	  }
        }

        auto result = calculateTotalEntriesAndError(X_tmp);
        double scaling_factor_run =  run_data_LT_Dummy[i].PS / (run_data_LT_Dummy[i].charge * run_data_LT_Dummy[i].Eff_Fid * run_data_LT_Dummy[i].Eff_Elec_LiveTime * run_data_LT_Dummy[i].Eff_Comp_LiveTime );
        outFile_runs_LT_Dummy << run_data_LT_Dummy[i].run_number<<"\t"<<run_data_LT_Dummy[i].trigger<< "\t"<< result.first*scaling_factor_run << "\t" <<result.second*scaling_factor_run <<"\n" ;

	X_Data_LT_Dummy->Add(X_tmp);
      }
      outFile_runs_LT_Dummy.close();

      std::cout<< "Production scale factors for ST, LT , and dummy: "<< production_scale_ST << " , "<< production_scale_LT << " , " << production_scale_LT_Dummy << std::endl;

      X_Data_ST->Scale(production_scale_ST);
      X_Data_LT->Scale(production_scale_LT);
      X_Data_LT_Dummy->Scale(production_scale_LT_Dummy);

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
      R =  thickness_LD2 / thickness_dummy;
      X_Data_LT_Dummy->Scale(R);

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


	if (lt_vars_mc.stop_id==0.0 && lt_vars_mc.hpsyptar < cuts_ph_max_lt_mc    && lt_vars_mc.hpsyptar > cuts_ph_min_lt_mc     &&
	                               lt_vars_mc.hpsxptar < cuts_th_max_lt_mc    && lt_vars_mc.hpsxptar > cuts_th_min_lt_mc     &&
				       lt_vars_mc.hpsytar  < cuts_ytar_max_lt_mc  && lt_vars_mc.hpsytar  > cuts_ytar_min_lt_mc ){

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
	  
	  double weight_ytarCorr , weight_MC_JacobianCorr , delta_new;
	  
	  if (ytar_corr_ON){weight_ytarCorr = weight_ytar_corr( lt_vars_mc.hpsytar , spectrometer);}
	  else weight_ytarCorr=1.0;
	  if (MC_Jacobian_corr_ON){weight_MC_JacobianCorr = weight_MC_Jacobian_corr ( lt_vars_mc.hpsxptar   , lt_vars_mc.hpsyptar,  spectrometer  );}
	  else weight_MC_JacobianCorr=1.0;
	  if (delta_corr_ON) {delta_new = weight_delta_corr(lt_vars_mc.hsdp , spectrometer );}
	  else delta_new = lt_vars_mc.hsdp;
	  double var_x_axis;
	  if      (var_dep=="delta")  { var_x_axis = delta_new;  }
	  else if (var_dep=="xb"   )  { var_x_axis = deltaToX( delta_new , p_spec , Einitial , angle);  }
	  
	  if (  delta_new < cuts_delta_max_lt_mc && delta_new > cuts_delta_min_lt_mc     ){	  X_MC_LT->Fill( var_x_axis,weight*weight_F2*weight_ytarCorr*weight_MC_JacobianCorr );   } 
	  
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

      cout<< "SIM ST:" <<endl;
      int weight_zero_ST_count=0;
      int progress_step_ST = nentries_st_mc / 100;
      for (int j=1; j<nentries_st_mc; j++) {
	t_st_mc->GetEntry(j);
	if (   st_vars_mc.stop_id==0.0  && 
				          st_vars_mc.hpsyptar < cuts_ph_max_st_mc    && st_vars_mc.hpsyptar > cuts_ph_min_st_mc     &&
            				  st_vars_mc.hpsxptar < cuts_th_max_st_mc    && st_vars_mc.hpsxptar > cuts_th_min_st_mc     &&
					  st_vars_mc.hpsytar  < cuts_ytar_max_st_mc  && st_vars_mc.hpsytar  > cuts_ytar_min_st_mc     ){
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
	  
	  double weight_ytarCorr , weight_MC_JacobianCorr , delta_new;
	  
	  if (ytar_corr_ON){weight_ytarCorr = weight_ytar_corr( st_vars_mc.hpsytar , spectrometer);}
	  else weight_ytarCorr=1.0;
	  if (MC_Jacobian_corr_ON){weight_MC_JacobianCorr = weight_MC_Jacobian_corr ( st_vars_mc.hpsxptar   , st_vars_mc.hpsyptar,  spectrometer  );}
	  else weight_MC_JacobianCorr=1.0;
	  if (delta_corr_ON) {delta_new = weight_delta_corr(st_vars_mc.hsdp , spectrometer );}
	  else delta_new = st_vars_mc.hsdp;
          double var_x_axis;
          if      (var_dep=="delta")  { var_x_axis = delta_new;  }
          else if (var_dep=="xb"   )  { var_x_axis = deltaToX( delta_new , p_spec , Einitial , angle);  }
	  //cout << "MC var_x_axis: " << var_x_axis << endl;
	  if ( delta_new < cuts_delta_max_st_mc && delta_new > cuts_delta_min_st_mc  ){X_MC_ST->Fill( var_x_axis , weight*weight_F2*weight_ytarCorr*weight_MC_JacobianCorr );}
	  
	  hsdp_values_ST.push_back(st_vars_mc.hsdp);
	  weight_values_ST.push_back(weight);

	  if (weight==0){  cout <<"E prime: " << Eprime << " , angle: "<< angle << " , delta is "<< st_vars_mc.hsdp << endl; weight_zero_ST_count=weight_zero_ST_count+1;}
	  //if (weight_zero_ST_count>3) { cout<< "exiting the code, because weight_zero_ST_count>3"<<endl;  exit(0);}
	}

	if (j % progress_step_ST == 0) {
	  int progress = (j * 100) / nentries_st_mc;
	  cout << "\rProgress: " << progress << "% completed , for ST" << flush;
	}
      }
      cout << "There were "<< weight_zero_ST_count << " zeroes in ST." <<endl;
      

      X_MC_ST->Scale(mc_scale_ST);
      X_MC_LT->Scale(mc_scale_LT);
      X_MC_ST_Unweighted->Scale(mc_scale_ST);
      X_MC_LT_Unweighted->Scale(mc_scale_LT);


      TCanvas *cc1;
      PlotWithRatio(cc1, X_Data_ST, X_Data_LT, X_MC_ST, X_MC_LT,X_MC_ST_Unweighted ,X_MC_LT_Unweighted , out_name);
      delete cc1; 
      plotHistograms_Dummy(X_Data_LT, X_Data_LT_Dummy, out_name);

      plottingWeightFactors( hsdp_values_ST, weight_values_ST ,hsdp_values_LT, weight_values_LT ,  out_name);

      cout<< "UNWEIGHTED:::\n\n\n";
      for (int i=0 ; i < 20 ; i++){  cout<<X_MC_ST_Unweighted->GetBinContent(1+i)<<endl;  }

  // DATA / MC , both properly scaled.

  X_Data_ST->Divide(X_MC_ST);
  X_Data_LT->Add(X_Data_LT_Dummy, -1); // performing dummy subtraction
  X_Data_LT->Divide(X_MC_LT);

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
  if (argc != 8) {
    std::cerr << "Usage: " << argv[0] << " <target_input> <angle_input> <momentum_spec_input> <spectrometer> <ytar_corr> <MCJacobian_corr> <delta_corr>" << std::endl;
    std::cerr << "Example: " << argv[0] << " C12 20 2p42 HMS ON ON ON" << std::endl;
    std::cerr << "Example of a TEST case: " << argv[0] << " C12 20 5p87 test ON ON ON" << std::endl;
    return 1;
  }

  cross_Sections_updated_workingProcess(argv[1], argv[2], argv[3], argv[4] , argv[5], argv[6], argv[7] );
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

void setBranchAddresses_data(TTree *tree, Variables &vars, TString spectrometer_option) {

//cout<< "INSIDE SETBRANCHADDRESS "<<spectrometer_option<<endl;

if (spectrometer_option=="SHMS"){cout<< "INSIDE SETBRANCHADDRESS INSIDE "<<spectrometer_option<<endl;}


if (spectrometer_option=="HMS" || spectrometer_option=="test"){
//cout<< "TAKING HMS BRANCHES"<<endl;
  tree->SetBranchAddress("H.gtr.dp",            &vars.ptardp);   //data delta    
  tree->SetBranchAddress("H.gtr.th",            &vars.ptarth);   //data target theta 
  tree->SetBranchAddress("H.gtr.ph",            &vars.ptarph);   //data target phi  
  tree->SetBranchAddress("H.gtr.x",             &vars.ptarx);    //data target x 
  tree->SetBranchAddress("H.gtr.y",             &vars.ptary);    //data target y  
  tree->SetBranchAddress("H.dc.x_fp",           &vars.dpsxfp) ;  //data focal plane x  
  tree->SetBranchAddress("H.dc.y_fp",           &vars.dpsyfp);   //data focal plane y  
  tree->SetBranchAddress("H.dc.xp_fp",          &vars.dpsxptar); //data focal plane xp   
  tree->SetBranchAddress("H.dc.yp_fp",          &vars.dpsyptar); //data focal plane yp    
  tree->SetBranchAddress("H.kin.W",             &vars.dW);       //data W   
  tree->SetBranchAddress("H.kin.Q2",            &vars.dQ2);      //data Q2  
  tree->SetBranchAddress("H.kin.nu",            &vars.dnu);      //data omega  
  tree->SetBranchAddress("H.kin.scat_ang_deg",  &vars.dtheta);   //data theta   
  tree->SetBranchAddress("H.cal.etottracknorm", &vars.decal);
  tree->SetBranchAddress("H.cer.npeSum",        &vars.npeSum);
}

else if (spectrometer_option=="SHMS"){
cout<< "TAKING SHMS BRANCHES"<<endl;
  tree->SetBranchAddress("P.gtr.dp",            &vars.ptardp);   //data delta    
  tree->SetBranchAddress("P.gtr.th",            &vars.ptarth);   //data target theta 
  tree->SetBranchAddress("P.gtr.ph",            &vars.ptarph);   //data target phi  
  tree->SetBranchAddress("P.gtr.x",             &vars.ptarx);    //data target x 
  tree->SetBranchAddress("P.gtr.y",             &vars.ptary);    //data target y  
  tree->SetBranchAddress("P.dc.x_fp",           &vars.dpsxfp) ;  //data focal plane x  
  tree->SetBranchAddress("P.dc.y_fp",           &vars.dpsyfp);   //data focal plane y  
  tree->SetBranchAddress("P.dc.xp_fp",          &vars.dpsxptar); //data focal plane xp   
  tree->SetBranchAddress("P.dc.yp_fp",          &vars.dpsyptar); //data focal plane yp    
  tree->SetBranchAddress("P.kin.W",             &vars.dW);       //data W   
  tree->SetBranchAddress("P.kin.Q2",            &vars.dQ2);      //data Q2  
  tree->SetBranchAddress("P.kin.nu",            &vars.dnu);      //data omega  
  tree->SetBranchAddress("P.kin.scat_ang_deg",  &vars.dtheta);   //data theta   
  tree->SetBranchAddress("P.cal.etottracknorm", &vars.decal);
  tree->SetBranchAddress("P.ngcer.npeSum",      &vars.npeSum);



}

}

double calculate_mc_scale_factor(int nentries_mc , double p_spec, double density , double length, double AM , double sim_charge , int delta_hi, int delta_lo ) {

  double domega           = (phi_hi - phi_lo)*(theta_hi-theta_lo) / 1000. /1000.; // differential solid angle in sr                  
  double ep_max           = p_spec*(1.+0.01*delta_hi);
  double ep_min           = p_spec*(1.+0.01*delta_lo);
  double deltaEprime      = ep_max-ep_min;
  double target_thickness = density*length;
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



void fillScaleFactors(const std::string& table_file, 
                      const std::string& target_input, 
                      const std::string& angle_input, 
                      const std::string& spectrometer, 
                      const std::string& momentum_spec_input, 
                      double p_spec,
                      std::vector<ScaleFactors>& run_data_LT, 
                      std::vector<ScaleFactors>& run_data_ST, 
                      std::vector<ScaleFactors>& run_data_LT_Dummy) {
  std::ifstream file(table_file);
  std::string line;

  if (file.is_open()) {
    // Clear existing data
    run_data_LT.clear();
    run_data_ST.clear();
    run_data_LT_Dummy.clear();

    while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#') {
	continue; // Skip empty lines or comments
      }

      std::stringstream ss(line);

      // Read the line and extract values
      std::string target, angle, spectr, momentum;
      double p_spec_in_file;
      ScaleFactors sf;
      ss >> target >> angle >> spectr >> momentum >> p_spec_in_file
	 >> sf.run_number >> sf.charge >> sf.PS >> sf.Eff_Fid 
	 >> sf.Eff_Elec_LiveTime >> sf.Eff_Comp_LiveTime >> sf.trigger;

      // Check if the current line matches the specified case
      if (target == target_input && angle == angle_input && spectr == spectrometer && 
	  momentum == momentum_spec_input && p_spec_in_file == p_spec) {
                
	// Determine which vector to fill based on the target
	if (target == "C12") {
	  run_data_ST.push_back(sf);
	} else if (target == "LD2") {
	  run_data_LT.push_back(sf);
	} else if (target == "Dummy") {
	  run_data_LT_Dummy.push_back(sf);
	}
      }
    }

    file.close();
  } else {
    std::cerr << "Unable to open file: " << table_file << std::endl;
  }
}











