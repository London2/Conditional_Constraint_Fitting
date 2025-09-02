#include <iostream>
#include <fstream>
#include <cmath>
#include "conditional_constraint.h"

  //Constructor
  Conditional_Constraint::Conditional_Constraint (int nbins, bool flag_reg, std::string file_reg, std::string path_MC, std::vector<std::string> v_fit_par_names) {
    (*this).configure(nbins, flag_reg, file_reg, path_MC, v_fit_par_names);
  }

  void Conditional_Constraint::configure (int nbins, bool flag_reg, std::string file_reg, std::string path_MC, std::vector<std::string> v_fit_par_names) {
    (*this).nbins = nbins;
    (*this).flag_reg = flag_reg;
    if (flag_reg) { Conditional_Constraint::read_reg(file_reg); }
    Conditional_Constraint::read_filelist(path_MC);
    Conditional_Constraint::resize_nbins();
    Conditional_Constraint::read_universes();
    Conditional_Constraint::setup_pars();
    Conditional_Constraint::set_fit_pars(v_fit_par_names);
  }

  //Build MC cov matrix and combine with data
  void Conditional_Constraint::build_cov () {

    //compute central values
    TMatrixD m_pred_cv(       (*this).nbins,    (*this).nbins);
    TMatrixD m_pred_cv_ext(   (*this).nbins_ext,(*this).nbins_ext);
    for (int ibin=0;ibin<(*this).nbins;ibin++) {
      for (int ithrow=0;ithrow<(*this).nthrows;ithrow++) { m_pred_cv(ibin,0) += (*this).v_pred[ithrow][ibin]; }
      m_pred_cv(    ibin,0) /= (*this).nthrows;
      m_pred_cv_ext(ibin,0)  = m_pred_cv(ibin,0);
    }
    std::vector<double> v_par_cv(npars,0);
    for (int i=0;i<(*this).npars;i++) {
      for (int j=0;j<(*this).nthrows;j++) { v_par_cv[i] += (*this).v_par_throw_vals[i][j]; }
      v_par_cv[i] /= (*this).nthrows;
      m_pred_cv_ext(nbins+i,0) = v_par_cv[i];
    }

    //compute MC cov
    for (int ibin=0;ibin<(*this).nbins;ibin++) {
      for (int jbin=0;jbin<(*this).nbins;jbin++) {
        (*this).cov_pred(ibin,jbin) = 0;
        for (int ithrow=0;ithrow<(*this).nthrows;ithrow++) { (*this).cov_pred(ibin,jbin) += ((*this).v_pred[ithrow][ibin]-m_pred_cv(ibin,0)) * ((*this).v_pred[ithrow][jbin]-m_pred_cv(jbin,0)); }
        (*this).cov_pred(ibin,jbin) /= (*this).nthrows;
      }
    }

    //compute cov between meas and parameters
    for (int ibin=0;ibin<(*this).nbins_ext;ibin++) {
      for (int jbin=0;jbin<(*this).nbins_ext;jbin++) {
        double ival, jval;
        (*this).cov_pred_ext(ibin,jbin) = 0;
        for (int ithrow=0;ithrow<(*this).nthrows;ithrow++) {
          if (ibin<nbins) { ival = (*this).v_pred[ithrow][ibin]-m_pred_cv(ibin,0); } else { ival = (*this).v_par_throw_vals[ibin-nbins][ithrow] - v_par_cv[ibin-nbins]; }
          if (jbin<nbins) { jval = (*this).v_pred[ithrow][jbin]-m_pred_cv(jbin,0); } else { jval = (*this).v_par_throw_vals[jbin-nbins][ithrow] - v_par_cv[jbin-nbins]; }
          (*this).cov_pred_ext(ibin,jbin) += ival * jval;
        }
        (*this).cov_pred_ext(ibin,jbin) /= (*this).nthrows;
      }
    }

    //combine data and MC cov matrices
    for (int ibin=0;ibin<(*this).nbins;ibin++) { for (int jbin=0;jbin<(*this).nbins;jbin++) { (*this).cov_data_ext(ibin,jbin) = (*this).cov_data(ibin,jbin); } }
    (*this).cov     = TMatrixD((*this).cov_data,     TMatrixD::kPlus, (*this).cov_pred);
    (*this).cov_ext = TMatrixD((*this).cov_data_ext, TMatrixD::kPlus, (*this).cov_pred_ext);
  }

  //Diagonalize covariance matrices
  void Conditional_Constraint::diagonalize_cov () {

    TMatrixDSym cov_sym((*this).nbins);
    for (int i=0;i<(*this).nbins;i++) { for (int j=0;j<(*this).nbins;j++) {
        cov_sym(i,j) = (*this).cov(i,j);
    } }
    TMatrixDSymEigen* cov_eigen = new TMatrixDSymEigen(cov_sym);
    (*this).cov_eigvals   = cov_eigen->GetEigenValues();
    (*this).cov_eigvecs   = cov_eigen->GetEigenVectors();      //eigenvectors are columns
    (*this).cov_eigvecs_T = cov_eigen->GetEigenVectors();
    (*this).cov_eigvecs_T.T();                                 //eigenvectors are rows here
    (*this).m_diff        = TMatrixD((*this).m_data,        TMatrixD::kMinus,(*this).m_pred_nom);
    (*this).m_diff_decomp = TMatrixD((*this).cov_eigvecs_T, TMatrixD::kMult, (*this).m_diff);
    for (int i=0;i<nbins;i++) { (*this).cov_decomp(i,i)    = (*this).cov_eigvals(i); }
    for (int i=0;i<nbins;i++) { (*this).m_chi2_decomp(i,0) = std::pow((*this).m_diff_decomp(i,0),2)/(*this).cov_decomp(i,i); }
  }

  //compute enlarged covariance using quantile mapping
  void Conditional_Constraint::scale_cov_QM () {
    Conditional_Constraint::diagonalize_cov ();
    std::map<double,double> map_chi2_quantiles = set_chi2_quantile_map((*this).m_chi2_decomp);
    for (int i=0;i<nbins;i++) { (*this).cov_decomp_QM(i,i) = (*this).cov_decomp(i,i) * std::max(1., (*this).m_chi2_decomp(i,0) / map_chi2_quantiles[(*this).m_chi2_decomp(i,0)]); }
    //convert back to standard basis
    (*this).cov      = TMatrixD( (*this).cov_eigvecs, TMatrixD::kMult,  TMatrixD((*this).cov_decomp_QM, TMatrixD::kMult, (*this).cov_eigvecs_T) );
    (*this).cov_data = TMatrixD( (*this).cov,      TMatrixD::kMinus, (*this).cov_pred );
    for (int i=0;i<(*this).nbins;i++) { for (int j=0;j<(*this).nbins;j++) { (*this).cov_ext(i,j) = (*this).cov(i,j); } }
  }

  void Conditional_Constraint::fit () {
    //Fit with all model parameters
    Conditional_Constraint::apply_constraint((*this).m_data,             (*this).m_pred_nom_ext, (*this).cov_ext,      false, (*this).m_pred_par_fit,     (*this).cov_par_fit);
    Conditional_Constraint::apply_constraint((*this).m_pred_par_fit,     (*this).m_pred_nom_ext, (*this).cov_pred_ext, true,  (*this).m_pred_xs_fit,      (*this).cov_xs_fit);
    //Fit with only selected parameters
    Conditional_Constraint::apply_constraint((*this).m_data_lim,         (*this).m_pred_nom_ext, (*this).cov_ext,      false, (*this).m_pred_par_lim_fit, (*this).cov_pred_par_lim_fit);
    Conditional_Constraint::apply_constraint((*this).m_pred_par_lim_fit, (*this).m_pred_nom_ext, (*this).cov_pred_ext, true,  (*this).m_pred_xs_lim_fit,  (*this).cov_pred_xs_lim_fit);
    (*this).m_pred_xs_lim_trim_fit = (*this).m_pred_xs_lim_fit.GetSub(0, (*this).nbins-1, 0, 0);
    //Calculate chi2
    (*this).chi2_nom     = Conditional_Constraint::compute_chi2((*this).m_data, (*this).m_pred_nom,             (*this).cov);
    (*this).chi2_fit     = Conditional_Constraint::compute_chi2((*this).m_data, (*this).m_pred_xs_fit,          (*this).cov);
    (*this).chi2_lim_fit = Conditional_Constraint::compute_chi2((*this).m_data, (*this).m_pred_xs_lim_trim_fit, (*this).cov);
  }

  //read in Ac txt file, select correct rows, fill m_reg
  void Conditional_Constraint::read_reg (std::string file_reg) {
    TFile* f_reg = TFile::Open(file_reg.c_str());
    TMatrixD* m_reg = (TMatrixD*)f_reg->Get("m_reg");
    (*this).m_reg.ResizeTo((*this).nbins,(*this).nbins);
    for (int i=0;i<(*this).nbins;i++) { for (int j=0;j<(*this).nbins;j++) { (*this).m_reg(i,j) = (*m_reg)(i,j); } }
  }

  //read in MC filelist
  void Conditional_Constraint::read_filelist (std::string path_MC) {
    (*this).filelist_MC.clear();
    std::ifstream infile_MC(path_MC+"filelist_throws.txt");
    std::string line_MC;
    if (infile_MC.is_open()) {
      while ( getline (infile_MC,line_MC) ) { (*this).filelist_MC.push_back(path_MC+line_MC); }
      infile_MC.close();
    }
    (*this).nfiles = (*this).filelist_MC.size();
  }

  void Conditional_Constraint::resize_nbins () {
    (*this).cov_data.ResizeTo(                 (*this).nbins,(*this).nbins);
    (*this).cov_pred.ResizeTo(                 (*this).nbins,(*this).nbins);
    (*this).cov.ResizeTo(                      (*this).nbins,(*this).nbins);
    (*this).cov_decomp.ResizeTo(               (*this).nbins,(*this).nbins);
    (*this).cov_xs_fit.ResizeTo(               (*this).nbins,(*this).nbins);
    (*this).cov_decomp_QM.ResizeTo(            (*this).nbins,(*this).nbins);
    (*this).cov_eigvecs.ResizeTo(              (*this).nbins,(*this).nbins);
    (*this).cov_eigvecs_T.ResizeTo(            (*this).nbins,(*this).nbins);
    (*this).cov_eigvals.ResizeTo(              (*this).nbins);
    (*this).m_data.ResizeTo(                   (*this).nbins,1);
    (*this).m_pred_nom.ResizeTo(               (*this).nbins,1);
    (*this).m_pred_xs_fit.ResizeTo(            (*this).nbins,1);
    (*this).m_pred_xs_lim_trim_fit.ResizeTo(   (*this).nbins,1);
    (*this).m_diff.ResizeTo(                   (*this).nbins,1);
    (*this).m_diff_decomp.ResizeTo(            (*this).nbins,1);
    (*this).m_chi2_decomp.ResizeTo(            (*this).nbins,1);
  }

  //read in universes
  void Conditional_Constraint::read_universes () {
    (*this).v_pred.clear();
    (*this).map_par_nom_val.clear();
    (*this).map_par_throw_vals.clear();
    //iterate files
    (*this).nthrows = 0;
    for (int ifile=0;ifile<(*this).nfiles;ifile++) {
      TFile* f_universe = TFile::Open((*this).filelist_MC[ifile].c_str());
      TList* key_dirs = f_universe->GetListOfKeys();
      //iterate folders within each file
      for (auto key_dir: *key_dirs) {
        std::string dir_name = key_dir->GetName();
        //get data, cov, and nominal MC
        if ((*this).nthrows==0 && dir_name=="nominal") {
          TH1D* h_MC       = (TH1D*)f_universe->Get("nominal/MicroBooNE_CC1Mu0pNp_XSec_XpEMuCosThetaMu_nu_MC");
          TH1D* h_data     = (TH1D*)f_universe->Get("nominal/MicroBooNE_CC1Mu0pNp_XSec_XpEMuCosThetaMu_nu_data");
          TH2D* h_data_cov = (TH2D*)f_universe->Get("nominal/MicroBooNE_CC1Mu0pNp_XSec_XpEMuCosThetaMu_nu_COV");
	  for (int ibin=0;ibin<(*this).nbins;ibin++) {
	    (*this).m_data(ibin,0)     = h_data->GetBinContent(ibin+1)*1e38;
	    (*this).m_pred_nom(ibin,0) = h_MC->GetBinContent(ibin+1)*1e38;
            for (int jbin=0;jbin<(*this).nbins;jbin++) { (*this).cov_data(ibin,jbin) = h_data_cov->GetBinContent(ibin+1,jbin+1); } //already in 1e38 units
          }
	  if ((*this).flag_reg) { (*this).m_pred_nom = TMatrixD((*this).m_reg, TMatrixD::kMult, (*this).m_pred_nom); }
        }
        //get MC throw, go into directory
        if ( dir_name=="nominal" || dir_name=="postfit" || dir_name=="error_iterations" ) { continue; }
        f_universe->cd();
        f_universe->cd(dir_name.c_str());
        (*this).v_pred.push_back({});
        TMatrixD m_pred((*this).nbins,1);
        TH1D* h_MC = (TH1D*)gDirectory->Get("MicroBooNE_CC1Mu0pNp_XSec_XpEMuCosThetaMu_nu_MC");
        for (int ibin=0;ibin<(*this).nbins;ibin++) { m_pred(ibin,0) = h_MC->GetBinContent(ibin+1)*1e38; }
        if ((*this).flag_reg) { m_pred = TMatrixD(m_reg, TMatrixD::kMult, m_pred); }
        for (int ibin=0;ibin<(*this).nbins;ibin++) { (*this).v_pred.back().push_back(m_pred(ibin,0)); }
        //Get model parameter values
        TList* key_pars = gDirectory->GetListOfKeys();
        for (auto key_par: *key_pars) {
          std::string par_name = key_par->GetName();
          if ( (*this).nthrows==0 && par_name.find("_t_data") != std::string::npos ) { (*this).map_par_nom_val[par_name] =             ((TH1D*)gDirectory->Get(par_name.c_str()))->GetBinContent(1); }
          if (                       par_name.find("_t_MC")   != std::string::npos ) { (*this).map_par_throw_vals[par_name].push_back( ((TH1D*)gDirectory->Get(par_name.c_str()))->GetBinContent(1) ); }
        }
        (*this).nthrows++;
      } //throws
    } //file
  }

  //Get list of parameters by iterating over map
  void Conditional_Constraint::setup_pars () {
    (*this).v_par_throw_vals.clear();
    (*this).v_par_names.clear();
    (*this).v_par_nom.clear();
    (*this).npars = (*this).map_par_nom_val.size();
    for (std::map<std::string,std::vector<double>>::iterator it = (*this).map_par_throw_vals.begin(); it != (*this).map_par_throw_vals.end(); it++) {
      std::string par_name = it->first;
      std::string par_name_nom = par_name.substr(0,par_name.find("_t_MC"))+"_t_data";
      std::vector<double> par_vals = it->second;
      (*this).v_par_throw_vals.push_back(par_vals);
      (*this).v_par_names.push_back(par_name);
      (*this).v_par_nom.push_back((*this).map_par_nom_val[par_name_nom]);
    }
  }

  //Count and arrange list of parameters
  void Conditional_Constraint::set_fit_pars (std::vector<std::string> v_fit_par_names) {
    (*this).nfit = v_fit_par_names.size();
    (*this).v_fit_par_names.clear();
    for (int i=0;i<(*this).nfit;i++) { (*this).v_fit_par_names.push_back(v_fit_par_names[i]); }
    //Reorder parameter list so that fit parameters are last (usefull for conditional constraint later)
    if ((*this).nfit!=0) {
      for (int i=(*this).npars-1;i>=0;i--) {  //go in descending order because we are modifying the back of the list as we iterate over it
        std::string par_name = (*this).v_par_names[i];
        for (int j=0;j<(*this).nfit;j++) { 
          if (par_name==(*this).v_fit_par_names[j]) {
            (*this).v_par_names.push_back(par_name);
            (*this).v_par_throw_vals.push_back((*this).v_par_throw_vals[i]);
            (*this).v_par_names.erase((*this).v_par_names.begin()+i);
            (*this).v_par_throw_vals.erase((*this).v_par_throw_vals.begin()+i);
          }
        }
      }
    } else { (*this).nfit = (*this).npars; }
    Conditional_Constraint::setup_joint_measurement();
  }

  //Extend number of bins to include measurement bins and parameter bins
  void Conditional_Constraint::setup_joint_measurement() {
    (*this).nbins_ext = (*this).nbins +   (*this).npars;
    (*this).m_pred_nom_ext.ResizeTo(      (*this).nbins_ext,1);
    (*this).m_data_lim.ResizeTo(          (*this).nbins_ext-(*this).nfit,1);    //exclude selected parameters to be fit
    (*this).m_pred_par_fit.ResizeTo(      (*this).npars,1);
    (*this).m_pred_par_lim_fit.ResizeTo(  (*this).nfit,1);
    (*this).m_pred_xs_lim_fit.ResizeTo(   (*this).nbins_ext-(*this).nfit,1);
    (*this).cov_data_ext.ResizeTo(        (*this).nbins_ext,(*this).nbins_ext);
    (*this).cov_pred_ext.ResizeTo(        (*this).nbins_ext,(*this).nbins_ext);
    (*this).cov_ext.ResizeTo(             (*this).nbins_ext,(*this).nbins_ext);
    (*this).cov_pred_xs_lim_fit.ResizeTo( (*this).nbins_ext-(*this).nfit,(*this).nbins_ext-(*this).nfit);
    (*this).cov_pred_par_lim_fit.ResizeTo((*this).nfit, (*this).nfit);
    (*this).cov_par_fit.ResizeTo(         (*this).npars,(*this).npars);
    for (int i=0;i<(*this).nbins;i++) {
      (*this).m_pred_nom_ext(i,0) = (*this).m_pred_nom(i,0);
      (*this).m_data_lim(    i,0) = (*this).m_data(i,0);
    }
    for (int i=(*this).nbins;i<(*this).nbins_ext;i++)              { (*this).m_pred_nom_ext(i,0) = (*this).v_par_nom[i-nbins]; }
    for (int i=(*this).nbins;i<(*this).nbins_ext-(*this).nfit;i++) { (*this).m_data_lim(i,0)     = (*this).v_par_nom[i-nbins]; }
  }

  //helper function that constructs a TGraph that maps from existing data-pred chi2 values to adjusted values that better follow a chi2 distribution
  //uses method of quantile mapping (also called probability integral transform)
  std::map<double,double> Conditional_Constraint::set_chi2_quantile_map (TMatrixD m_chi2_meas) {
    int nbins_meas = m_chi2_meas.GetNrows();
    std::vector<double> v_chi2_val;  //start with a chi2 value of 0 that will be assigned to CDF=0
    for (int i=0;i<nbins_meas;i++) { v_chi2_val.push_back(m_chi2_meas(i,0)); }
    std::sort(v_chi2_val.begin(), v_chi2_val.end());

    std::map<double,double> map_chi2_quantiles;
    for (int i=0;i<nbins_meas;i++) {
      double cdf = (i+0.5)/nbins_meas;
      map_chi2_quantiles[v_chi2_val[i]] = std::pow(TMath::NormQuantile((1+cdf)/2),2);
    }
    return map_chi2_quantiles;
  }

  //helper function to invert a matrix
  TMatrixD Conditional_Constraint::invert (TMatrixD cov) {
    int ncols = cov.GetNcols();
    Double_t cov_array[ncols*ncols];
    TMatrixDSym cov_sym(ncols);
    TMatrixDSym cov_inv(ncols);

    for (int i=0;i<ncols;i++) { for (int j=0;j<ncols;j++) { cov_array[i*ncols+j] = cov(i,j); } }
    cov_sym.SetMatrixArray(cov_array);
    TDecompChol cho(cov_sym);
    cho.Decompose();
    cho.Invert(cov_inv);
    //TMatrixDSymEigen cov_eigen(cov_sym);
    //TMatrixDSymEigen cov_inv_eigen(cov_inv);

    return (TMatrixD)cov_inv;
  }

  //helper function to compute the posterior prediction given a prior prediction (CV and cov) and measurement over a subset of those bins
  //data and mc_full are Nx1 matrices just to allow easy computation
  //x is the measurement space, y is the remaining space unconstrained by the measurement
  void Conditional_Constraint::apply_constraint(TMatrixD data, TMatrixD mc_full, TMatrixD cov, bool yfront, TMatrixD &mc_y_post, TMatrixD &mc_y_post_cov) {
  
    int nx = data.GetNrows();
    int ny = cov.GetNrows() - nx;
    TMatrixD mc_x(nx,1);
    TMatrixD mc_y(ny,1);
    TMatrixD cov_xx(nx,nx);
    TMatrixD cov_yy(ny,ny);
  
    if (yfront) {
      mc_x = mc_full.GetSub(ny, nx+ny-1, 0,  0);
      mc_y = mc_full.GetSub(0,  ny-1,    0,  0);
      cov_xx = cov.GetSub(  ny, nx+ny-1, ny, nx+ny-1);
      cov_yy = cov.GetSub(  0,  ny-1,    0,  ny-1);
    } else {
      mc_x = mc_full.GetSub(0,  nx-1,    0,  0);
      mc_y = mc_full.GetSub(nx, nx+ny-1, 0,  0);
      cov_xx = cov.GetSub(  0,  nx-1,    0,  nx-1);
      cov_yy = cov.GetSub(  nx, nx+ny-1, nx, nx+ny-1);
    }
  
    TMatrixD diff_x = TMatrixD(data, TMatrixD::kMinus, mc_x);
    TMatrixD cov_xy = TMatrixD(nx, ny);
    TMatrixD cov_yx = TMatrixD(ny, nx);
    for (int i=0;i<nx;i++) {
      for (int j=0;j<ny;j++) {
        int ix = i + ny *   yfront;
        int jy = j + nx * (!yfront);
        cov_xy(i,j) = cov(ix,jy);
        cov_yx(j,i) = cov(ix,jy);
      }
    }
    TMatrixD cov_xx_inv = invert(cov_xx);
  
    TMatrixD temp1 = TMatrixD(cov_yx, TMatrixD::kMult, cov_xx_inv);
    TMatrixD temp2 = TMatrixD(temp1,  TMatrixD::kMult, diff_x);
   TMatrixD temp3 = TMatrixD(temp1,  TMatrixD::kMult, cov_xy);
    mc_y_post      = TMatrixD(mc_y,   TMatrixD::kPlus,  temp2);
    mc_y_post_cov  = TMatrixD(cov_yy, TMatrixD::kMinus, temp3);
  }


  //helper function to compute the chi2 between two vectors given a covariance matrix
  double Conditional_Constraint::compute_chi2 (TMatrixD m1, TMatrixD m2, TMatrixD cov) {
    TMatrixD md  = m2 - m1;
    TMatrixD mdT = m2 - m1;
    mdT.T();
    TMatrixD cov_inv = invert(cov);
    TMatrixD mret = mdT * cov_inv * md;
    return mret(0,0); 
  }

