#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompChol.h"
#include "TMatrixDSymEigen.h"
#include "src/conditional_constraint.h"
#include "src/conditional_constraint.cxx"

int main (int argc, char *argv[]) {

  //inputs
  std::string file_reg = "mc_universes/regularization_matrix.root";
  std::string path_MC = "mc_universes/";
  std::vector<std::string> v_fit_par_names = {};
  int nbins = 69;
  bool flag_reg = false;

  //Build class and fit
  Conditional_Constraint cc(nbins, flag_reg, file_reg, path_MC, v_fit_par_names);
  cc.build_cov();
  cc.fit();
  std::cout << "fit without regularization_matrix ==========" << std::endl;
  std::cout << "    nominal chi2/ndf = " << cc.chi2_nom     << " / " << nbins << std::endl;
  std::cout << "        fit chi2/ndf = " << cc.chi2_fit     << " / " << nbins << std::endl;
  TMatrixD m_data     = cc.m_data;
  TMatrixD m_pred_nom = cc.m_pred_nom;
  TMatrixD m_pred_fit = cc.m_pred_xs_fit;
  TMatrixD cov        = cc.cov;

  //limit the number of fit parameters and refit
  v_fit_par_names = {"MaCCQE_t_MC","RPA_CCQE_t_MC","NormCCMEC_t_MC","XSecShape_CCMEC_t_MC"};
  cc.set_fit_pars(v_fit_par_names);
  cc.build_cov();
  cc.fit();
  std::cout << "limited fit chi2/ndf = " << cc.chi2_lim_fit << " / " << nbins << std::endl;
  TMatrixD m_pred_lim_fit = cc.m_pred_xs_lim_trim_fit;

  //Apply quantile mapping and refit
  cc.scale_cov_QM();
  cc.build_cov();
  cc.fit();
  std::cout << "    nominal QM chi2/ndf = " << cc.chi2_nom     << " / " << nbins << std::endl;
  std::cout << "        fit QM chi2/ndf = " << cc.chi2_fit     << " / " << nbins << std::endl;
  std::cout << "limited fit QM chi2/ndf = " << cc.chi2_lim_fit << " / " << nbins << std::endl;
  TMatrixD m_pred_fit_QM     = cc.m_pred_xs_fit;
  TMatrixD m_pred_lim_fit_QM = cc.m_pred_xs_lim_trim_fit;
  TMatrixD cov_QM            = cc.cov;

  //Apply regularization matrix, requires reconfiguring
  flag_reg = true;
  cc.configure(nbins, flag_reg, file_reg, path_MC, v_fit_par_names);
  cc.build_cov();
  cc.fit();
  std::cout << "fit with regularization_matrix ==========" << std::endl;
  std::cout << "    nominal chi2/ndf = " << cc.chi2_nom     << " / " << nbins << std::endl;
  std::cout << "        fit chi2/ndf = " << cc.chi2_fit     << " / " << nbins << std::endl;
  std::cout << "limited fit chi2/ndf = " << cc.chi2_lim_fit << " / " << nbins << std::endl;
  TMatrixD m_data_reg         = cc.m_data;
  TMatrixD m_pred_reg_nom     = cc.m_pred_nom;
  TMatrixD m_pred_reg_fit     = cc.m_pred_xs_fit;
  TMatrixD m_pred_reg_lim_fit = cc.m_pred_xs_lim_trim_fit;
  TMatrixD cov_reg            = cc.cov;

  //Prepare plots
  gStyle->SetOptStat(0);
  TH1D* h_data             = new TH1D("h_data",            "h_data",            nbins,0,nbins);
  TH1D* h_data_QM          = new TH1D("h_data_QM",         "h_data_QM",         nbins,0,nbins);
  TH1D* h_data_reg         = new TH1D("h_data_reg",        "h_data_reg",        nbins,0,nbins);
  TH1D* h_pred_nom         = new TH1D("h_pred_nom",        "h_pred_nom",        nbins,0,nbins);
  TH1D* h_pred_fit         = new TH1D("h_pred_fit",        "h_pred_fit",        nbins,0,nbins);
  TH1D* h_pred_lim_fit     = new TH1D("h_pred_lim_fit",    "h_pred_lim_fit",    nbins,0,nbins);
  TH1D* h_pred_fit_QM      = new TH1D("h_pred_fit_QM",     "h_pred_fit_QM",     nbins,0,nbins);
  TH1D* h_pred_lim_fit_QM  = new TH1D("h_pred_lim_fit_QM", "h_pred_lim_fit_QM", nbins,0,nbins);
  TH1D* h_pred_reg_nom     = new TH1D("h_pred_reg_nom",    "h_pred_reg_nom",    nbins,0,nbins);
  TH1D* h_pred_reg_fit     = new TH1D("h_pred_reg_fit",    "h_pred_reg_fit",    nbins,0,nbins);
  TH1D* h_pred_reg_lim_fit = new TH1D("h_pred_reg_lim_fit","h_pred_reg_lim_fit",nbins,0,nbins);
  for (int i=0;i<nbins;i++) {
    h_data->SetBinContent(i+1,            m_data(i,0));
    h_data->SetBinError(i+1,              std::sqrt(cov(i,i)));
    h_data_QM->SetBinContent(i+1,         m_data(i,0));
    h_data_QM->SetBinError(i+1,           std::sqrt(cov_QM(i,i)));
    h_data_reg->SetBinContent(i+1,        m_data_reg(i,0));
    h_data_reg->SetBinError(i+1,          std::sqrt(cov_reg(i,i)));
    h_pred_nom->SetBinContent(i+1,        m_pred_nom(i,0));
    h_pred_fit->SetBinContent(i+1,        m_pred_fit(i,0));
    h_pred_lim_fit->SetBinContent(i+1,    m_pred_lim_fit(i,0));
    h_pred_fit_QM->SetBinContent(i+1,     m_pred_fit_QM(i,0));
    h_pred_lim_fit_QM->SetBinContent(i+1, m_pred_lim_fit_QM(i,0));
    h_pred_reg_nom->SetBinContent(i+1,    m_pred_reg_nom(i,0));
    h_pred_reg_fit->SetBinContent(i+1,    m_pred_reg_fit(i,0));
    h_pred_reg_lim_fit->SetBinContent(i+1,m_pred_reg_lim_fit(i,0));
  }

  //Data vs nominal, fit, and limited fit models
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->cd();
  h_data->SetLineWidth(2);
  h_pred_nom->SetLineWidth(2);
  h_pred_fit->SetLineWidth(2);
  h_pred_lim_fit->SetLineWidth(2);
  h_data->SetLineColor(kBlack);
  h_pred_nom->SetLineColor(kBlue);
  h_pred_fit->SetLineColor(kRed);
  h_pred_lim_fit->SetLineColor(kGreen);
  h_data->Draw();
  h_pred_nom->Draw("same");
  h_pred_fit->Draw("same");
  h_pred_lim_fit->Draw("same");
  c1->SaveAs("figs/nominal_fit.pdf");

  //Apply quantile mapping
  TCanvas* c2 = new TCanvas("c2","c2");
  c2->cd();
  h_data_QM->SetLineWidth(2);
  h_pred_nom->SetLineWidth(2);
  h_pred_fit_QM->SetLineWidth(2);
  h_pred_lim_fit_QM->SetLineWidth(2);
  h_data_QM->SetLineColor(kBlack);
  h_pred_nom->SetLineColor(kBlue);
  h_pred_fit_QM->SetLineColor(kRed);
  h_pred_lim_fit_QM->SetLineColor(kGreen);
  h_data_QM->Draw();
  h_pred_nom->Draw("same");
  h_pred_fit_QM->Draw("same");
  h_pred_lim_fit_QM->Draw("same");
  c2->SaveAs("figs/quantile_mapped_fit.pdf");

  //Apply regularization matrix
  TCanvas* c3 = new TCanvas("c3","c3");
  c3->cd();
  h_data_reg->SetLineWidth(2);
  h_pred_reg_nom->SetLineWidth(2);
  h_pred_reg_fit->SetLineWidth(2);
  h_pred_reg_lim_fit->SetLineWidth(2);
  h_data_reg->SetLineColor(kBlack);
  h_pred_reg_nom->SetLineColor(kBlue);
  h_pred_reg_fit->SetLineColor(kRed);
  h_pred_reg_lim_fit->SetLineColor(kGreen);
  h_data_reg->Draw();
  h_pred_reg_nom->Draw("same");
  h_pred_reg_fit->Draw("same");
  h_pred_reg_lim_fit->Draw("same");
  c3->SaveAs("figs/regularized_fit.pdf");

}



