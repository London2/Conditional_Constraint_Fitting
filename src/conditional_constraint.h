#pragma once
#include <cmath>

class Conditional_Constraint{
  public:
    int nbins;
    double chi2_nom, chi2_fit, chi2_lim_fit;
    mutable TMatrixD m_reg, cov, cov_data, cov_pred, cov_decomp, cov_decomp_QM, cov_xs_fit, cov_par_fit;
    mutable TMatrixD m_data, m_pred_nom, m_pred_xs_fit, m_pred_xs_lim_trim_fit, m_pred_par_fit, m_pred_par_lim_fit, m_chi2_decomp;
    std::vector<std::string> v_par_names, v_fit_par_names;

    Conditional_Constraint (int nbins, bool flag_reg, std::string file_reg, std::string path_MC, std::vector<std::string> v_fit_par_names);
    void set_fit_pars (std::vector<std::string> v_fit_par_names);
    void configure (int nbins, bool flag_reg, std::string file_reg, std::string path_MC, std::vector<std::string> v_fit_par_names);
    void build_cov ();
    void scale_cov_QM ();
    void fit ();

  protected:
    bool flag_reg;
    int nfiles, npars, nbins_ext, nfit, nthrows;
    std::vector<std::string> filelist_MC;
    std::vector<std::vector<double>> v_pred, v_par_throw_vals;
    std::vector<double> v_par_nom;
    std::map<std::string,double> map_par_nom_val;
    std::map<std::string,std::vector<double>> map_par_throw_vals;
    mutable TVectorD cov_eigvals;
    mutable TMatrixD cov_ext, cov_data_ext, cov_pred_ext, cov_pred_xs_lim_fit, cov_pred_par_lim_fit, cov_eigvecs, cov_eigvecs_T;
    mutable TMatrixD m_data_lim, m_pred_nom_ext, m_pred_xs_lim_fit, m_diff, m_diff_decomp;

    void read_reg (std::string file_reg);
    void read_filelist (std::string path_MC);
    void resize_nbins ();
    void read_universes ();
    void setup_pars ();
    void setup_joint_measurement();
    void diagonalize_cov ();
    void apply_constraint(TMatrixD data, TMatrixD mc_full, TMatrixD cov, bool yfront, TMatrixD &mc_y_post, TMatrixD &mc_y_post_cov);
    TMatrixD invert (TMatrixD cov);
    std::map<double,double> set_chi2_quantile_map (TMatrixD m_chi2_meas);
    double compute_chi2 (TMatrixD m1, TMatrixD m2, TMatrixD cov);
};
