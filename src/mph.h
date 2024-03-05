#ifndef MPH_H_
#define MPH_H_

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif // EIGEN_USE_MKL_ALL

#include <cmath>
#include <cstdio>
#include <ctime>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifndef EIGEN_MKL_LIB
#define EIGEN_MKL_LIB
#ifndef EIGEN_LIB
#define EIGEN_LIB
#include "Eigen/Dense"
#endif // EIGEN_LIB
#include <mkl.h>
#endif // EIGEN_MKL_LIB

#include <omp.h>

#include "StrFunc.h"

using namespace Eigen;

typedef Matrix<signed char, Dynamic, Dynamic> MatrixXc;
typedef Matrix<signed char, Dynamic, 1> VectorXc;

inline bool file_check (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

// global variables
extern int num_threads;

// input data
std::map<std::string, std::pair<float, float> > read_phenotype_file(std::string phenotype_file, std::string trait_name, std::string error_variance_weight_name);
std::map<std::string, VectorXf> read_covariate_file(std::string covariate_file, std::vector<std::string>& covariate_names);
void read_marker_info_file(
	std::string marker_info_file, 
	std::string group_header, 
	std::string weight_name, 
	std::map<std::string, std::pair<std::string, float> >& marker2group_weight);
void get_subject_set(
	const std::string subject_set_file,
	const std::string plink_fam_file,
	std::vector<bool>& bindi,
	std::vector<std::string>& indi_keep );
void get_marker_set(
	const std::map<std::string, std::pair<std::string, float> >& marker2group_weight, 
	const std::string plink_bim_file,
	std::vector<bool>& bmarker,
	std::vector<double>& gvec );
void calc_hwep_geno012(const Ref<MatrixXc>& geno, const int midp, Ref<VectorXd> hwep);

// GRM
void calc_grm_by_subset(const char grm_type, const float min_maf, const float min_hwep, 
	Ref<MatrixXd> SS, const Ref<MatrixXc>& kmat, const Ref<VectorXd>& hwep, Ref<VectorXd> swt, 
	Ref<MatrixXd> grm, double &sumwt, int &post_qc_marker_num);
void write_grm_into_file(const std::string bin_file_prefix, const std::vector<std::string>& indi_keep, const float sum2pq, const Ref<MatrixXf> grm);
void write_dGRM_into_file(const std::string bin_file_prefix, const std::vector<std::string>& indi_keep, const float sum2pq, const Ref<MatrixXd> grm);
int get_grm_dim(const std::string binary_grm_file_prefix);
void get_grm_indi(const std::string indi_file, std::vector<std::string>& indi_keep);
void read_grm_from_file(const std::string binary_grm_file_prefix, const std::vector<std::string>& indi_keep, float& sum2pq, Ref<MatrixXf> grm, bool verbose=true);
void read_dGRM_from_file(const std::string binary_grm_file_prefix, const std::vector<std::string>& indi_keep, float& sum2pq, Ref<MatrixXd> grm);
void merge_grms(const std::string binary_grm_list_file, const std::vector<std::string>& indi_keep, float& out_sum2pq, Ref<MatrixXf> out_grm);
void deduct_grms(const std::string grm_list_file, const std::vector<std::string>& indi_keep, float& out_sum2pq, Ref<MatrixXf> out_grm);
void sum_grms(const std::string binary_grm_list_file, const std::vector<std::string>& indi_keep, float& gvar, Ref<MatrixXf> out_grm);
void make_core(const std::string grm_list_file, const std::vector<std::string>& indi_keep, const std::string output_file);
void make_fore(const std::string grm_list_file, const std::vector<std::string>& indi_keep, const std::string output_file);

#endif // MPH_H_
