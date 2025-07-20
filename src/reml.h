#ifndef REML_H_
#define REML_H_

#include "mph.h"

// REML/MINQUE core functions
void run_reml_minque(
	const std::vector<std::map<std::string, std::pair<float, float>>>& indi2pheno_weight,
	const std::map<std::string, VectorXf>& indi2covar,
	const std::vector<std::string>& trait_names,
	std::vector<std::string>& covar_names,
	const std::string& grm_list_file,
	const std::string& binary_grm_file,
	const std::string& keep_file,
	const std::string& output_file,
	const std::vector<float>& h2,
	int minq_iter,
	float minq_tol,
	int n_rand_vec,
	int n_grms,
	const std::string& rdist,
	int seed,
	bool save_mem,
	bool verbose
);

// Helper functions
void project_xmat(const MatrixXf &HinvX, const LDLT<MatrixXf>& lltXtHinvX, const MatrixXf &X, Ref<MatrixXf> Y);
float dogleg(const VectorXf &grad, const MatrixXf &FI, const float delta, Ref<VectorXf> p, float &tau_nr, float &tau_sd);
void optim_quad(const float a, const float b, const float c, const float xlb, const float xub, float &xmin, float &xmax, float &ymin, float &ymax);

inline double getOneDfChisqPval (const double x) {
	if(x <= 0) return(1.0);
	return erfc(sqrt(x) * M_SQRT1_2);
}

#endif // REML_H_
