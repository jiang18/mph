
#include "mph.h"

void write_snp_into_files(
	const std::string output_file_prefix,
	std::map<std::string, double>& marker2maf,
	std::map<std::string, std::pair<std::string, int> >& marker2pos,
	std::map<std::string, std::pair<std::string, std::string> >& marker2alleles,
	std::map<std::string, std::pair<std::string, double> >& marker2group_weight,
	const Ref<VectorXi> ncols,
	const std::vector<std::string>& marker_keep,
	const std::vector<VectorXd>& uN) 
{
	unsigned i=0, j=0, k=0;
	std::ofstream ofs;
	// snp effects
	ofs.open(output_file_prefix + ".snp.csv");
	ofs<<"SNP,Chr,Pos,Allele1,Allele2,MAF,Group,Weight,Effect"<<std::endl;
	for(i=0; i<ncols.size(); ++i) {
		for(j=0; j<ncols[i]; ++j) {
			std::string this_marker = marker_keep[k];
			ofs<<this_marker<<','<<marker2pos[this_marker].first<<','<<marker2pos[this_marker].second<<',';
			ofs<<marker2alleles[this_marker].first<<','<<marker2alleles[this_marker].second<<','<<marker2maf[this_marker]<<',';
			ofs<<marker2group_weight[this_marker].first<<','<<marker2group_weight[this_marker].second<<',';
			ofs<<uN[i][j]<<std::endl;
			++k;
		}
	}
	ofs.close();
}


void write_set_into_files(
	const std::string output_file_prefix,
	const Ref<VectorXi> ncols,
	const std::vector<std::string>& group_keep,
	const Ref<VectorXd> varu ) 
{
	unsigned i=0;
	std::ofstream ofs;
	// group parameter estimates
	ofs.open(output_file_prefix + ".set.csv");
	ofs<<"Set,Size,Var"<<std::endl;
	for(i=0; i<group_keep.size(); ++i) {
		ofs<<group_keep[i]<<','<<ncols[i]<<','<<varu[i]<<std::endl;
	}
	ofs.close();
}


void write_covar_into_files(
	const std::string output_file_prefix,
	const std::vector<std::string>& covar_names,
	const Ref<VectorXd> bN ) 
{
	unsigned i=0;
	std::ofstream ofs;
	// covariates
	ofs.open(output_file_prefix + ".covar.csv");
	for(i=0; i<covar_names.size(); ++i) 
		ofs<<covar_names[i]<<','<<bN[i]<<std::endl;
	ofs.close();
}


void write_model_into_files(
	const std::string output_file_prefix,
	const int indi_keep_num,
	const Ref<VectorXi> ncols,
	const double vare,
	bool hyper_opt_flag,
	const Ref<VectorXd> hyper ) 
{
	std::ofstream ofs;
	// free energy and error parameter estimates
	ofs.open(output_file_prefix + ".model.csv");
	ofs<<"pop_size,"<<indi_keep_num<<std::endl;
	ofs<<"num_markers,"<<ncols.sum()<<std::endl;
	ofs<<"err_var,"<<vare<<std::endl;
	ofs<<"hyper_opt,"<<hyper_opt_flag<<std::endl;
	ofs<<"hyper[0],"<<hyper[0]<<std::endl;
	ofs.close();
}


void write_output_into_files(
	const std::string output_file_prefix,
	const int indi_keep_num,
	std::map<std::string, double>& marker2maf,
	std::map<std::string, std::pair<std::string, int> >& marker2pos,
	std::map<std::string, std::pair<std::string, std::string> >& marker2alleles,
	std::map<std::string, std::pair<std::string, double> >& marker2group_weight,
	const Ref<VectorXi> ncols,
	const std::vector<std::string>& group_keep,
	const std::vector<std::string>& marker_keep,
	const std::vector<std::string>& covar_names,
	const std::vector<VectorXd>& uN,
	const Ref<VectorXd> kshapeN,
	const Ref<VectorXd> kscaleN,
	const Ref<VectorXd> bN,
	const double& cN,
	const double& dN,
	const double& fe) 
{
	unsigned i=0, j=0, k=0;
	std::ofstream ofs;
	// snp effects
	ofs.open(output_file_prefix + ".snp.csv");
	ofs<<"SNP,Chr,Pos,Allele1,Allele2,MAF,Group,Weight,Effect"<<std::endl;
	for(i=0; i<group_keep.size(); ++i) {
		for(j=0; j<ncols[i]; ++j) {
			std::string this_marker = marker_keep[k];
			ofs<<this_marker<<','<<marker2pos[this_marker].first<<','<<marker2pos[this_marker].second<<',';
			ofs<<marker2alleles[this_marker].first<<','<<marker2alleles[this_marker].second<<','<<marker2maf[this_marker]<<',';
			ofs<<marker2group_weight[this_marker].first<<','<<marker2group_weight[this_marker].second<<',';
			ofs<<uN[i][j]<<std::endl;
			++k;
		}
	}
	ofs.close();
	// group parameter estimates
	ofs.open(output_file_prefix + ".set.csv");
	ofs<<"Set,Size,Inv_Var,Var_Posterior_Shape,Var_Posterior_Scale,Var"<<std::endl;
	for(i=0; i<group_keep.size(); ++i) {
		ofs<<group_keep[i]<<','<<ncols[i]<<','<<kshapeN[i]/kscaleN[i]<<','<<kshapeN[i]<<','<<kscaleN[i];
		if(kshapeN[i] > 1.) ofs<<','<<kscaleN[i]/(kshapeN[i]-1.);
		ofs<<std::endl;
	}
	ofs.close();
	// covariates
	ofs.open(output_file_prefix + ".covar.csv");
	for(i=0; i<covar_names.size(); ++i) 
		ofs<<covar_names[i]<<','<<bN[i]<<std::endl;
	ofs.close();
	// free energy and error parameter estimates
	ofs.open(output_file_prefix + ".f.csv");
	ofs<<"pop_size,"<<indi_keep_num<<std::endl;
	ofs<<"num_markers,"<<ncols.sum()<<std::endl;
	ofs<<"free_energy,"<<fe<<std::endl;
	ofs<<"inv_err_var,"<<cN/dN<<std::endl;
	ofs<<"err_var_posterior_shape,"<<cN<<std::endl;
	ofs<<"err_var_posterior_scale,"<<dN<<std::endl;
	ofs<<"err_var,"<<dN/(cN-1)<<std::endl;
	ofs.close();
}

void write_cv_into_file (std::map<int, VectorXd>& pred_eval, const std::string output_file_prefix)
{
	unsigned i=0;
	std::map<int, VectorXd>::iterator it;
	std::ofstream ofs;
	ofs.open(output_file_prefix + ".cv.csv");
	unsigned num_rep=pred_eval.begin()->second.size();
	VectorXd reps (num_rep);
	ofs<<"rep";
	for(it=pred_eval.begin(); it!=pred_eval.end(); ++it) ofs<<","<<it->first;
	ofs<<std::endl;
	for (i=0; i<num_rep; ++i) {
		ofs<<i+1;
		for(it=pred_eval.begin(); it!=pred_eval.end(); ++it) ofs<<","<<it->second[i];
		ofs<<std::endl;
	}
	ofs<<"mean";
	for(it=pred_eval.begin(); it!=pred_eval.end(); ++it) ofs<<","<<it->second.mean();
	ofs<<std::endl;
	ofs<<"stderr";
	for(it=pred_eval.begin(); it!=pred_eval.end(); ++it) {
		reps = it->second;
		reps.array() -= reps.mean();
		double stderr = reps.norm()/std::sqrt(num_rep-1)/std::sqrt(num_rep);
		ofs<<","<<stderr;
	}
	ofs<<std::endl;
}
