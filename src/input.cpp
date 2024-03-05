#include <set>
#include "mph.h"
#include "snphwe.h"

const std::string default_group ("NULL");
const float default_marker_effvar_weight = 1.0;
const float default_error_variance_weight = 1.0;
const char delimiter = ',';

void calc_hwep_geno012(const Ref<MatrixXc>& geno, const int midp, Ref<VectorXd> hwep)
{
	#pragma omp parallel for
	for(int i=0; i<geno.cols(); ++i) {
		int n0 = (geno.col(i).array() == 0).count();
		int n1 = (geno.col(i).array() == 1).count();
		int n2 = (geno.col(i).array() == 2).count();

		hwep(i) = SNPHWE2(n1, n0, n2, midp);
	}
}

// Remember to take advantage of Return Value Optimization in the main function.
std::map<std::string, std::pair<float, float> > read_phenotype_file(std::string phenotype_file, std::string trait_name, std::string error_variance_weight_name)
{
	std::map<std::string, std::pair<float, float> > indi2pheno_weight;
	std::ifstream ifs;
	std::string line;
	
	ifs.open(phenotype_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<phenotype_file<<std::endl;
		exit(1);
	}
	
	std::getline(ifs, line);
	if(line[line.size()-1] == '\r') line.erase(line.size()-1);
	std::vector<std::string> header = StrFunc::split(line, delimiter);
	
	int trait_col = StrFunc::find_token_in_header(trait_name, header, phenotype_file);
	int error_variance_weight_col = StrFunc::find_token_in_header(error_variance_weight_name, header, phenotype_file);
	// Beware of the situation where the last column is missing.
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::vector<std::string> cols = StrFunc::split(line, delimiter);
		if(!(header.size() == cols.size() || (line.back() == delimiter && header.size() == cols.size()+1))) {
			std::cout<<"\nError: "<<line<<" in "<<phenotype_file<<std::endl;
			exit(1);
		}
		if(error_variance_weight_col == na_col && cols.size() > trait_col && !cols[trait_col].empty())
			indi2pheno_weight[cols[0]] = std::make_pair(std::stof(cols[trait_col]), default_error_variance_weight);
		else if(cols.size() > error_variance_weight_col && !cols[error_variance_weight_col].empty() && cols.size() > trait_col && !cols[trait_col].empty())
			indi2pheno_weight[cols[0]] = std::make_pair(std::stof(cols[trait_col]), std::stof(cols[error_variance_weight_col]));
	}
	ifs.close();
	
	return indi2pheno_weight;
}

std::map<std::string, VectorXf> read_covariate_file(std::string covariate_file, std::vector<std::string>& covariate_names)
{
	std::map<std::string, VectorXf > indi2covar;
	int i;
	std::ifstream ifs;
	std::string line;
	
	ifs.open(covariate_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<covariate_file<<std::endl;
		exit(1);
	}
	
	std::map<std::string, int> name2col; 
	std::getline(ifs, line);
	if(line[line.size()-1] == '\r') line.erase(line.size()-1);
	std::vector<std::string> header = StrFunc::split(line, delimiter);
	std::cout<<"Covariates in covariate file: "<<std::endl;
	for (i=1; i<header.size(); ++i) {
		std::cout<<" "<<header[i];
		name2col[header[i]] = i;
	}
	std::cout<<std::endl;
	std::vector<int> covariate_cols;
	std::cout<<"Covariates in command line option: "<<std::endl;
	for (i=0; i<covariate_names.size(); ++i) {
		std::cout<<" "<<covariate_names[i];
		if(name2col.find(covariate_names[i]) == name2col.end()) {
			std::cout<<"\n"<<covariate_names[i]<<" not in covariate file"<<std::endl;
			exit(1);
		} else covariate_cols.push_back(name2col[covariate_names[i]]);
	}
	std::cout<<std::endl;
	VectorXf covars (covariate_cols.size());
	// Beware of the situation where the last column is missing.
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::vector<std::string> cols = StrFunc::split(line, delimiter);
		if(!(header.size() == cols.size() || (line.back() == delimiter && header.size() == cols.size()+1))) {
			std::cout<<"\nError: "<<line<<" in "<<covariate_file<<std::endl;
			exit(1);
		}
		if(covariate_cols.empty()) indi2covar[cols[0]] = VectorXf::Ones(1);
		else {
			for (i=0; i<covariate_cols.size(); ++i) {
				if(cols.size() > covariate_cols[i] && !cols[covariate_cols[i]].empty()) covars[i] = std::stof(cols[covariate_cols[i]]);
				else goto nextline;
			}
			indi2covar[cols[0]] = covars;
		}
		nextline:;
	}
	ifs.close();
	if(covariate_cols.empty()) covariate_names.push_back("intercept");
	return indi2covar;
}

void read_marker_info_file(
	std::string marker_info_file, 
	std::string group_header, 
	std::string weight_name, 
	std::map<std::string, std::pair<std::string, float> >& marker2group_weight)
{
	marker2group_weight.clear();
	
	std::ifstream ifs;
	std::string line;
	
	ifs.open(marker_info_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<marker_info_file<<std::endl;
		exit(1);
	}
	
	std::getline(ifs, line);
	if(line[line.size()-1] == '\r') line.erase(line.size()-1);
	std::vector<std::string> header = StrFunc::split(line, delimiter);
	
	int group_col = StrFunc::find_token_in_header(group_header, header, marker_info_file);
	if(!group_header.empty() && group_col == na_col) {
		std::cout<<"Marker group column "<<group_header<<" not found"<<std::endl;
		exit(1);
	}
	int weight_col = StrFunc::find_token_in_header(weight_name, header, marker_info_file);
	if(!weight_name.empty() && weight_col == na_col) {
		std::cout<<"Marker weight column "<<weight_name<<" not found"<<std::endl;
		exit(1);
	}
	
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::vector<std::string> cols = StrFunc::split(line, delimiter);
		if(!(header.size() == cols.size() || (line.back() == delimiter && header.size() == cols.size()+1))) {
			std::cout<<"\nError: "<<line<<" in "<<marker_info_file<<std::endl;
			exit(1);
		}
		if(weight_col == na_col) {
				marker2group_weight[cols[0]] = std::make_pair(default_group, default_marker_effvar_weight);
		} else if(cols.size() > weight_col && !cols[weight_col].empty()) {
			marker2group_weight[cols[0]] = std::make_pair(default_group, std::stof(cols[weight_col]));
		}
	}
	ifs.close();
}

void get_subject_set(
	const std::string subject_set_file,
	const std::string plink_fam_file,
	std::vector<bool>& bindi,
	std::vector<std::string>& indi_keep )
{
	bindi.clear();
	indi_keep.clear();
	
	std::ifstream ifs;
	std::string line;
	std::set<std::string> keep;
	
	if(subject_set_file != "") {
		ifs.open(subject_set_file);
		if(!ifs.is_open()) {
			std::cout<<"Cannot open "<<subject_set_file<<std::endl;
			exit(1);
		}
		while( std::getline(ifs, line) ){
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			if(line.size() < 1) continue;
			keep.insert(line);
		}
		ifs.clear();
		ifs.close();
	}
	
	ifs.open(plink_fam_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<plink_fam_file<<std::endl;
		exit(1);
	}
	std::stringstream ss;
	std::string iid;
	while( std::getline(ifs, line) ){
	// line content in plink fam: FID IID IID-father IID-mother Sex Phenotype
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		if(line.size() < 1) continue;
		ss.str(line);
		ss >> iid >> iid;
		
		if(keep.empty() || (keep.find(iid) != keep.end())) {
			bindi.push_back(true);
			indi_keep.push_back(iid);
		} else bindi.push_back(false);
	}
	ifs.clear();
	ifs.close();
}

void get_marker_set(
	const std::map<std::string, std::pair<std::string, float> >& marker2group_weight, 
	const std::string plink_bim_file,
	std::vector<bool>& bmarker,
	std::vector<double>& gvec )
{
	bmarker.clear();
	gvec.clear();
	
	std::ifstream ifs;
	ifs.open(plink_bim_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<plink_bim_file<<std::endl;
		exit(1);
	}
	std::string line;
	std::stringstream ss;
	std::string snp;
	while( std::getline(ifs, line) ){
	// line content in plink bim: Chr SNP CM BP A1 A2
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		if(line.size() < 1) continue;
		ss.str(line);
		ss >> snp >> snp;
		
		if(marker2group_weight.find(snp) != marker2group_weight.end()) {
			bmarker.push_back(true);
			gvec.push_back(marker2group_weight.at(snp).second);
		} else bmarker.push_back(false);
	}
	ifs.clear();
	ifs.close();
}

