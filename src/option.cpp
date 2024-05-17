#include <bitset>
#include <getopt.h>
#include <random>
#include <thread>
#include "mph.h"

int num_threads;

void print_help();
void option(int option_num, char **option_str);
void project_xmat(const MatrixXf &HinvX,const LDLT<MatrixXf>& lltXtHinvX,const MatrixXf &X,Ref<MatrixXf> Y);
float dogleg(const VectorXf &grad, const MatrixXf &FI, const float delta, Ref<VectorXf> p, float &tau_nr, float &tau_sd);
inline double getOneDfChisqPval (const double x) {
	if(x <= 0) return(1.0);
	return erfc(sqrt(x) * M_SQRT1_2);
}

int main(int argc, char **argv)
{   
	std::cout<<"***************************************************************************"<<std::endl;
	std::cout<<"* MPH by Jicai Jiang"<<std::endl;
	std::cout<<"* MINQUE for Partitioning Heritability"<<std::endl;
	std::cout<<"* Version 0.53.2 (May 16, 2024)"<<std::endl;
	std::cout<<"* (C) 2021-present, Jicai Jiang, NC State University"<<std::endl;
	std::cout<<"***************************************************************************"<<std::endl;
    
	time_t start=std::time(nullptr);
	std::cout<<"Analysis started: "<<std::ctime(&start)<<std::endl;
	try{ option(argc, argv); }
	catch(const std::string &err_msg){ std::cerr<<"\n"<<err_msg<<std::endl; }
	catch(const char *err_msg){ std::cerr<<"\n"<<err_msg<<std::endl; }
	time_t end=std::time(nullptr);
	std::cout<<"\nAnalysis finished: "<<std::ctime(&end);
	long time_used=end-start;
	std::cout<<"Computational time: "<<time_used/3600<<":"<<(time_used%3600)/60<<":"<<time_used%60<<std::endl;

	return 0;
}

void option(int option_num, char **option_str) {
	if(option_num < 3) {
		print_help();
		exit(1);
	}
	// variables
	bool verbose = false;
	bool minq = false;
	bool reml = false;
	bool simu_pheno = false;
	bool pred_indi = false;
	bool make_grm = false;
	bool merge = false;
	bool deduct = false;
	bool zero_grm = false;
	bool is_core_set = false;
	bool is_fore_set = false;
	bool save_mem = false;
	
	char grm_type = 'A';
	int minq_iter = 20;
	float minq_tol = 0.01;
	int n_rand_vec = 100;
	int n_pheno = 100;
	int seed = 0;
	std::string rdist = "Rademacher";
	int n_grms = 0;
	float zero_cutoff = 0;
	std::string subject_set_file="";
	std::string phen_file="", covar_file="", str_traits="", str_error_variance_weights="", str_covar_names="";
	std::vector<std::string> trait_names;
	std::vector<std::string> error_variance_weight_names;
	std::vector<std::string> covar_names;
	std::string marker_info_file="", marker_group="", marker_effvar_weight_name="";
	std::string binary_genotype_file="";
	float min_maf=0;
	std::string output_file="";
	num_threads=1;

	int num_traits = 1;
	
	std::string str_h2 = "0.5";
	std::vector<float> h2 = {0.5};
	std::string binary_grm_file = "";
	std::string grm_list_file = "";
	std::string mq_file = "";
	float min_hwep = 0;
	const int hwe_midp = 0;
	
	// parse command-line arguments
	static struct option long_options[] =
	{
		{"verbose",  no_argument, 0, 'V'},
		{"minque",  no_argument, 0, 'M'},
		{"save_memory",  no_argument, 0, 'A'},
		{"reml",  no_argument, 0, 'R'},
		{"simulate",  no_argument, 0, 'S'},
		{"predict",  no_argument, 0, 'P'},
		{"make_grm",  no_argument, 0, 'G'},
		{"merge_grms",  no_argument, 0, 'H'},
		{"deduct_grms",  no_argument, 0, 'I'},
		{"make_core",  no_argument, 0, 'C'},
		{"make_fore",  no_argument, 0, 'F'},
		{"dominance",  no_argument, 0, 'D'},
		{"num_phenotypes",  required_argument, 0, 'x'},
		{"num_iterations",  required_argument, 0, 'p'},
		{"tolerance",  required_argument, 0, 'q'},
		{"num_random_vectors",  required_argument, 0, 'a'},
		{"seed",  required_argument, 0, 'k'},
		{"distribution",  required_argument, 0, 's'},
		{"num_grms",  required_argument, 0, 'd'},
		{"zero_grm",  required_argument, 0, 'u'},
		{"heritability",  required_argument, 0, 'h'},
		{"subject_set",  required_argument, 0, 'i'},
		{"binary_grm_file",  required_argument, 0, '2'},
		{"min_hwe_pval",  required_argument, 0, '4'},
		{"grm_list",  required_argument, 0, '6'},
		{"mq_file",  required_argument, 0, '8'},
		{"min_maf",  required_argument, 0, 'j'},
		{"phenotype_file",  required_argument, 0, 'f'},
		{"trait_names",  required_argument, 0, 't'},
		{"error_weight_names",  required_argument, 0, 'e'},
		{"covariate_file",  required_argument, 0, 'c'},
		{"covariate_names",  required_argument, 0, 'v'},
		{"snp_info_file",  required_argument, 0, 'm'},
		{"snp_set_name",  required_argument, 0, 'r'},
		{"snp_weight_name",  required_argument, 0, 'w'},
		{"binary_genotype_file",  required_argument, 0, 'b'},
		{"num_threads",    required_argument, 0, 'n'},
		{"output_file",    required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};
	int option_index = 0;
	int opt;
	while ((opt = getopt_long(option_num, option_str, "VMARSPGHICFDx:p:q:a:k:s:d:u:h:i:2:4:6:8:j:f:t:e:c:v:m:r:w:b:n:o:",long_options, &option_index)) != -1)
	{
		switch (opt) {
			case 'V':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				verbose = true;
				break;
			case 'M':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				minq = true;
				break;
			case 'A':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				save_mem = true;
				break;
			case 'R':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				reml = true;
				break;
			case 'S':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				simu_pheno = true;
				break;
			case 'P':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				pred_indi = true;
				break;
			case 'G':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				make_grm = true;
				break;
			case 'H':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				merge = true;
				break;
			case 'I':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				deduct = true;
				break;
			case 'C':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				is_core_set = true;
				break;
			case 'F':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				is_fore_set = true;
				break;
			case 'D':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				grm_type = 'D';
				break;
			case 'x':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				n_pheno = std::stoi(optarg);
				break;
			case 'p':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				minq_iter = std::stoi(optarg);
				break;
			case 'q':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				minq_tol = std::stof(optarg);
				break;
			case 'a':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				n_rand_vec = std::stoi(optarg);
				break;
			case 'k':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				seed = std::stoi(optarg);
				break;
			case 's':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				rdist = optarg;
				break;
			case 'd':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				n_grms = std::stoi(optarg);
				break;
			case 'u':
			{
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				zero_grm = true;
				zero_cutoff = std::stof(optarg);
				break;
			}
			case 'h':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				str_h2 = optarg;
				break;
			case 'i':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				subject_set_file = optarg;
				break;
			case '2':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				binary_grm_file = optarg;
				break;
			case '4':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				min_hwep = std::stof(optarg);
				break;
			case '6':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				grm_list_file = optarg;
				break;
			case '8':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				mq_file = optarg;
				break;
			case 'j':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				min_maf = std::stof(optarg);
				break;
			case 'f':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				phen_file = optarg;
				break;
			case 't':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				str_traits = optarg;
				break;
			case 'e':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				str_error_variance_weights = optarg;
				break;
			case 'c':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				covar_file = optarg;
				break;
			case 'v':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				str_covar_names = optarg;
				break;
			case 'm':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				marker_info_file = optarg;
				break;
			case 'r':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				marker_group = optarg;
				break;
			case 'w':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				marker_effvar_weight_name = optarg;
				break;
			case 'b':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				binary_genotype_file = optarg;
				break;
			case 'n':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				num_threads = std::stoi(optarg);
				break;
			case 'o':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				output_file = optarg;
				break;
			default:
				print_help();
				exit(1);
		}
	}
	
	// check options and associated options
	if(zero_cutoff < -1 || zero_cutoff > 1) throw("\nError: --zero_grm must be within [-1, 1].\n");
	if(min_maf < 0 || min_maf >= 0.5) throw("\nError: min MAF is out of [0, 0.5).\n");
	if(min_hwep < 0 || min_hwep >= 1) throw("\nError: min HWE P value is out of [0, 1).\n");
	if(num_threads > 100 || num_threads < 1) throw("\nError: number of threads is out of [1,100].\n");
	if(output_file.empty()) throw("\nError: output file is required for analysis.\n");
	if(reml || minq) {
		if(phen_file.empty()) throw("\nError: --phenotype is required for --reml/--minque.\n");
		if(str_traits.empty()) throw("\nError: --trait is required for --reml/--minque.\n");
		trait_names = StrFunc::split(str_traits,',');
		if(!str_covar_names.compare("all") && covar_file.empty()) 
			throw("\nError: --covariate_file is required when --covariate_names all is given.\n");
		if(covar_file.empty()) covar_file = phen_file;
		if(!str_covar_names.compare("all")) 
			covar_names = StrFunc::get_token_names(covar_file);
		else if(!str_covar_names.empty())
			covar_names = StrFunc::split(str_covar_names,',');

		num_traits = trait_names.size();
		for(int i=0; i<num_traits-1; ++i) 
			for(int j=i+1; j<num_traits; ++j) 
				if(trait_names[i] == trait_names[j]) throw("\nError: trait names must be unique in --trait.\n");
	}
	std::vector<std::string> str_vec = StrFunc::split(str_h2,',');
	if(str_vec.size() != 1 && str_vec.size() != num_traits) 
		throw("\nError: the number of values in --heritability must be 1 or correspond to --trait.\n");
	h2.resize(str_vec.size());
	for(int i=0; i<h2.size(); ++i)  h2[i] = std::stof(str_vec[i]);
	for(int i=h2.size(); i<num_traits; ++i)  h2[i] = h2[0];
	for(int i=0; i<num_traits; ++i)
		if(h2[i] <= 0 || h2[i] > 1) throw("\nError: --heritability is out of (0, 1].\n");

	if(!str_error_variance_weights.empty()) {
		error_variance_weight_names = StrFunc::split(str_error_variance_weights,',');
		if(error_variance_weight_names.size() != num_traits)
			throw("\nError: the number of values in --error_weight_names must correspond to --trait if specified.\n");
	} else error_variance_weight_names.resize(num_traits, "");

	if(reml || minq) {
		if(grm_list_file == "") throw("\nError: --grm_list must be set for --reml/--minque.\n");
		if(n_grms <0) throw("\nError: --num_grms must be >=0.\n");
		if(minq_tol < 1e-4) minq_tol = 1e-4;
		if(minq_iter < 1) throw("\nError: --num_iterations must be >0.\n");
		if(rdist != "Gaussian" && rdist != "Rademacher")
			throw("\nError: --distribution must be either Gaussian or Rademacher.\n");
	}
	if(make_grm && marker_info_file == "")
		throw("\nError: --snp_info_file must be set for --make_grm.\n");
	if(grm_list_file.empty()) {
		if(merge) throw("\nError: --grm_list is required for --merge_grms.\n");
		if(deduct) throw("\nError: --grm_list is required for --deduct_grms.\n");
		if(is_core_set) throw("\nError: --grm_list is required for --make_core.\n");
		if(is_fore_set) throw("\nError: --grm_list is required for --make_fore.\n");
		if(simu_pheno) throw("\nError: --grm_list is required for --simulate.\n");
	}
	
	// set threads
	omp_set_num_threads(num_threads);
	
	if(zero_grm) {
		std::cout<<"\nZeroing out GRM off-diag elements <"<<zero_cutoff<<std::endl;
		if(binary_grm_file.empty()) throw("\nError: --binary_grm_file is required.\n");
		std::string indi_file = binary_grm_file + ".grm.iid";
		std::vector<std::string> indi_keep;
		get_grm_indi(indi_file, indi_keep);
		MatrixXf grm (indi_keep.size(), indi_keep.size());
		float sumwt;
		read_grm_from_file(binary_grm_file, indi_keep, sumwt, grm);
		grm.triangularView<Lower>() /= sumwt;
		grm.triangularView<StrictlyUpper>() = grm.transpose();
		
		VectorXf gdiag = grm.diagonal();
		grm.diagonal().setOnes();
		int64_t num_zeros = (grm.array() < zero_cutoff).count();
		grm = (grm.array() < zero_cutoff).select(0, grm);
		grm.diagonal() = gdiag;
		sumwt = 1.0f;
		std::cout<<num_zeros<<" GRM off-diag elements were zeroed out."<<std::endl;
		std::cout<<"\nWriting new GRM to "<<output_file<<std::endl;
		write_grm_into_file(output_file, indi_keep, sumwt, grm);
		std::cout<<"Written new GRM to file."<<std::endl;
		
		return;
	}
	
	// merge grms
	if(merge) {
		std::cout<<"\nCombing a list of GRMs into a new binary GRM file...\n"<<std::endl;
		std::string indi_file;
		if(subject_set_file == "") {
			std::ifstream ifs;
			std::string line = "";
			// read individual list file
			ifs.open(grm_list_file);
			if(!ifs.is_open()) {
				std::cout<<"Cannot open "<<grm_list_file<<std::endl;
				exit(1);
			}
			std::getline(ifs, line);
			ifs.clear();
			ifs.close();
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			if(line == "") throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = subject_set_file;
		}
		std::vector<std::string> indi_keep;
		get_grm_indi(indi_file, indi_keep);
		MatrixXf mgrm (indi_keep.size(), indi_keep.size());
		float sumwt;
		merge_grms(grm_list_file, indi_keep, sumwt, mgrm);
		std::cout<<"GRM combining finished.\n"<<std::endl;
		std::cout<<"Writing GRM to "<<output_file<<std::endl;
		write_grm_into_file(output_file, indi_keep, sumwt, mgrm);
		std::cout<<"Written GRM to file."<<std::endl;
		
		return;
	}
	
	// deduct grms
	if(deduct) {
		std::cout<<"\nDeducting GRMs 2+ from GRM 1...\n"<<std::endl;
		std::string indi_file;
		if(subject_set_file == "") {
			std::ifstream ifs;
			std::string line = "";
			// read individual list file
			ifs.open(grm_list_file);
			if(!ifs.is_open()) {
				std::cout<<"Cannot open "<<grm_list_file<<std::endl;
				exit(1);
			}
			std::getline(ifs, line);
			ifs.clear();
			ifs.close();
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			if(line == "") throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = subject_set_file;
		}
		std::vector<std::string> indi_keep;
		get_grm_indi(indi_file, indi_keep);
		MatrixXf mgrm (indi_keep.size(), indi_keep.size());
		float sumwt;
		deduct_grms(grm_list_file, indi_keep, sumwt, mgrm);
		std::cout<<"GRM deduction finished.\n"<<std::endl;
		std::cout<<"Writing GRM to "<<output_file<<std::endl;
		write_grm_into_file(output_file, indi_keep, sumwt, mgrm);
		std::cout<<"Written GRM to file."<<std::endl;
		
		return;
	}
	
	// make core relationship matrices
	if(is_core_set) {
		std::cout<<"\nMaking relationship matrices for COvariance btw Random Effects:\n"<<std::endl;
		std::string indi_file;
		if(subject_set_file == "") {
			std::ifstream ifs;
			std::string line = "";
			// read individual list file
			ifs.open(grm_list_file);
			if(!ifs.is_open()) {
				std::cout<<"Cannot open "<<grm_list_file<<std::endl;
				exit(1);
			}
			std::getline(ifs, line);
			ifs.clear();
			ifs.close();
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			if(line == "") throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = subject_set_file;
		}
		std::vector<std::string> indi_keep;
		get_grm_indi(indi_file, indi_keep);
		
		make_core(grm_list_file, indi_keep, output_file);
		
		return;
	}
	
	// make fore relationship matrices
	if(is_fore_set) {
		std::cout<<"\nMaking relationship matrices for First-Order interactions btw Random Effects:\n"<<std::endl;
		std::string indi_file;
		if(subject_set_file == "") {
			std::ifstream ifs;
			std::string line = "";
			// read individual list file
			ifs.open(grm_list_file);
			if(!ifs.is_open()) {
				std::cout<<"Cannot open "<<grm_list_file<<std::endl;
				exit(1);
			}
			std::getline(ifs, line);
			ifs.clear();
			ifs.close();
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			if(line == "") throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = subject_set_file;
		}
		std::vector<std::string> indi_keep;
		get_grm_indi(indi_file, indi_keep);
		
		make_fore(grm_list_file, indi_keep, output_file);
		
		return;
	}
	
	// simulate phenotypes based on a list of covariance matrices
	if(simu_pheno) {
		std::cout<<"\nSimulating phenotypes based on a list of GRMs...\n"<<std::endl;
		std::string indi_file;
		if(subject_set_file == "") {
			std::ifstream ifs;
			std::string line = "";
			// read individual list file
			ifs.open(grm_list_file);
			if(!ifs.is_open()) {
				std::cout<<"Cannot open "<<grm_list_file<<std::endl;
				exit(1);
			}
			std::getline(ifs, line);
			ifs.clear();
			ifs.close();
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			if(line == "") throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = subject_set_file;
		}
		std::vector<std::string> indi_keep;
		get_grm_indi(indi_file, indi_keep);
		MatrixXf sgrm (indi_keep.size(), indi_keep.size());
		float gvar;
		sum_grms(grm_list_file, indi_keep, gvar, sgrm);
		std::cout<<"\nComputed S = sum(G_i*var_i).\n"<<std::endl;
		
		// code block added June-20-2022
		VectorXf dgg = sgrm.diagonal();
		sgrm.triangularView<Upper>() = sgrm.transpose();
		
		sgrm.diagonal().array() += 1e-4 * gvar;
		LLT<Ref<MatrixXf> > llt(sgrm);
		if(llt.info() == NumericalIssue) {
			throw("S appears not to be positive definite!\n");
		}
		std::cout<<"Completed Cholesky decomposition (L) of S."<<std::endl;
		
		std::mt19937 rng(seed+54321);
		std::normal_distribution<float> rnorm(0.0,1.0);
		MatrixXf rnum(indi_keep.size(), n_pheno*2);
		for(int j=0; j<rnum.cols(); ++j) {
			for(int i=0; i<rnum.rows(); ++i)
				rnum(i,j) = rnorm(rng);
		}
		std::cout<<"Generated IID random vectors (x and y) from Gaussian."<<std::endl;
		
		MatrixXf gv = llt.matrixL() * rnum.leftCols(n_pheno);
		MatrixXf pheno = gv + std::sqrt((1.0-h2[0])/h2[0]*gvar) * rnum.rightCols(n_pheno);
		std::cout<<"Computed p = L*x + sqrt(var_e)*y."<<std::endl;
		
		// code block edited June-21-2022
		int trn_size = sgrm.rows() * 0.9;
		sgrm.diagonal() = dgg.array() + (1.0-h2[0])/h2[0]*gvar;
		sgrm.triangularView<Lower>() = sgrm.transpose();
		std::cout<<"\nComputed V = S + var_e*I.\n";
		MatrixXf py = sgrm.topLeftCorner(trn_size,trn_size).llt().solve(pheno.topRows(trn_size));
		sgrm.diagonal() = dgg;
		MatrixXf egv = sgrm.topRightCorner(trn_size, sgrm.rows()-trn_size).transpose() * py;
		std::cout<<"Completed S(validation x training) * solve(V(training), p(training)).\n";
		
		std::ofstream ofs (output_file + ".sim.csv", std::ofstream::out);
		ofs<<"iid";
		for(int j=0; j<n_pheno; ++j) ofs<<","<<j+1;
		for(int j=0; j<n_pheno; ++j) ofs<<",gv"<<j+1;  // genomic values
		for(int j=0; j<n_pheno; ++j) ofs<<",egv"<<j+1;  // optimal estimated genomic values (added June-20-2022)
		ofs<<"\n";
		for(int i=0; i<indi_keep.size(); ++i) {
			ofs<<indi_keep[i];
			for(int j=0; j<n_pheno; ++j) ofs<<","<<pheno(i,j);
			for(int j=0; j<n_pheno; ++j) ofs<<","<<gv(i,j);  // genomic values
			// optimal estimated genomic values (edited June-21-2022)
			if(i < trn_size) ofs<<std::string(n_pheno, ',');
			else for(int j=0; j<n_pheno; ++j) ofs<<","<<egv(i-trn_size,j);
			ofs<<"\n";
		}
		std::cout<<"Written simulated phenotypes into "<<output_file + ".sim.csv"<<std::endl;
		
		return;
	}
	
	// empirical best linear unbiased predictions
	if(pred_indi) {
		std::cout<<"\nStarted empirical best linear unbiased predictions."<<std::endl;
		
		std::string vc_file = mq_file + ".mq.vc.csv";
		std::vector<std::string> vc_name;
		std::vector<float> vc_init;
		std::ifstream ifs;
		std::string line;
		// read MQ VC file
		ifs.open(vc_file);
		if(!ifs.is_open()) {
			std::cout<<"Cannot open "<<vc_file<<std::endl;
			exit(1);
		}
		std::cout<<"\nReading the MQ VC file from ["+vc_file+"]."<<std::endl;
		std::getline(ifs, line);
		while (std::getline(ifs, line)) {
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			std::vector<std::string> temp = StrFunc::split(line, ',');
			vc_name.push_back(temp[2]);
			vc_init.push_back(std::stof(temp[4]));
		}
		ifs.close();
		if(n_grms == 0) n_grms = vc_name.size()-1;
		if(n_grms > vc_name.size()-1) throw("\nError: --num_grms you specified is too big.\n");
		std::cout<<"There are "<<n_grms<<" GRMs (or general rel. matrices) in the VC file."<<std::endl;
		
		std::string indi_file = vc_name[0] + ".grm.iid";
		std::cout<<"\nReading GRM individual IDs from "<<indi_file<<std::endl;
		std::vector<std::string> indi_grm;
		get_grm_indi(indi_file, indi_grm);
		std::cout<<"There are "<<indi_grm.size()<<" individuals in "<<indi_file<<std::endl;
		
		// read MQ PY file
		std::string py_file = mq_file + ".mq.py.csv";
		std::map<std::string, float> indi2py;
		ifs.open(py_file);
		if(!ifs.is_open()) {
			std::cout<<"Cannot open "<<py_file<<std::endl;
			exit(1);
		}
		std::cout<<"\nReading the MQ PY file from ["+py_file+"]."<<std::endl;
		std::getline(ifs, line);
		while (std::getline(ifs, line)) {
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			std::vector<std::string> temp = StrFunc::split(line, ',');
			indi2py[temp[0]] = std::stof(temp[1]);
		}
		ifs.close();
		std::cout<<"There are "<<indi2py.size()<<" individuals in the PY file."<<std::endl;
		
		int num_training = 0;
		VectorXf Py(indi_grm.size()); Py.setZero();
		for(int i=0; i<indi_grm.size(); ++i) {
			if(indi2py.find(indi_grm[i]) != indi2py.end()) {
				Py(i) = indi2py[indi_grm[i]];
				num_training ++;
			}
		}
		std::cout<<"There are "<<num_training<<" training individuals.\n"<<std::endl;
		
		MatrixXf grm(indi_grm.size(), indi_grm.size());
		MatrixXf blup(indi_grm.size(),n_grms);
		float sumwt;
		for (int i=0; i<n_grms; ++i) {
			line = vc_name[i];
			read_grm_from_file(line, indi_grm, sumwt, grm);
			grm.triangularView<Lower>() /= sumwt;
			std::cout<<"GRM ["<<line<<"] has been obtained."<<std::endl;
			
			blup.col(i) = grm.selfadjointView<Lower>() * (Py * vc_init[i]);
			std::cout<<"EBLUP for ["<<line<<"] has been computed.\n"<<std::endl;
		}
		std::cout<<"Completed empirical best linear unbiased predictions."<<std::endl;
		
		std::ofstream ofs (output_file + ".mq.blup.csv", std::ofstream::out);
		ofs<<"IID,sum";
		for (int i=0; i<n_grms; ++i) ofs<<","<<vc_name[i];
		ofs<<"\n";
		for(int i=0; i<indi_grm.size(); ++i) {
			ofs<<indi_grm[i]<<","<<blup.row(i).sum();
			for (int j=0; j<n_grms; ++j) ofs<<","<<blup(i,j);
			ofs<<"\n";
		}
		ofs.close();
		
		return;
	}
	
	// read phenotype file and covariate file
	std::vector<std::map<std::string, std::pair<float, float> > > indi2pheno_weight(num_traits);
	std::map<std::string, VectorXf > indi2covar;
	if(reml || minq) {
		for(int i=0; i<num_traits; ++i)
			indi2pheno_weight[i] = read_phenotype_file(phen_file, trait_names[i], error_variance_weight_names[i]);
		indi2covar = read_covariate_file(covar_file, covar_names);
	}
	
	if(reml || minq) {
		std::cout<<"\nStarted REML/MINQUE."<<std::endl;

		int i=0, j=0;
		std::cout<<"\n"<<num_traits<< (num_traits > 1 ? " traits are" : " trait is") <<" included in the analysis:\n";
		for(i=0; i<num_traits; ++i) std::cout<<"["<<i<<"] = "<<trait_names[i]<<"\n";
		
		std::vector<std::string> vc_name;
		std::vector<float> vc_init;
		std::ifstream ifs;
		std::string line;
		// read GRM list file
		ifs.open(grm_list_file);
		if(!ifs.is_open()) {
			std::cout<<"Cannot open "<<grm_list_file<<std::endl;
			exit(1);
		}
		std::cout<<"\nReading the GRM list from ["+grm_list_file+"]."<<std::endl;
		while (std::getline(ifs, line)) {
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			std::stringstream  ss(line);
			std::string vc_i;
			float wt = 0.0;
			ss >> vc_i >> wt;
			vc_name.push_back(vc_i);
			vc_init.push_back(wt);
		}
		ifs.close();
		std::cout<<vc_name.size()<< (vc_name.size() > 1 ? " matrices are" : " matrix is") <<" in the list."<<std::endl;
		if(n_grms == 0) n_grms = vc_name.size();
		if(n_grms > vc_name.size()) throw("\nError: --num_grms you specified is > the total number in --grm_list\n");
		if(std::accumulate(vc_init.begin(), vc_init.end(), 0) < 0) throw("\nError: the sum of initial VCs in --grm_list is <0.\n");
		
		std::string indi_file;
		if(subject_set_file == "") {
			indi_file = vc_name[0] + ".grm.iid";
		} else {
			indi_file = subject_set_file;
		}
		
		std::cout<<"\nReading GRM individual IDs from ["+indi_file+"]."<<std::endl;
		std::vector<std::string> indi_grm;
		get_grm_indi(indi_file, indi_grm);
		std::vector<std::string> indi_keep;
		for(i=0; i<indi_grm.size(); ++i) {
			bool any_missing = false;
			if(indi2covar.find(indi_grm[i]) == indi2covar.end())
				any_missing = true;
			for(j=0; j<num_traits; ++j)
				if(indi2pheno_weight[j].find(indi_grm[i]) == indi2pheno_weight[j].end()) any_missing = true;
			if(!any_missing) indi_keep.push_back(indi_grm[i]);
		}
		const int keep_num = indi_keep.size();
		std::cout<<"Non-missing analysis set contains "<<keep_num<<" individuals.\n"<<std::endl;
		if(keep_num < 10) throw("\nError: there are <10 individuals in analysis.\n"); 
		
		VectorXf yvec (keep_num*num_traits);
		VectorXf rvec (keep_num*num_traits);
		for(i=0; i<num_traits; ++i) {
			for(j=0; j<keep_num; ++j) {
				yvec[i*keep_num+j] = indi2pheno_weight[i][indi_keep[j]].first;
				rvec[i*keep_num+j] = indi2pheno_weight[i][indi_keep[j]].second;
			}
		}
		VectorXf pheno_sd(num_traits);
		VectorXf pheno_mean(num_traits);
		for(i=0; i<num_traits; ++i) {
			pheno_mean(i) = yvec.segment(i*keep_num, keep_num).mean();
			yvec.segment(i*keep_num, keep_num).array() -= pheno_mean(i);
			pheno_sd(i) = yvec.segment(i*keep_num, keep_num).norm()/std::sqrt(keep_num);
			yvec.segment(i*keep_num, keep_num) /= pheno_sd(i);
		}
		std::cout<<"Phenotypic standard deviation(s):\n";
		for(i=0; i<num_traits; ++i) std::cout<<"["<<i<<"] = "<<pheno_sd(i)<<"\n";
		std::cout<<"\n";

		const int n_vcs = vc_name.size() + 1;
		const int per_vc = num_traits * (num_traits+1) / 2;
		const int n_params = n_vcs * per_vc;

		// index traits
		MatrixXi param_index(per_vc, 2);
		VectorXi var_index(num_traits);
		i = 0;
		for(int g=0; g<num_traits; ++g) {
			for(int k=0; k<=g; ++k) {
				param_index(i,0) = g;
				param_index(i,1) = k;
				if(g == k) var_index(g) = i;
				i++;
			}
		}

		MatrixXf rmat(keep_num, per_vc);
		for(i=0; i<per_vc; ++i) {
			int g = param_index(i,0);
			int k = param_index(i,1);
			if(g == k) rmat.col(i) = rvec.segment(g*keep_num, keep_num);
			else rmat.col(i) = rvec.segment(g*keep_num, keep_num).array().sqrt() * rvec.segment(k*keep_num, keep_num).array().sqrt();
		}

		VectorXf phen_corr(per_vc);
		for(i=0; i<per_vc; ++i) {
			int g = param_index(i,0);
			int k = param_index(i,1);
			if(g == k) phen_corr(i) = 1;
			else phen_corr(i) = yvec.segment(g*keep_num, keep_num).dot(yvec.segment(k*keep_num, keep_num)) / keep_num;
		}
		if(num_traits > 1) {
			std::cout<<"Phenotypic correlation(s):\n";
			for(i=0; i<per_vc; ++i) {
				int g = param_index(i,0);
				int k = param_index(i,1);
				if(g == k) continue;
				std::cout<<"["<<g<<"]x["<<k<<"] = "<<phen_corr(i)<<"\n";
			}
			std::cout<<"\n";
		}

		MatrixXf xmat (keep_num, covar_names.size());
		for(j=0; j<keep_num; ++j) xmat.row(j) = indi2covar[indi_keep[j]];
		// QR decomposition of covariate matrix to keep independent columns
		ColPivHouseholderQR<MatrixXf> qr(xmat);
		qr.setThreshold(1e-5);
		if(qr.rank() < xmat.cols()) {
			std::cout<<"Warning: the covariate matrix is not of full rank.\n";
			MatrixXf indep_xmat(xmat.rows(), qr.rank());
			std::vector<std::string> indep_covar_names(qr.rank());
			VectorXi indep_indices = qr.colsPermutation().indices().head(qr.rank());
			std::sort(indep_indices.data(), indep_indices.data()+indep_indices.size());
			for(i=0; i<indep_indices.size(); ++i) {
				indep_xmat.col(i) = xmat.col(indep_indices(i));
				indep_covar_names[i] = covar_names[indep_indices(i)];
			}
			
			xmat = indep_xmat;
			covar_names = indep_covar_names;
			std::cout<<"Kept only independent columns.\n\n";
		}
		if(num_traits>1) {
			xmat.conservativeResize(keep_num*num_traits, covar_names.size()*num_traits);
			xmat.rightCols(covar_names.size()*(num_traits-1)).setZero();
			xmat.bottomRows(keep_num*(num_traits-1)).setZero();
			for(i=1; i<num_traits; ++i) xmat.block(keep_num*i, covar_names.size()*i, keep_num, covar_names.size()) = xmat.topLeftCorner(keep_num, covar_names.size());
		}
		
		const int n_packed = (vc_name.size() - 1)/2 + 1;
		MatrixXf dd;
		std::vector<MatrixXf> gg(n_packed);
		MatrixXf grm(keep_num, keep_num);
		VectorXf wt(vc_name.size());
		float sumwt;
		if(!save_mem) {
			dd.resize(keep_num, vc_name.size());
			for (i=0; i<vc_name.size(); ++i) {
				read_grm_from_file(vc_name[i], indi_keep, sumwt, grm);
				
				dd.col(i) = grm.diagonal();
				wt(i) = sumwt;
				if(i % 2 == 0) {
					gg[i/2].resize(keep_num, keep_num);
					#pragma omp parallel for
					for(int j=1; j<keep_num; ++j) 
						gg[i/2].col(j).head(j) = grm.row(j).head(j);
				} else {
					#pragma omp parallel for
					for(int j=0; j<keep_num; ++j) 
						gg[i/2].col(j).tail(keep_num-j) = grm.col(j).tail(keep_num-j);
				}
				
				std::cout<<"GRM ["<<vc_name[i]<<"] has been obtained."<<std::endl;
			}
			#pragma omp parallel for
			for(int i=0; i<vc_name.size(); ++i) {
				dd.col(i) /= wt(i);
				if(i % 2 == 0) gg[i/2].triangularView<StrictlyUpper>() /= wt(i);
				else gg[i/2].triangularView<StrictlyLower>() /= wt(i);
			}
		} else {
			std::cout<<"GRMs will be read on the fly."<<std::endl;
		}
		
		// set initial VCs
		VectorXf param(n_params);
		param.setZero();
		for(j=0; j<per_vc; ++j) {
			if(param_index(j,0) != param_index(j,1)) continue;
			for(i=0; i<n_vcs-1; ++i) param(i+j*n_vcs) = vc_init[i];
			if( param.segment(j*n_vcs,n_vcs-1).sum() >0 ) param.segment(j*n_vcs,n_vcs-1) *= h2[param_index(j,0)]/param.segment(j*n_vcs,n_vcs-1).sum();
			param(j*n_vcs+n_vcs-1) = 1 - param.segment(j*n_vcs,n_vcs-1).sum();
		}
		for(i=0; i<per_vc; ++i) {
			int g = param_index(i,0);
			int k = param_index(i,1);
			if(g == k) continue;
			param.segment(i*n_vcs, n_vcs) = (param.segment(var_index(g)*n_vcs, n_vcs).array() * param.segment(var_index(k)*n_vcs, n_vcs).array()).sqrt() * phen_corr(i);
		}
		std::cout<<"\nInitial VC values:\n";
		for(i=0; i<per_vc; ++i) {
			int g = param_index(i,0);
			int k = param_index(i,1);
			std::cout<<"["<<g<<"]x["<<k<<"] = "<<param.segment(i*n_vcs, n_vcs).transpose()<<"\n";
		}
		
		// logLL terms of the null model
		MatrixXf xrx = xmat.transpose() * rvec.cwiseInverse().asDiagonal() * xmat;
		LLT<Ref<MatrixXf> > chol_xrx(xrx);
		VectorXf xry = xmat.transpose() * ( yvec.cwiseQuotient(rvec) );
		float varp = ( yvec.dot(yvec.cwiseQuotient(rvec)) - xry.dot(chol_xrx.solve(xry)) ) / (yvec.size()-xmat.cols());
		const double null_ldet_r = rvec.array().log().sum() + yvec.size() * std::log(varp);
		const double null_ldet_xrx = 2 * chol_xrx.matrixLLT().diagonal().array().cast<double>().log().sum() - xmat.cols() * std::log(varp);
		const double null_ypy = yvec.size()-xmat.cols();
		xrx.resize(0,0);
		
		// sample random vectors from normal dist
		// Rademacher samples have a poorer performance.
		std::cout<<"\nGenerating "<<n_rand_vec<<" random vectors from "<<rdist<<std::endl;
		std::mt19937 rng(seed+54321);
		MatrixXf Vrand(keep_num*num_traits, n_rand_vec);
		if(rdist == "Gaussian") {
			std::normal_distribution<float> rnorm(0.0,1.0);
			for(int j=0; j<n_rand_vec; ++j) {
				for(int i=0; i<keep_num*num_traits; ++i)
					Vrand(i,j) = rnorm(rng);
			}
		} else {
			Vrand.setOnes();
			std::uniform_real_distribution<float> runif(0.0,1.0);
			for(int j=0; j<n_rand_vec; ++j) {
				for(int i=0; i<keep_num*num_traits; ++i)
					if(runif(rng) < 0.5) Vrand(i,j) = -1.0f;
			}
		}
		
		// projection matrices
		MatrixXf HinvX;
		MatrixXf XtHinvX;
		LDLT<MatrixXf> lltXtHinvX;
		// trace matrices
		MatrixXf PYstar(keep_num*num_traits, n_rand_vec + 1);
		std::vector<MatrixXf> GPYstar(n_params);
		std::vector<MatrixXf> PGPYstar(n_params);
		// indices of lower triangle of S
		std::vector<int> lowerS_idx;
		for(i=0; i<n_params; ++i) {
			for( j=0; j<=i; ++j) {
				lowerS_idx.push_back(i);
				lowerS_idx.push_back(j);
			}
		}
		// MINQUE matrices
		MatrixXf S(n_params, n_params);               // MINQUE LHS
		VectorXf q(n_params);                         // MINQUE RHS
		VectorXf grad(n_params);                      // gradient
		VectorXf grad_fst(n_params);
		VectorXf grad_new(n_params);
		MatrixXf FI(n_params, n_params);              // Fisher information matrix
		// trust region method
		float eta1 = 0.0001, eta2 = 0.99;             // TODO: better values?
		float alpha1 = 0.25, alpha2 = 3;
		float delta = 1e10;
		float rho;
		
		// first iteration
		std::cout<<"\nStarted REML/MINQUE iteration 1"<<std::endl;
		std::cout<<"Updating H = sum(relationship-matrix * varcomp)."<<std::endl;
		grm.resize(keep_num*num_traits, keep_num*num_traits);
		grm.setZero();
		MatrixXf igg;
		if(!save_mem) {
			for(i=0; i<vc_name.size(); ++i) {
				if(i % 2 == 0) {
					for(int l=0; l<per_vc; ++l) {
						int g = param_index(l,0);
						int k = param_index(l,1);
						#pragma omp parallel for
						for(int j=1; j<keep_num; ++j) {
							grm.block(g*keep_num,k*keep_num,keep_num,keep_num).col(j).head(j) += param(l*n_vcs+i) * gg[i/2].col(j).head(j);
						}
					}
				} else {
					for(int l=0; l<per_vc; ++l) {
						int g = param_index(l,0);
						int k = param_index(l,1);
						#pragma omp parallel for
						for(int j=0; j<keep_num-1; ++j) {
							grm.block(g*keep_num,k*keep_num,keep_num,keep_num).col(j).tail(keep_num-j-1) += param(l*n_vcs+i) * gg[i/2].col(j).tail(keep_num-j-1);
						}
					}
				}
			}
			#pragma omp parallel for
			for(int l=0; l<per_vc; ++l) {
				int g = param_index(l,0);
				int k = param_index(l,1);
				grm.block(g*keep_num,k*keep_num,keep_num,keep_num).triangularView<StrictlyLower>() += grm.block(g*keep_num,k*keep_num,keep_num,keep_num).transpose();
				if(g != k) grm.block(g*keep_num,k*keep_num,keep_num,keep_num).triangularView<StrictlyUpper>() = grm.block(g*keep_num,k*keep_num,keep_num,keep_num).transpose();

				grm.block(g*keep_num,k*keep_num,keep_num,keep_num).diagonal() = param(l*n_vcs+(n_vcs-1))*rmat.col(l) + dd * param.segment(l*n_vcs, n_vcs-1);
			}
		} else {
			igg.resize(keep_num, keep_num);
			for (i=0; i<vc_name.size(); ++i) {
				read_grm_from_file(vc_name[i], indi_keep, sumwt, igg, false);
				wt(i) = sumwt;
				for(int l=0; l<per_vc; ++l) {
					int g = param_index(l,0);
					int k = param_index(l,1);
					#pragma omp parallel for
					for(int j=0; j<keep_num; ++j) {
						grm.block(g*keep_num,k*keep_num,keep_num,keep_num).col(j).tail(keep_num-j) += igg.col(j).tail(keep_num-j) * (param(l*n_vcs+i)/sumwt);
					}
				}
			}
			#pragma omp parallel for
			for(int l=0; l<per_vc; ++l) {
				int g = param_index(l,0);
				int k = param_index(l,1);
				if(g != k) grm.block(g*keep_num,k*keep_num,keep_num,keep_num).triangularView<StrictlyUpper>() = grm.block(g*keep_num,k*keep_num,keep_num,keep_num).transpose();
				grm.block(g*keep_num,k*keep_num,keep_num,keep_num).diagonal() += param(l*n_vcs+(n_vcs-1))*rmat.col(l);
			}
		}
		std::cout<<"Trying Cholesky(H)."<<std::endl;
		LLT<Ref<MatrixXf> > llt(grm);
		if(llt.info() == NumericalIssue) {
			throw("\nH is not positive definite. Change initial VCs to rerun.");
		}
		std::cout<<"Cholesky(H) completed normally."<<std::endl;
		if(verbose) std::cout<<"The reciprocal condition number of H = "<<llt.rcond()<<std::endl;
		HinvX = llt.solve(xmat);
		XtHinvX.noalias() = xmat.transpose() * HinvX;
		lltXtHinvX.compute(XtHinvX);
		std::cout<<"Constructed matrices for projecting out fixed effects."<<std::endl;
		PYstar.col(0) = llt.solve(yvec);
		PYstar.rightCols(n_rand_vec) = llt.matrixU().solve(Vrand);
		project_xmat(HinvX,lltXtHinvX,xmat,PYstar);
		std::vector<MatrixXf> mvGPY(num_traits);
		for(int l=0; l<n_params; ++l) {
			GPYstar[l].setZero(num_traits*keep_num, n_rand_vec + 1);
		}
		for(i=0; i<vc_name.size(); ++i) {
			if(save_mem) {
				read_grm_from_file(vc_name[i], indi_keep, sumwt, igg, false);
				if(per_vc == 1) {
					GPYstar[i].noalias() = igg.selfadjointView<Lower>() * (PYstar/sumwt);
				} else {
					for(j=0; j<num_traits; ++j) mvGPY[j].noalias() = igg.selfadjointView<Lower>() * (PYstar.middleRows(j*keep_num,keep_num)/sumwt);
					#pragma omp parallel for
					for(int l=0; l<per_vc; ++l) {
						int g = param_index(l,0);
						int k = param_index(l,1);
						if(g == k) GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num) = mvGPY[k];
						else {
							GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num) = mvGPY[k];
							GPYstar[l*n_vcs+i].middleRows(k*keep_num, keep_num) = mvGPY[g];
						}
					}
				}
			} else {
				gg[i/2].diagonal() = dd.col(i);
				if(per_vc == 1) {
					if(i % 2 == 0) GPYstar[i].noalias() = gg[i/2].selfadjointView<Upper>() * PYstar;
					else GPYstar[i].noalias() = gg[i/2].selfadjointView<Lower>() * PYstar;
				} else {
					for(j=0; j<num_traits; ++j) {
						if(i % 2 == 0) mvGPY[j].noalias() = gg[i/2].selfadjointView<Upper>() * PYstar.middleRows(j*keep_num,keep_num);
						else mvGPY[j].noalias() = gg[i/2].selfadjointView<Lower>() * PYstar.middleRows(j*keep_num,keep_num);
					} 
					#pragma omp parallel for
					for(int l=0; l<per_vc; ++l) {
						int g = param_index(l,0);
						int k = param_index(l,1);
						if(g == k) GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num) = mvGPY[k];
						else {
							GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num) = mvGPY[k];
							GPYstar[l*n_vcs+i].middleRows(k*keep_num, keep_num) = mvGPY[g];
						}
					}
				}
			}
		}
		if(per_vc == 1) {
			GPYstar[i].noalias() = rmat.col(0).asDiagonal() * PYstar;
		} else {
			#pragma omp parallel for
			for(int l=0; l<per_vc; ++l) {
				int g = param_index(l,0);
				int k = param_index(l,1);
				if(g == k) GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num).noalias() = rmat.col(l).asDiagonal() * PYstar.middleRows(k*keep_num,keep_num);
				else {
					GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num).noalias() = rmat.col(l).asDiagonal() * PYstar.middleRows(k*keep_num,keep_num);
					GPYstar[l*n_vcs+i].middleRows(k*keep_num, keep_num).noalias() = rmat.col(l).asDiagonal() * PYstar.middleRows(g*keep_num,keep_num);
				}
			}
		}

		#pragma omp parallel for
		for(i=0; i<n_params; ++i) {
			q[i] = PYstar.col(0).dot(GPYstar[i].col(0));
		}
		std::cout<<"Constructed the RHS of MINQUE equations."<<std::endl;
		
		std::ofstream ofs (output_file + ".mq.iter.csv", std::ofstream::out);
		ofs<<"iter,num_traits,sample_size,num_GRMs,logLL,dLLpred,dogleg_Newton"<<std::endl;
		
		for(int iter=1; iter<=minq_iter; ++iter) {
			if(iter >1) {
				std::cout<<"\nStarted REML/MINQUE iteration "<<iter<<std::endl;
				std::cout<<"Constructed matrices for projecting out fixed effects."<<std::endl;
				std::cout<<"Constructed the RHS of MINQUE equations."<<std::endl;
			}
			for(i=0; i<n_params; ++i) {
				PGPYstar[i] = llt.solve(GPYstar[i]);
				project_xmat(HinvX,lltXtHinvX,xmat,PGPYstar[i]);
			}
			
			#pragma omp parallel for
			for(int k=0; k<lowerS_idx.size(); k+=2) {
				int i = lowerS_idx[k];
				int j = lowerS_idx[k+1];
				// Fisher information
				S(i,j) = (GPYstar[i].rightCols(n_rand_vec).cwiseProduct(PGPYstar[j].rightCols(n_rand_vec))).colwise().sum().mean();
				// Average information
				// S(i,j) = GPYstar[i].col(0).dot(PGPYstar[j].col(0));
			}
			S.triangularView<Upper>() = S.transpose();
			std::cout<<"Constructed the LHS of MINQUE equations."<<std::endl;
			if(S.llt().info() == NumericalIssue) {
				throw("\nFisher information matrix is not positive definite. Change initial VCs to rerun.");
			}
			if(iter == 1) grad = 0.5 * (q - S*param);
			else grad = grad_new;
			FI = 0.5 * S;
			if(iter == 2) grad_fst = grad;
			
			// gradient clipping
			// if(grad.stableNorm() > grad_fst.stableNorm()) grad = grad_fst.stableNorm() / grad.stableNorm() * grad;
			
			std::cout<<"Gradient norm = "<<grad.norm()<<std::endl;
			if(verbose) {
				std::cout<<"Gradient:"<<std::endl;
				std::cout<<grad<<std::endl;
			}
			// must use double for computing determinant
			float logLL = 0.5*null_ldet_r - llt.matrixLLT().diagonal().array().cast<double>().log().sum() + 0.5*(null_ldet_xrx - lltXtHinvX.vectorD().array().cast<double>().log().sum()) + 0.5*(null_ypy - PYstar.col(0).cast<double>().dot(yvec.cast<double>()));
			std::cout<<"Actual logLL = "<<logLL<<std::endl;
			
		trm_step:
			std::cout<<"Try a trust-region step with radius "<<delta<<std::endl;
			VectorXf delta_param_trust(n_params);
			float tau_nr, tau_sd;
			float dLLpred = dogleg(grad, FI, delta, delta_param_trust, tau_nr, tau_sd );
			std::cout<<"Predicted delta logLL = "<<dLLpred<<std::endl;
			param += delta_param_trust;
			// The part of NR must be close to 1 for convergence. (added Sept-19-2022)
			// Changed 0.9999 to 0.995 (Nov-5-2023)
			if((dLLpred < minq_tol && tau_nr>=0.995) || iter == minq_iter) { 
				std::cout<<"Norm(delta VC)/Norm(VC) = "<<delta_param_trust.norm()/param.norm()<<std::endl;
				std::cout<<"Sum(VCs): ";
				for(i=0; i<per_vc; ++i) {
					int g = param_index(i,0);
					int k = param_index(i,1);
					std::cout<<"["<<g<<"]x["<<k<<"] = "<<param.segment(i*n_vcs, n_vcs).sum();
					if(i < per_vc-1) std::cout<<", ";
					else std::cout<<".\n";
				}
				std::cout<<"Completed REML/MINQUE iteration "<<iter<<std::endl;
				
				ofs<<iter<<","<<num_traits<<","<<keep_num<<","<<n_grms<<","<<logLL<<","<<dLLpred<<","<<tau_nr<<std::endl;
				
				if(verbose) {
					std::cout<<"RHS:"<<std::endl;
					std::cout<<q<<std::endl;
					std::cout<<"LHS:"<<std::endl;
					std::cout<<S<<std::endl;
					std::cout<<"VC:"<<std::endl;
					std::cout<<param<<std::endl;
				}
				
				if(dLLpred > minq_tol) {
					std::cout<<"\nWarning: not converged after "<<minq_iter<<" iterations."<<std::endl;
				}
				
				if( iter > 1 && grad.norm() > grad_fst.norm() ) { // TODO: better way?
					grad.cwiseAbs().maxCoeff(&i);
					std::cout<<"\nWarning:"<<std::endl;
					std::cout<<(i % n_vcs == n_vcs-1 ? "err" : vc_name[i % n_vcs])<<" might have exploding gradient."<<std::endl;
					std::cout<<"Initial = "<<grad_fst(i)<<" vs Ending = "<<grad(i)<<std::endl;
					std::cout<<"Can re-run without it."<<std::endl;
				}
				
				break;
			}
			
			std::cout<<"Updating H = sum(relationship-matrix * varcomp)."<<std::endl;
			grm.setZero();
			if(!save_mem) {
				for(i=0; i<vc_name.size(); ++i) {
					if(i % 2 == 0) {
						for(int l=0; l<per_vc; ++l) {
							int g = param_index(l,0);
							int k = param_index(l,1);
							#pragma omp parallel for
							for(int j=1; j<keep_num; ++j) {
								grm.block(g*keep_num,k*keep_num,keep_num,keep_num).col(j).head(j) += param(l*n_vcs+i) * gg[i/2].col(j).head(j);
							}
						}
					} else {
						for(int l=0; l<per_vc; ++l) {
							int g = param_index(l,0);
							int k = param_index(l,1);
							#pragma omp parallel for
							for(int j=0; j<keep_num-1; ++j) {
								grm.block(g*keep_num,k*keep_num,keep_num,keep_num).col(j).tail(keep_num-j-1) += param(l*n_vcs+i) * gg[i/2].col(j).tail(keep_num-j-1);
							}
						}
					}
				}
				#pragma omp parallel for
				for(int l=0; l<per_vc; ++l) {
					int g = param_index(l,0);
					int k = param_index(l,1);
					grm.block(g*keep_num,k*keep_num,keep_num,keep_num).triangularView<StrictlyLower>() += grm.block(g*keep_num,k*keep_num,keep_num,keep_num).transpose();
					if(g != k) grm.block(g*keep_num,k*keep_num,keep_num,keep_num).triangularView<StrictlyUpper>() = grm.block(g*keep_num,k*keep_num,keep_num,keep_num).transpose();

					grm.block(g*keep_num,k*keep_num,keep_num,keep_num).diagonal() = param(l*n_vcs+(n_vcs-1))*rmat.col(l) + dd * param.segment(l*n_vcs, n_vcs-1);
				}
			} else {
				for (i=0; i<vc_name.size(); ++i) {
					read_grm_from_file(vc_name[i], indi_keep, sumwt, igg, false);
					for(int l=0; l<per_vc; ++l) {
						int g = param_index(l,0);
						int k = param_index(l,1);
						#pragma omp parallel for
						for(int j=0; j<keep_num; ++j) {
							grm.block(g*keep_num,k*keep_num,keep_num,keep_num).col(j).tail(keep_num-j) += igg.col(j).tail(keep_num-j) * (param(l*n_vcs+i)/sumwt);
						}
					}
				}
				#pragma omp parallel for
				for(int l=0; l<per_vc; ++l) {
					int g = param_index(l,0);
					int k = param_index(l,1);
					if(g != k) grm.block(g*keep_num,k*keep_num,keep_num,keep_num).triangularView<StrictlyUpper>() = grm.block(g*keep_num,k*keep_num,keep_num,keep_num).transpose();
					grm.block(g*keep_num,k*keep_num,keep_num,keep_num).diagonal() += param(l*n_vcs+(n_vcs-1))*rmat.col(l);
				}
			}
			
			std::cout<<"Trying Cholesky(H)."<<std::endl;
			llt.compute(grm);
			if(llt.info() == Success) {
				std::cout<<"Cholesky(H) completed normally."<<std::endl;
				if(verbose) std::cout<<"The reciprocal condition number of H = "<<llt.rcond()<<std::endl;
				
				HinvX = llt.solve(xmat);
				XtHinvX.noalias() = xmat.transpose() * HinvX;
				lltXtHinvX.compute(XtHinvX);
				
				PYstar.col(0) = llt.solve(yvec);
				PYstar.rightCols(n_rand_vec) = llt.matrixU().solve(Vrand);
				project_xmat(HinvX,lltXtHinvX,xmat,PYstar);

				for(i=0; i<vc_name.size(); ++i) {
					if(save_mem) {
						read_grm_from_file(vc_name[i], indi_keep, sumwt, igg, false);
						if(per_vc == 1) {
							GPYstar[i].noalias() = igg.selfadjointView<Lower>() * (PYstar/sumwt);
						} else {
							for(j=0; j<num_traits; ++j) mvGPY[j].noalias() = igg.selfadjointView<Lower>() * (PYstar.middleRows(j*keep_num,keep_num)/sumwt);
							#pragma omp parallel for
							for(int l=0; l<per_vc; ++l) {
								int g = param_index(l,0);
								int k = param_index(l,1);
								if(g == k) GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num) = mvGPY[k];
								else {
									GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num) = mvGPY[k];
									GPYstar[l*n_vcs+i].middleRows(k*keep_num, keep_num) = mvGPY[g];
								}
							}
						}
					} else {
						gg[i/2].diagonal() = dd.col(i);
						if(per_vc == 1) {
							if(i % 2 == 0) GPYstar[i].noalias() = gg[i/2].selfadjointView<Upper>() * PYstar;
							else GPYstar[i].noalias() = gg[i/2].selfadjointView<Lower>() * PYstar;
						} else {
							for(j=0; j<num_traits; ++j) {
								if(i % 2 == 0) mvGPY[j].noalias() = gg[i/2].selfadjointView<Upper>() * PYstar.middleRows(j*keep_num,keep_num);
								else mvGPY[j].noalias() = gg[i/2].selfadjointView<Lower>() * PYstar.middleRows(j*keep_num,keep_num);
							} 
							#pragma omp parallel for
							for(int l=0; l<per_vc; ++l) {
								int g = param_index(l,0);
								int k = param_index(l,1);
								if(g == k) GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num) = mvGPY[k];
								else {
									GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num) = mvGPY[k];
									GPYstar[l*n_vcs+i].middleRows(k*keep_num, keep_num) = mvGPY[g];
								}
							}
						}
					}
				}
				if(per_vc == 1) {
					GPYstar[i].noalias() = rmat.col(0).asDiagonal() * PYstar;
				} else {
					#pragma omp parallel for
					for(int l=0; l<per_vc; ++l) {
						int g = param_index(l,0);
						int k = param_index(l,1);
						if(g == k) GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num).noalias() = rmat.col(l).asDiagonal() * PYstar.middleRows(k*keep_num,keep_num);
						else {
							GPYstar[l*n_vcs+i].middleRows(g*keep_num, keep_num).noalias() = rmat.col(l).asDiagonal() * PYstar.middleRows(k*keep_num,keep_num);
							GPYstar[l*n_vcs+i].middleRows(k*keep_num, keep_num).noalias() = rmat.col(l).asDiagonal() * PYstar.middleRows(g*keep_num,keep_num);
						}
					}
				}
				#pragma omp parallel for
				for(i=0; i<n_params; ++i) {
					q(i) = PYstar.col(0).dot(GPYstar[i].col(0));
					grad_new(i) = 0.5 * ( q(i) - PYstar.rightCols(n_rand_vec).cwiseProduct(GPYstar[i].rightCols(n_rand_vec)).colwise().sum().mean() );
				}
				std::cout<<"Computed gradient for the step."<<std::endl;
				float dLLapprox = 0.5*(grad+grad_new).dot(delta_param_trust);
				rho = dLLapprox / dLLpred;
				if(grad_new.norm() > 2.0f * grad.norm() ) {    // TODO: decide better values
					std::cout<<"Gradient norm becomes too large given the step."<<std::endl;
					rho = -1;
				}
				if(verbose) {
					std::cout<<"Gradient provided the step:"<<std::endl;
					std::cout<<grad_new<<std::endl;
				}
			} else {
				std::cout<<"H is not PSD given the step."<<std::endl;
				rho = -1;
			}
			
			if(rho < eta1) {
				delta = alpha1 * FI.diagonal().cwiseProduct(delta_param_trust).norm();
			} else if(rho < eta2) {
				
			} else {
				delta = std::max(delta, alpha2 * FI.diagonal().cwiseProduct(delta_param_trust).norm());
			}
			
			if(rho > eta1) {
				std::cout<<"Accept the step."<<std::endl;
			} else {
				std::cout<<"Reject the step."<<std::endl;
				std::cout<<"Reduced trust-region radius.\n"<<std::endl;
				param -= delta_param_trust;
				goto trm_step;
			}
			
			std::cout<<"Approx(delta logLL)/predicted(delta logLL) = "<<rho<<std::endl;
			std::cout<<"Norm(delta VC)/Norm(VC) = "<<delta_param_trust.norm()/param.norm()<<std::endl;
			std::cout<<"Sum(VCs): ";
			for(i=0; i<per_vc; ++i) {
				int g = param_index(i,0);
				int k = param_index(i,1);
				std::cout<<"["<<g<<"]x["<<k<<"] = "<<param.segment(i*n_vcs, n_vcs).sum();
				if(i < per_vc-1) std::cout<<", ";
				else std::cout<<".\n";
			}
			std::cout<<"Completed REML/MINQUE iteration "<<iter<<std::endl;
			
			ofs<<iter<<","<<num_traits<<","<<keep_num<<","<<n_grms<<","<<logLL<<","<<dLLpred<<","<<tau_nr<<std::endl;
			
			if(verbose) {
				std::cout<<"RHS:"<<std::endl;
				std::cout<<q<<std::endl;
				std::cout<<"LHS:"<<std::endl;
				std::cout<<S<<std::endl;
				std::cout<<"VC:"<<std::endl;
				std::cout<<param<<std::endl;
			}
		}
		ofs.close();
		grm.resize(0,0);
		igg.resize(0,0);
		for(i=0; i<gg.size(); ++i) gg[i].resize(0,0);
		
		MatrixXf covar_est = FI.ldlt().solve(MatrixXf::Identity(n_params,n_params)) * (1+1.0f/n_rand_vec);
		VectorXf mm = wt.head(n_grms);
		VectorXf enrichment_score(n_params);
		MatrixXf enrichment_var(n_params, n_params);
		VectorXf pve(n_params);
		MatrixXf pve_var(n_params, n_params);
		enrichment_score.setZero();
		enrichment_var.setZero();
		pve.setZero();
		pve_var.setZero();

		for(i=0; i<per_vc; ++i) {
			VectorXf tmpvc = param.segment(n_vcs*i, n_grms);
			enrichment_score.segment(n_vcs*i, n_grms) = mm.sum()/tmpvc.sum()*tmpvc.array()/mm.array();
			MatrixXf tmpvar = mm.sum()/tmpvc.sum() * mm.asDiagonal().inverse() * (MatrixXf::Identity(n_grms,n_grms) - tmpvc/tmpvc.sum() * RowVectorXf::Ones(n_grms));
			if(n_grms > 1) enrichment_var.block(n_vcs*i, n_vcs*i, n_grms, n_grms) = tmpvar * covar_est.block(n_vcs*i, n_vcs*i, n_grms, n_grms) * tmpvar.transpose();

			tmpvc = param.segment(n_vcs*i, n_vcs);
			pve.segment(n_vcs*i, n_vcs) = tmpvc/tmpvc.sum();
			tmpvar = (MatrixXf::Identity(n_vcs,n_vcs) - tmpvc/tmpvc.sum() * RowVectorXf::Ones(n_vcs))/tmpvc.sum();
			pve_var.block(n_vcs*i, n_vcs*i, n_vcs, n_vcs) = tmpvar * covar_est.block(n_vcs*i, n_vcs*i, n_vcs, n_vcs) * tmpvar.transpose();
		}
		
		// scale param and covar_est
		VectorXf scaling(n_params);
		for(i=0; i<per_vc; ++i) {
			int g = param_index(i, 0);
			int k = param_index(i, 1);
			scaling.segment(i*n_vcs, n_vcs).array() = pheno_sd(g) * pheno_sd(k);
		}
		param.array() *= scaling.array();
		covar_est = scaling.asDiagonal() * covar_est * scaling.asDiagonal();
		// write VC into a file
		ofs.open(output_file + ".mq.vc.csv", std::ofstream::out);
		ofs<<"trait_x,trait_y,vc_name,m,var,seV,pve,seP,enrichment,seE";
		for(i=0; i<n_grms; ++i) ofs<<","<<vc_name[i];  // sampling variance of enrichment scores
		for(j=0; j<per_vc; ++j) {
			for(i=0; i<vc_name.size(); ++i) ofs<<","<<vc_name[i];  // sampling variance of VC estimates
			ofs<<",err";
		}
		ofs<<"\n";
		for(int k=0; k<per_vc; ++k) {
			for(i=0; i<n_vcs; ++i) {
				ofs<<trait_names[param_index(k, 0)]<<","<<trait_names[param_index(k, 1)]<<",";
				if(i == n_vcs-1) ofs<<"err,NA";
				else ofs<<vc_name[i]<<","<<std::setprecision(8)<<wt(i);
				int l = i+k*n_vcs;
				ofs<<std::setprecision(-1)<<","<<param(l)<<","<<std::sqrt(covar_est(l,l))<<","<<pve(l)<<","<<std::sqrt(pve_var(l,l));
				if(enrichment_score(l) == 0) ofs<<",NA,NA";
				else if(enrichment_var(l,l) == 0) ofs<<","<<enrichment_score(l)<<",NA";
				else ofs<<","<<enrichment_score(l)<<","<<std::sqrt(enrichment_var(l,l));

				for( j=0; j<n_grms; ++j) {
					if(enrichment_var(l,j+k*n_vcs) == 0) ofs<<",NA";
					else ofs<<","<<enrichment_var(l,j+k*n_vcs);
				}
				for(int g=0; g<per_vc; ++g) {
					for( j=0; j<n_vcs; ++j) 
						ofs<<","<<covar_est(l,j+g*n_vcs);
				}
				ofs<<"\n";
			}
		}
		ofs.close();

		// correlations between traits
		if(num_traits > 1) {
			float x1, x2, x12, var_x1, var_x2, var_x12, cov_x1_x2, cov_x1_x12, cov_x2_x12;
			float cor_val, cor_se;

			ofs.open(output_file + ".mq.cor.csv", std::ofstream::out);
			ofs<<"vc_name,trait_x,trait_y,cor,se\n";
			for(j=0; j<n_vcs; ++j) {
				for(i=0; i<per_vc; ++i) {
					int g = var_index(param_index(i, 0));
					int k = var_index(param_index(i, 1));
					if(g == k) continue;

					x1 = param(g*n_vcs+j);
					x2 = param(k*n_vcs+j);
					x12 = param(i*n_vcs+j);
					var_x1 = covar_est(g*n_vcs+j, g*n_vcs+j);
					var_x2 = covar_est(k*n_vcs+j, k*n_vcs+j);
					var_x12 = covar_est(i*n_vcs+j, i*n_vcs+j);
					cov_x1_x2 = covar_est(g*n_vcs+j, k*n_vcs+j);
					cov_x1_x12 = covar_est(g*n_vcs+j, i*n_vcs+j);
					cov_x2_x12 = covar_est(k*n_vcs+j, i*n_vcs+j);
					
					cor_val = x12 / sqrt(x1 * x2);
					cor_se = std::sqrt( cor_val*cor_val * ( var_x1/(4*x1*x1) + var_x2/(4*x2*x2) + var_x12/(x12*x12) + cov_x1_x2/(2*x1*x2) - cov_x1_x12/(x1*x12) - cov_x2_x12/(x2*x12) ) );
					ofs<<(j == n_vcs-1 ? "err" : vc_name[j])<<","<<trait_names[param_index(i, 0)]<<","<<trait_names[param_index(i, 1)]<<","<<cor_val<<","<<cor_se<<"\n";
				}
			}
			// genetic correlation using all n_grms VCs
			for(i=0; i<per_vc; ++i) {
				int g = var_index(param_index(i, 0));
				int k = var_index(param_index(i, 1));
				if(g == k) continue;
				
				x1 = param.segment(g*n_vcs, n_grms).sum();
				x2 = param.segment(k*n_vcs, n_grms).sum();
				x12 = param.segment(i*n_vcs, n_grms).sum();
				var_x1 = covar_est.block(g*n_vcs, g*n_vcs, n_grms, n_grms).sum();
				var_x2 = covar_est.block(k*n_vcs, k*n_vcs, n_grms, n_grms).sum();
				var_x12 = covar_est.block(i*n_vcs, i*n_vcs, n_grms, n_grms).sum();
				cov_x1_x2 = covar_est.block(g*n_vcs, k*n_vcs, n_grms, n_grms).sum();
				cov_x1_x12 = covar_est.block(g*n_vcs, i*n_vcs, n_grms, n_grms).sum();
				cov_x2_x12 = covar_est.block(k*n_vcs, i*n_vcs, n_grms, n_grms).sum();
				
				cor_val = x12 / sqrt(x1 * x2);
				cor_se = std::sqrt( cor_val*cor_val * ( var_x1/(4*x1*x1) + var_x2/(4*x2*x2) + var_x12/(x12*x12) + cov_x1_x2/(2*x1*x2) - cov_x1_x12/(x1*x12) - cov_x2_x12/(x2*x12) ) );
				ofs<<"G,"<<trait_names[param_index(i, 0)]<<","<<trait_names[param_index(i, 1)]<<","<<cor_val<<","<<cor_se<<"\n";
			}
			ofs.close();
		}
		
		// scale Py by pheno_sd
		// write Py into a file
		for(i=0; i<num_traits; ++i) PYstar.col(0).segment(i*keep_num,keep_num) /= pheno_sd(i);
		ofs.open(output_file + ".mq.py.csv", std::ofstream::out);
		ofs<<"IID";
		for(i=0; i<num_traits; ++i) ofs<<","<<trait_names[i];
		ofs<<"\n";
		for(i=0; i<keep_num; ++i) {
			ofs<<indi_keep[i];
			for(j=0; j<num_traits; ++j) ofs<<","<<PYstar(j*keep_num+i,0);
			ofs<<"\n";
		}
		ofs.close();
		
		// calculate shift adjustment
		xmat.conservativeResize(keep_num, covar_names.size());
		VectorXf adj_mean = (xmat.transpose() * xmat).llt().solve(xmat.transpose() * VectorXf::Ones(keep_num));
		VectorXf adj_pred = xmat * adj_mean;
		bool shifting = false;
		if(fabs(adj_pred.maxCoeff()-1.0) < 0.0001 && fabs(adj_pred.minCoeff()-1.0) < 0.0001) shifting = true;

		// scale BLUE and SE by pheno_sd
		// write BLUE into a file
		VectorXf blue = lltXtHinvX.solve(HinvX.transpose() * yvec);
		MatrixXf blue_var = lltXtHinvX.solve( MatrixXf::Identity(covar_names.size()*num_traits,covar_names.size()*num_traits) );
		ofs.open(output_file + ".mq.blue.csv", std::ofstream::out);
		ofs<<"trait,covar,blue,se,pval";
		for(i=0; i<num_traits; ++i) 
			for(j=0; j<covar_names.size(); ++j)  ofs<<","<<trait_names[i]<<"."<<covar_names[j];
		ofs<<"\n";
		IOFormat CommaInitFmt(FullPrecision, 0, ",");
		RowVectorXf row_scale(blue.size());
		for(i=0; i<num_traits; ++i) row_scale.segment(i*covar_names.size(), covar_names.size()).setConstant(pheno_sd(i));
		for(i=0; i<num_traits; ++i) {
			if(shifting) blue.segment(i*covar_names.size(), covar_names.size()) += adj_mean * pheno_mean(i) / pheno_sd(i);
			for(j=0; j<covar_names.size(); ++j) {
				int l = i*covar_names.size()+j;
				ofs<<trait_names[i]<<","<<covar_names[j]<<","<<blue(l)*pheno_sd(i)<<","<<std::sqrt(blue_var(l,l))*pheno_sd(i)<<","<<getOneDfChisqPval(blue(l)*blue(l)/blue_var(l,l))<<",";

				RowVectorXf row = blue_var.row(l).array() * row_scale.array() * pheno_sd(i);
				ofs<<row.format(CommaInitFmt)<<"\n";
			}
		}
		ofs.close();

		std::cout<<"\nCompleted REML/MINQUE."<<std::endl;
		
		return;
	}
	
	// make GRM
	if(make_grm) {
		if(grm_type=='A') {
			std::cout<<"\nStarted to make additive GRM.\n";
			std::cout<<"GRM = ZWZ', where Z = standardized genotypes and W = weights.\n"<<std::endl;
		} else {
			std::cout<<"\nStarted to make dominance GRM.\n";
			std::cout<<"GRM = ZWZ', where Z = standardized genotypes for dominance deviation.\n"<<std::endl;
		}
		
		if( !(file_check(binary_genotype_file + ".fam") && file_check(binary_genotype_file + ".bim") && file_check(binary_genotype_file + ".bed")) )
			throw("\nError: PLINK bed/bim/fam files [" + binary_genotype_file + "] not found\n");
		
		if(subject_set_file != "") 
			std::cout<<"Set the subject set file ["<<subject_set_file<<"]\n";
		
		// get subjects from --subject_set (if any) and plink fam
		std::vector<bool> bindi;
		std::vector<std::string> indi_keep;
		get_subject_set(subject_set_file, binary_genotype_file+".fam", bindi, indi_keep);
		int keep_num = indi_keep.size();
		std::cout<<bindi.size()<<" individuals are found in ["<<binary_genotype_file<<".fam].\n";
		std::cout<<keep_num<<" individuals are retained for making GRM."<<std::endl;
		if(keep_num == 0) throw("\nError: subject set retained for analysis is empty.\n");
		
		// get SNPs from --snp_info_file and plink bim
		// read marker prior file
		std::map<std::string, std::pair<std::string, float> > marker2group_weight;
		read_marker_info_file(marker_info_file, marker_group, marker_effvar_weight_name, marker2group_weight);
		std::vector<bool> bmarker;
		std::vector<double> gvec;
		get_marker_set(marker2group_weight, binary_genotype_file+".bim", bmarker, gvec);
		Eigen::Map<VectorXd> gwt(gvec.data(), gvec.size());
		int pre_qc_marker_num = gwt.size();
		std::cout<<bmarker.size()<<" variants are found in ["<<binary_genotype_file<<".bim].\n";
		std::cout<<pre_qc_marker_num<<" variants are retained before quality control."<<std::endl;
		if(pre_qc_marker_num == 0) throw("\nError: marker set retained for analysis is empty.\n");
		
		// variables for making GRM
		int subset_size = 50e6/keep_num;
		const int num_sets = pre_qc_marker_num/subset_size;
		const int rest_size = pre_qc_marker_num % subset_size;
		MatrixXc kmat(keep_num, subset_size);
		MatrixXd SS = MatrixXd::Zero(keep_num, subset_size);
		std::cout<<"\nSplit SNPs into "<<num_sets + (rest_size>0 ? 1:0)<<" subsets.\n";
		MatrixXd grm = MatrixXd::Zero(keep_num, keep_num);
		double sumwt = 0;
		int post_qc_marker_num = 0;
		VectorXd swt(subset_size);
		VectorXd hwep(subset_size);
		
		//read binary genotype file
		std::string plink_bed_file = binary_genotype_file + ".bed";
		std::ifstream ifs;
		ifs.open(plink_bed_file, std::ifstream::binary);
		std::cout<<"Reading PLINK BED file from ["+plink_bed_file+"] in SNP-major format."<<std::endl;
		unsigned snp_num = bmarker.size();
		unsigned indi_num = bindi.size();
		unsigned i, j, k;
		char ch[1];
		std::bitset<8> b;
		for(i=0; i<3; i++) ifs.read(ch,1); // skip the first three bytes
		unsigned snp_indx=0, indi_indx=0;
		unsigned km_idx = 0;
		for(j=0, snp_indx=0; j<snp_num; j++){
		// Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
			if(!bmarker[j]){
				for(i=0; i<indi_num; i+=4) ifs.read(ch,1);
				continue;
			}
			for(i=0, indi_indx=0; i<indi_num;){
				ifs.read(ch,1);
				if(!ifs) throw("\nError: problem with the BED file. Has the FAM/BIM file been changed?\n");
				b=ch[0];
				k=0;
				while(k < 7 && i < indi_num){
					if(!bindi[i]) k+=2;
					else{
						kmat(indi_indx, km_idx) = b[k] + b[k + 1];
						k += 2;
						indi_indx++;
					}
					i++;
				}
			}
			km_idx++;
			snp_indx++;
			if(km_idx == subset_size) {
				calc_hwep_geno012(kmat, hwe_midp, hwep);
				swt = gwt.segment(snp_indx-subset_size,subset_size);
				calc_grm_by_subset(grm_type, min_maf, min_hwep, SS, kmat, hwep, swt, grm, sumwt, post_qc_marker_num);
				
				km_idx = 0;
				
				std::cout<<"Completed subset "<<snp_indx/subset_size<<std::endl; 
			}
			if(snp_indx==pre_qc_marker_num) break;
		}
		ifs.clear();
		ifs.close();
		
		if(rest_size>0) {
			kmat.conservativeResize(NoChange, rest_size);
			SS.resize(NoChange, rest_size);
			hwep.resize(rest_size);
			
			calc_hwep_geno012(kmat, hwe_midp, hwep);
			swt = gwt.tail(rest_size);
			calc_grm_by_subset(grm_type, min_maf, min_hwep, SS, kmat, hwep, swt, grm, sumwt, post_qc_marker_num);
			
			std::cout<<"Completed subset "<<num_sets+1<<std::endl; 
		}
		SS.resize(0,0);
		kmat.resize(0,0);
		
		std::cout<<"\nCompleted making "<<(grm_type=='A' ? "additive" : "dominance")<<" GRM.\n";
		std::cout<<post_qc_marker_num<<" post-QC variants are used for making GRM.\n";
		
		std::cout<<"\nWriting GRM to file..."<<std::endl;
		write_dGRM_into_file((grm_type=='A' ? output_file : output_file+".d"), indi_keep, (float)sumwt, grm);
		std::cout<<"Written GRM to file."<<std::endl;
		
		return;
	}
}

void print_help() {
	printf("\nPlease see online documentation at https://jiang18.github.io/mph/\n");
}

void project_xmat(
	const MatrixXf &HinvX,
	const LDLT<MatrixXf>& lltXtHinvX,
	const MatrixXf &X,
	Ref<MatrixXf> Y)
{
	Y = Y - HinvX * lltXtHinvX.solve(X.transpose() * Y);
}

// y = a*x^2 + b*x + c
void optim_quad (
	const float a,
	const float b,
	const float c,
	const float xlb,
	const float xub,
	float &xmin,
	float &xmax,
	float &ymin,
	float &ymax
 )
{
	float x0 = -b/(2*a);
	float y0 = a*x0*x0 + b*x0 + c;
	float ylb = a*xlb*xlb + b*xlb + c;
	float yub = a*xub*xub + b*xub + c;
	
	if(xlb < x0 && x0 < xub) {
		if(y0 < ylb) {
			xmin = x0;
			ymin = y0;
			if(ylb < yub) {
				xmax = xub;
				ymax = yub;
			} else {
				xmax = xlb;
				ymax = ylb;
			}
		} else {
			xmax = x0;
			ymax = y0;
			if(ylb < yub) {
				xmin = xlb;
				ymin = ylb;
			} else {
				xmin = xub;
				ymin = yub;
			}
		}
		return;
	}
	
	if(ylb < yub) {
		xmin = xlb;
		xmax = xub;
		ymin = ylb;
		ymax = yub;
	} else {
		xmin = xub;
		xmax = xlb;
		ymin = yub;
		ymax = ylb;
	}
}

// maximize p' * grad - 0.5 * p' * FI *p
// s.t. ||diag(FI) * p|| <= delta
float dogleg (
	const VectorXf &grad,
	const MatrixXf &FI,
	const float delta,
	Ref<VectorXf> p,
	float &tau_nr,
	float &tau_sd )
{
	float ret_f;
// Given diag = diag(FI), 
// p = diag * p
// g = inv(diag)*(-grad)
// B = inv(diag)*FI*inv(diag)
// minimize p' * g + 0.5 * p' * B *p
// s.t. ||p|| <= delta
	VectorXf diag = FI.diagonal();
	VectorXf g = -grad.cwiseQuotient(diag);
	MatrixXf B = diag.cwiseInverse().asDiagonal() * FI * diag.cwiseInverse().asDiagonal();
	
	VectorXf pu = -g.dot(g)/g.dot(B*g) * g;
	VectorXf pb = -B.fullPivLu().solve(g);
	VectorXf pd = pb - pu;
	
// 0<=tau<=1, funcU_tau = tau*pu.dot(g) + 0.5 * tau^2 * pu.dot(B * pu)
	float tauU_min, tauU_max, funcU_min, funcU_max;
	optim_quad(0.5 * pu.dot(B * pu), pu.dot(g), 0, 0, std::min(1.0f, delta/pu.norm()), tauU_min, tauU_max, funcU_min, funcU_max);
	
// 0<=tau<=1, p_tau = pu + tau * pd
// tau^2 * pd.dot(pd) + 2*tau*pd.dot(pu) + pu.dot(pu) - delta^2 <= 0
// funcB_tau = 0.5 * tau^2 * pd.dot(B * pd) + tau*[pu.dot(B*pd)+pd.dot(g)] + 0.5 * pu.dot(B * pu) + pu.dot(g)
	float constrain_a = pd.dot(pd);
	float constrain_b = 2*pd.dot(pu);
	float constrain_c = pu.dot(pu) - delta*delta;
	float constrain_delta = std::pow(constrain_b, 2) - 4*constrain_a*constrain_c;
	
	float funcB_a = 0.5 * pd.dot(B * pd);
	float funcB_b = pu.dot(B*pd)+pd.dot(g);
	float funcB_c = 0.5 * pu.dot(B * pu) + pu.dot(g);
	float tauB_min, funcB_min, tauB_max, funcB_max;
	if(constrain_delta < 0){
		p = tauU_min * pu.cwiseQuotient(diag);
		ret_f = -funcU_min;
		tau_sd = tauU_min;
		tau_nr = 0;
	} else if(constrain_delta == 0) {
		tauB_min = -constrain_b/(2*constrain_a);
		funcB_min = tauB_min*tauB_min*funcB_a + tauB_min*funcB_b + funcB_c;
		if(funcB_min < funcU_min) {
			p = (pu + tauB_min*pd).cwiseQuotient(diag);
			ret_f = -funcB_min;
			tau_sd = 1-tauB_min;
			tau_nr = tauB_min;
		} else {
			p = tauU_min * pu.cwiseQuotient(diag);
			ret_f = -funcU_min;
			tau_sd = tauU_min;
			tau_nr = 0;
		}
	} else {
		float tauB_lb = (-constrain_b - std::sqrt(constrain_delta))/(2*constrain_a);
		float tauB_ub = (-constrain_b + std::sqrt(constrain_delta))/(2*constrain_a);
		
		if(tauB_ub <= 0 || tauB_lb > 1) {
			p = tauU_min * pu.cwiseQuotient(diag);
			ret_f = -funcU_min;
			tau_sd = tauU_min;
			tau_nr = 0;
		} else {
			optim_quad(funcB_a, funcB_b, funcB_c, std::max(0.0f,tauB_lb), std::min(1.0f, tauB_ub), tauB_min, tauB_max, funcB_min, funcB_max);
			
			if(funcB_min < funcU_min) {
				p = (pu + tauB_min*pd).cwiseQuotient(diag);
				ret_f = -funcB_min;
				tau_sd = 1-tauB_min;
				tau_nr = tauB_min;
			} else {
				p = tauU_min * pu.cwiseQuotient(diag);
				ret_f = -funcU_min;
				tau_sd = tauU_min;
				tau_nr = 0;
			}
		}
	}
	std::cout<<"Step = "<<tau_sd<<" steepest descent + "<<tau_nr<<" Newton-Raphson."<<std::endl;
	return(ret_f);
}
