#include <bitset>
#include <getopt.h>
#include <random>
#include "mph.h"
#include "reml.h"

int num_threads;

void print_help();
void option(int option_num, char **option_str);

int main(int argc, char **argv)
{   
	std::cout<<"***************************************************************************"<<std::endl;
	std::cout<<"* MPH by Jicai Jiang"<<std::endl;
	std::cout<<"* MINQUE for Partitioning Heritability"<<std::endl;
	std::cout<<"* Version 0.55.0 (July 20, 2025)"<<std::endl;
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
	bool subset_grm = false;
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
	std::string keep_file="";
	std::string phen_file="", covar_file="", str_traits="", str_error_variance_weights="", str_covar_names="";
	std::vector<std::string> trait_names;
	std::vector<std::string> error_variance_weight_names;
	std::vector<std::string> covar_names;
	std::string marker_info_file="", marker_effvar_weight_name="", str_geno_coding="";
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
		{"subset_grm",  required_argument, 0, 'U'},
		{"num_phenotypes",  required_argument, 0, 'x'},
		{"num_iterations",  required_argument, 0, 'p'},
		{"tolerance",  required_argument, 0, 'q'},
		{"num_random_vectors",  required_argument, 0, 'a'},
		{"seed",  required_argument, 0, 'k'},
		{"distribution",  required_argument, 0, 's'},
		{"num_grms",  required_argument, 0, 'd'},
		{"heritability",  required_argument, 0, 'h'},
		{"keep",  required_argument, 0, 'i'},
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
		{"snp_genotype_coding",  required_argument, 0, 'r'},
		{"snp_weight_name",  required_argument, 0, 'w'},
		{"binary_genotype_file",  required_argument, 0, 'g'},
		{"bfile",  required_argument, 0, 'b'},
		{"num_threads",    required_argument, 0, 'n'},
		{"output_file",    required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};
	int option_index = 0;
	int opt;
	while ((opt = getopt_long(option_num, option_str, "VMARSPGHICFDUx:p:q:a:k:s:d:h:i:2:4:6:8:j:f:t:e:c:v:m:r:w:g:b:n:o:",long_options, &option_index)) != -1)
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
			case 'U':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				subset_grm = true;
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
			case 'h':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				str_h2 = optarg;
				break;
			case 'i':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				keep_file = optarg;
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
				str_geno_coding = optarg;
				break;
			case 'w':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				marker_effvar_weight_name = optarg;
				break;
			case 'g':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				binary_genotype_file = optarg;
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
	if(min_maf < 0 || min_maf >= 0.5) throw("\nError: min MAF is out of [0, 0.5).\n");
	if(min_hwep < 0 || min_hwep >= 1) throw("\nError: min HWE P value is out of [0, 1).\n");
	if(num_threads > 100 || num_threads < 1) throw("\nError: number of threads is out of [1,100].\n");
	if(output_file.empty()) throw("\nError: output file is required for analysis.\n");
	if(subset_grm) {
		if(binary_grm_file.empty()) throw("\nError: --binary_grm_file is required for --subset_grm.\n");
		if(keep_file.empty()) throw("\nError: --keep is required for --subset_grm.\n");
	}
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
	
		if(grm_list_file.empty() && binary_grm_file.empty()) 
			throw("\nError: --grm_list and/or --binary_grm_file must be set for --reml/--minque.\n");
		if(n_grms <0) throw("\nError: --num_grms must be >=0.\n");
		if(minq_tol < 1e-4) minq_tol = 1e-4;
		if(minq_iter < 1) throw("\nError: --num_iterations must be >0.\n");
		if(rdist != "Gaussian" && rdist != "Rademacher")
			throw("\nError: --distribution must be either Gaussian or Rademacher.\n");
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

	if(make_grm) {
		if(marker_info_file.empty())
			throw("\nError: --snp_info_file must be set for --make_grm.\n");
		if(binary_genotype_file.empty())
			throw("\nError: --bfile must be set for --make_grm.\n");
	}
	if(grm_list_file.empty()) {
		if(merge) throw("\nError: --grm_list is required for --merge_grms.\n");
		if(deduct) throw("\nError: --grm_list is required for --deduct_grms.\n");
		if(is_core_set) throw("\nError: --grm_list is required for --make_core.\n");
		if(is_fore_set) throw("\nError: --grm_list is required for --make_fore.\n");
		if(simu_pheno) throw("\nError: --grm_list is required for --simulate.\n");
	}
	
	// set threads
	omp_set_num_threads(num_threads);
	
	if(subset_grm) {
		std::cout<<"\nSubsetting GRM ["<<binary_grm_file<<"]...\n"<<std::endl;
		
		std::vector<std::string> indi_keep;
		get_grm_indi(keep_file, indi_keep);
		std::cout<<"Keeping "<<indi_keep.size()<<" individuals from ["<<keep_file<<"] in the specified order.\n"<<std::endl;
		
		MatrixXf grm (indi_keep.size(), indi_keep.size());
		float sumwt;
		read_grm_from_file(binary_grm_file, indi_keep, sumwt, grm);
		
		std::cout<<"Writing subset GRM to "<<output_file<<std::endl;
		write_grm_into_file(output_file, indi_keep, sumwt, grm);
		std::cout<<"Subset GRM written to file."<<std::endl;
		
		return;
	}
	
	// merge grms
	if(merge) {
		std::cout<<"\nCombing a list of GRMs into a new binary GRM file...\n"<<std::endl;
		std::string indi_file;
		if(keep_file.empty()) {
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
			if(line.empty()) throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = keep_file;
		}
		std::vector<std::string> indi_keep;
		get_grm_indi(indi_file, indi_keep);
		MatrixXf mgrm (indi_keep.size(), indi_keep.size());
		float sumwt;
		merge_grms(grm_list_file, indi_keep, sumwt, mgrm);
		std::cout<<"GRM combining finished.\n"<<std::endl;
		std::cout<<"Writing GRM to "<<output_file<<std::endl;
		write_grm_into_file(output_file, indi_keep, sumwt, mgrm);
		std::cout<<"GRM written to file."<<std::endl;
		
		return;
	}
	
	// deduct grms
	if(deduct) {
		std::cout<<"\nDeducting GRMs 2+ from GRM 1...\n"<<std::endl;
		std::string indi_file;
		if(keep_file.empty()) {
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
			if(line.empty()) throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = keep_file;
		}
		std::vector<std::string> indi_keep;
		get_grm_indi(indi_file, indi_keep);
		MatrixXf mgrm (indi_keep.size(), indi_keep.size());
		float sumwt;
		deduct_grms(grm_list_file, indi_keep, sumwt, mgrm);
		std::cout<<"GRM deduction finished.\n"<<std::endl;
		std::cout<<"Writing GRM to "<<output_file<<std::endl;
		write_grm_into_file(output_file, indi_keep, sumwt, mgrm);
		std::cout<<"GRM written to file."<<std::endl;
		
		return;
	}
	
	// make core relationship matrices
	if(is_core_set) {
		std::cout<<"\nMaking relationship matrices for COvariance btw Random Effects:\n"<<std::endl;
		std::string indi_file;
		if(keep_file.empty()) {
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
			if(line.empty()) throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = keep_file;
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
		if(keep_file.empty()) {
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
			if(line.empty()) throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = keep_file;
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
		if(keep_file.empty()) {
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
			if(line.empty()) throw("\nError: the first line of the GRM list file is empty.\n");
			std::stringstream ss(line);
			std::string vc_i;
			ss >> vc_i;
			indi_file = vc_i + ".grm.iid";
		} else {
			indi_file = keep_file;
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
	
	// REML/MINQUE
	if(reml || minq) {
		// read phenotype file and covariate file
		std::vector<std::map<std::string, std::pair<float, float> > > indi2pheno_weight(num_traits);
		std::map<std::string, VectorXf > indi2covar;
		for(int i=0; i<num_traits; ++i)
			indi2pheno_weight[i] = read_phenotype_file(phen_file, trait_names[i], error_variance_weight_names[i]);
		indi2covar = read_covariate_file(covar_file, covar_names);
	
		run_reml_minque(indi2pheno_weight, indi2covar, trait_names, covar_names, 
			grm_list_file, binary_grm_file, keep_file, output_file, h2, minq_iter, minq_tol, 
			n_rand_vec, n_grms, rdist, seed, save_mem, verbose);
		
		return;
	}
	
	// make GRM
	if(make_grm) {
		if(str_geno_coding.empty()) {
			if(grm_type=='A') {
				std::cout<<"\nStarted to make additive GRM.\n";
				std::cout<<"GRM = ZWZ', where Z = standardized genotypes and W = weights.\n"<<std::endl;
			} else {
				std::cout<<"\nStarted to make dominance GRM.\n";
				std::cout<<"GRM = ZWZ', where Z = standardized genotypes for dominance deviation.\n"<<std::endl;
			}
		} else {
			std::cout<<"\nStarted to make GRM with custom genotype coding.\n";
			std::cout<<"GRM = ZWZ', where Z = custom genotype codes and W = weights.\n"<<std::endl;
		}
		
		if( !(file_check(binary_genotype_file + ".fam") && file_check(binary_genotype_file + ".bim") && file_check(binary_genotype_file + ".bed")) )
			throw("\nError: PLINK bed/bim/fam files [" + binary_genotype_file + "] not found\n");
		
		if(!keep_file.empty()) 
			std::cout<<"Set the subject set file ["<<keep_file<<"]\n";
		
		// get subjects from --keep (if any) and PLINK fam
		std::vector<bool> bindi;
		std::vector<std::string> indi_keep;
		get_subject_set(keep_file, binary_genotype_file+".fam", bindi, indi_keep);
		int keep_num = indi_keep.size();
		std::cout<<bindi.size()<<" individuals are found in ["<<binary_genotype_file<<".fam].\n";
		std::cout<<keep_num<<" individuals are retained for making GRM."<<std::endl;
		if(keep_num == 0) throw("\nError: subject set retained for analysis is empty.\n");
		
		// get SNPs from --snp_info_file and PLINK bim
		// read SNP info file
		std::vector<bool> bmarker;
		std::vector<double> gvec;
		std::vector<Vector3d> code_vec;
		std::map<std::string, float> marker2weight;
		get_marker_weight(marker_info_file, marker_effvar_weight_name, marker2weight);
		if(str_geno_coding.empty()) {
			get_marker_set_by_weight(marker2weight, binary_genotype_file+".bim", bmarker, gvec);
		} else {
			std::vector<std::string> coding_names = StrFunc::split(str_geno_coding,',');
			std::map<std::string, Vector3d > marker2codes = get_geno_coding(marker_info_file, coding_names);
			get_marker_set_by_codes(marker2weight, marker2codes, binary_genotype_file+".bim", bmarker, gvec, code_vec);
		}
		Eigen::Map<VectorXd> gwt(gvec.data(), gvec.size());
		int pre_qc_marker_num = gwt.size();
		std::cout<<bmarker.size()<<" variants are found in ["<<binary_genotype_file<<".bim].\n";
		std::cout<<pre_qc_marker_num<<" variants are retained before quality control."<<std::endl;
		if(pre_qc_marker_num == 0) throw("\nError: marker set retained for analysis is empty.\n");
		// Map is not used for the safety of memory continuity.
		Matrix3Xd code_mat(3, code_vec.size());
		for (unsigned i = 0; i < code_vec.size(); ++i) {
			code_mat.col(i) = code_vec[i];
		}
		std::vector<Vector3d>().swap(code_vec);
		std::map<std::string, float>().swap(marker2weight);
		
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
		Matrix3Xd code_sm;
		
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
				if(str_geno_coding.empty()) {
					calc_grm_by_subset(grm_type, min_maf, min_hwep, SS, kmat, hwep, swt, grm, sumwt, post_qc_marker_num);
				} else {
					code_sm = code_mat.middleCols(snp_indx-subset_size,subset_size);
					calc_custom_grm_by_subset(min_maf, min_hwep, SS, kmat, hwep, swt, code_sm, grm, sumwt, post_qc_marker_num);
				}
				
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
			if(str_geno_coding.empty()) {
				calc_grm_by_subset(grm_type, min_maf, min_hwep, SS, kmat, hwep, swt, grm, sumwt, post_qc_marker_num);
			} else {
				code_sm = code_mat.rightCols(rest_size);
				calc_custom_grm_by_subset(min_maf, min_hwep, SS, kmat, hwep, swt, code_sm, grm, sumwt, post_qc_marker_num);
			}
			
			std::cout<<"Completed subset "<<num_sets+1<<std::endl; 
		}
		SS.resize(0,0);
		kmat.resize(0,0);
		
		if(str_geno_coding.empty()) {
			std::cout<<"\nCompleted making "<<(grm_type=='A' ? "additive" : "dominance")<<" GRM.\n";
		} else {
			std::cout<<"\nCompleted making GRM with custom genotype coding.\n";
		}
		std::cout<<post_qc_marker_num<<" post-QC variants are used for making GRM.\n";
		
		std::cout<<"\nWriting GRM to file..."<<std::endl;
		write_dGRM_into_file(output_file + (str_geno_coding.empty() ? (grm_type == 'A' ? "" : ".d") : ".custom"), indi_keep, (float)sumwt, grm);
		std::cout<<"GRM written to file."<<std::endl;
		
		return;
	}
}

void print_help() {
	printf("\nPlease see online documentation at https://jiang18.github.io/mph/\n");
}
