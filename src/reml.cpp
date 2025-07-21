#include "reml.h"
#include <random>
#include <sys/stat.h>

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
	bool verbose)
{
	std::cout<<"\nStarted REML/MINQUE."<<std::endl;
	
	int num_traits = trait_names.size();
	
	int i=0, j=0;
	std::cout<<"\n"<<num_traits<< (num_traits > 1 ? " traits are" : " trait is") <<" included in the analysis:\n";
	for(i=0; i<num_traits; ++i) std::cout<<"["<<i<<"] = "<<trait_names[i]<<"\n";
	
	std::vector<std::string> vc_name;
	std::vector<float> vc_init;
	if(!grm_list_file.empty()) {
		// read GRM list file
		std::ifstream ifs;
		std::string line;
		ifs.open(grm_list_file);
		if(!ifs.is_open()) {
			std::cout<<"Cannot open "<<grm_list_file<<std::endl;
			exit(1);
		}
		std::cout<<"\nReading the GRM list from ["+grm_list_file+"]."<<std::endl;
		while (std::getline(ifs, line)) {
			if(line[line.size()-1] == '\r') line.erase(line.size()-1);
			if(line.empty()) continue;
			std::stringstream  ss(line);
			std::string vc_i;
			float wt = 1.0;
			ss >> vc_i >> wt;
			vc_name.push_back(vc_i);
			vc_init.push_back(wt);
		}
		ifs.close();
		
		std::cout<<vc_name.size()<< (vc_name.size() > 1 ? " matrices are" : " matrix is") <<" in the list."<<std::endl;
		if(std::accumulate(vc_init.begin(), vc_init.end(), 0) <= 0) throw("\nError: the sum of initial VCs in --grm_list must be positive.\n");
	}
	if(!binary_grm_file.empty()) {
		if(!vc_name.empty()) {
			std::cout<<"\n["<<binary_grm_file<<"] specified by --binary_grm_file and appended to the list.\n";
			float avg_vc = std::accumulate(vc_init.begin(), vc_init.end(), 0.0f) / vc_init.size();
			vc_init.push_back(avg_vc);
		} else {
			std::cout<<"\n["<<binary_grm_file<<"] specified by --binary_grm_file as the only GRM.\n";
			vc_init.push_back(1.0);
		}
		vc_name.push_back(binary_grm_file);
	}
	
	// Check if all GRM files exist
	for(const auto& grm_name : vc_name) {
		if(!file_check(grm_name + ".grm.bin") || !file_check(grm_name + ".grm.iid")) {
			throw("\nError: GRM files [" + grm_name + " (.grm.bin/.grm.iid)] not found.\n");
		}
	}
	// Check for duplicates using file comparison
	for (i=0; i<vc_name.size()-1; ++i) {
		for (j=i+1; j<vc_name.size(); ++j) {
			struct stat stat1, stat2;
			stat((vc_name[i] + ".grm.bin").c_str(), &stat1);  // No need to check == 0
			stat((vc_name[j] + ".grm.bin").c_str(), &stat2);  // Files guaranteed to exist
			
			if(stat1.st_dev == stat2.st_dev && stat1.st_ino == stat2.st_ino) {
				throw("\nError: [" + vc_name[i] + "] and [" + vc_name[j] + "] refer to the same GRM.\n");
			}
		}
	}
	if(n_grms == 0) n_grms = vc_name.size();
	if(n_grms > vc_name.size()) {
		throw("\nError: --num_grms exceeds the number of available GRMs (" + std::to_string(vc_name.size()) + ").\n");
	}
	
	std::string indi_file;
	if(keep_file.empty()) {
		indi_file = vc_name[0] + ".grm.iid";
	} else {
		indi_file = keep_file;
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
			yvec[i*keep_num+j] = indi2pheno_weight[i].at(indi_keep[j]).first;
			rvec[i*keep_num+j] = indi2pheno_weight[i].at(indi_keep[j]).second;
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
	for(j=0; j<keep_num; ++j) xmat.row(j) = indi2covar.at(indi_keep[j]);
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
