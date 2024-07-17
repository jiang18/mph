#include "mph.h"
#include <fcntl.h>
#include <unistd.h>

void calc_grm_by_subset(const char grm_type, const float min_maf, const float min_hwep, 
	Ref<MatrixXd> SS, const Ref<MatrixXc>& kmat, const Ref<VectorXd>& hwep, Ref<VectorXd> swt, 
	Ref<MatrixXd> grm, double &sumwt, int &post_qc_marker_num)
{
	if(grm_type == 'A') {
		#pragma omp parallel for
		for(int s=0; s<kmat.cols(); ++s) {
			SS.col(s) = kmat.col(s).cast<double>();
			double pfreq = SS.col(s).mean()/2.;
			if(hwep(s)<min_hwep || (float)pfreq>=1-min_maf || (float)pfreq<=min_maf || swt(s)<=0) {
				swt(s) = 0;
				SS.col(s).setZero();
			} else {
				SS.col(s).array() -= 2*pfreq;
				SS.col(s) *= std::sqrt(swt(s)/(2*pfreq*(1-pfreq)));
			}
		}
		// std::cout<<"GRM += ZWZ', where Z=subset(standardized genotypes)."<<std::endl;
		grm.selfadjointView<Lower>().rankUpdate(SS);
		post_qc_marker_num += (swt.array()!=0).count();
		sumwt += swt.sum();
	} else if(grm_type == 'D') {
		#pragma omp parallel for
		for(int s=0; s<kmat.cols(); ++s) {
			double pfreq = kmat.col(s).cast<double>().mean()/2.;
			double npq_quot = pfreq/(pfreq - 1.0);
			double nqp_quot = (pfreq - 1.0)/pfreq;
			
			if(hwep(s)<min_hwep || (float)pfreq>=1-min_maf || (float)pfreq<=min_maf || swt(s)<=0) {
				swt(s) = 0;
				SS.col(s).setZero();
			} else {
				for(int l=0; l<kmat.rows(); ++l) {
					if(kmat(l,s) == 0) SS(l,s) = npq_quot;
					else if(kmat(l,s) == 1) SS(l,s) = 1.0;
					else SS(l,s) = nqp_quot;
				}
				if(swt(s)>1.000001 || swt(s)<0.999999) SS.col(s) *= std::sqrt(swt(s));
			}
		}
		// std::cout<<"GRM += ZWZ', where Z=subset(standardized genotypes for dominance deviation)."<<std::endl;
		grm.selfadjointView<Lower>().rankUpdate(SS);
		post_qc_marker_num += (swt.array()!=0).count();
		sumwt += swt.sum();
	}
}

void calc_custom_grm_by_subset(const float min_maf, const float min_hwep, 
	Ref<MatrixXd> SS, const Ref<MatrixXc>& kmat, const Ref<VectorXd>& hwep, Ref<VectorXd> swt, const Ref<Matrix3Xd>& code_sm, 
	Ref<MatrixXd> grm, double &sumwt, int &post_qc_marker_num)
{
	#pragma omp parallel for
	for(int s=0; s<kmat.cols(); ++s) {
		double pfreq = kmat.col(s).cast<double>().mean()/2.;
		if(hwep(s)<min_hwep || (float)pfreq>=1-min_maf || (float)pfreq<=min_maf || swt(s)<=0) {
			swt(s) = 0;
			SS.col(s).setZero();
		} else {
			for(int l=0; l<kmat.rows(); ++l) {
				if(kmat(l,s) == 0) SS(l,s) = code_sm(0,s);
				else if(kmat(l,s) == 1) SS(l,s) = code_sm(1,s);
				else SS(l,s) = code_sm(2,s);
			}
			if(swt(s)>1.000001 || swt(s)<0.999999) SS.col(s) *= std::sqrt(swt(s));
		}
	}
	// std::cout<<"GRM += ZWZ', where Z=subset(custom genotype codes)."<<std::endl;
	grm.selfadjointView<Lower>().rankUpdate(SS);
	post_qc_marker_num += (swt.array()!=0).count();
	sumwt += swt.sum();
}

void write_grm_into_file(const std::string bin_file_prefix, const std::vector<std::string>& indi_keep, const float sum2pq, const Ref<MatrixXf> grm)
{
	int num = indi_keep.size();
	std::cout<<"There are "<<num<<" individuals in GRM."<<std::endl;

	int i;

	std::string indi_file = bin_file_prefix + ".grm.iid";
	std::ofstream ofs_indi (indi_file, std::ofstream::out);
	for(i=0; i<num; ++i)
		ofs_indi<<indi_keep[i]<<std::endl;
	ofs_indi.close();
	
	std::string bin_file = bin_file_prefix + ".grm.bin";
	std::ofstream ofs_bin (bin_file, std::ofstream::binary);
	ofs_bin.write((char*)(&num), sizeof(int));
	ofs_bin.write((char*)(&sum2pq), sizeof(float));
	for(i=0; i<num; ++i) 
		ofs_bin.write((char*)grm.col(i).tail(num-i).data(), (num-i)*sizeof(float));
	ofs_bin.close();
}

void write_dGRM_into_file(const std::string bin_file_prefix, const std::vector<std::string>& indi_keep, const float sum2pq, const Ref<MatrixXd> grm)
{
	int num = indi_keep.size();
	std::cout<<"There are "<<num<<" individuals in GRM."<<std::endl;

	int i;

	std::string indi_file = bin_file_prefix + ".grm.iid";
	std::ofstream ofs_indi (indi_file, std::ofstream::out);
	for(i=0; i<num; ++i)
		ofs_indi<<indi_keep[i]<<std::endl;
	ofs_indi.close();
	
	std::string bin_file = bin_file_prefix + ".grm.bin";
	std::ofstream ofs_bin (bin_file, std::ofstream::binary);
	ofs_bin.write((char*)(&num), sizeof(int));
	ofs_bin.write((char*)(&sum2pq), sizeof(float));
	VectorXf colvec(num);
	for(i=0; i<num; ++i) {
		colvec = grm.col(i).cast<float>();
		ofs_bin.write((char*)colvec.tail(num-i).data(), (num-i)*sizeof(float));
	}
	ofs_bin.close();
}

int get_grm_dim(const std::string binary_grm_file_prefix)
{
	std::ifstream ifs;
	std::string line;
	int i=0;
	// read individual list file
	std::string indi_file = binary_grm_file_prefix + ".grm.iid";
	ifs.open(indi_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<indi_file<<std::endl;
		exit(1);
	}
	while (std::getline(ifs, line))
		i++;
	ifs.close();
	return(i);
}

void get_grm_indi(const std::string indi_file, std::vector<std::string>& indi_keep)
{
	std::ifstream ifs;
	std::string line;
	// read individual list file
	ifs.open(indi_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<indi_file<<std::endl;
		exit(1);
	}
	indi_keep.clear();
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		indi_keep.push_back(line);
	}
	ifs.close();
}

/*
1) grm's dim is the same as the size of indi_keep.
2) Only the lower triangular part of grm is referenced.
*/
void read_grm_from_file(const std::string binary_grm_file_prefix, const std::vector<std::string>& indi_keep, float& sum2pq, Ref<MatrixXf> grm, bool verbose)
{
	// quick check 1
	if(indi_keep.size() != grm.cols()) 
		throw("\nInconsistency of GRM dimension found in read_grm_from_file().\n");
	std::ifstream ifs;
	std::string line;
	int i=0;
	// read individual list file
	std::string indi_file = binary_grm_file_prefix + ".grm.iid";
	ifs.open(indi_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<indi_file<<std::endl;
		exit(1);
	}
	if(verbose) std::cout<<"Reading GRM IID file from ["+indi_file+"]."<<std::endl;
	
	std::map<std::string, int> indi2idx;
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		indi2idx[line] = i++;
	}
	ifs.close();
	if(verbose) std::cout<<"There are "<<i<<" individuals in GRM file."<<std::endl;
	
	//read binary GRM file
	std::string grm_file = binary_grm_file_prefix + ".grm.bin";
	int pfd = open(grm_file.c_str(), O_RDONLY);
	if(pfd == -1) {
		std::cout<<"Cannot open "<<grm_file<<std::endl;
		exit(1);
	}
	if(verbose) std::cout<<"Reading GRM BIN file from ["+grm_file+"]."<<std::endl;
	int num;
	read(pfd, (void*)(&num), sizeof(int));
	read(pfd, (void*)(&sum2pq), sizeof(float));
	if(num != i) throw("\nInconsistency between GRM IID file and GRM BIN file found regarding population size.\n");

	const int grm_dim = indi_keep.size();
	ArrayXi vidx(grm_dim);
	for(i=0; i<grm_dim; ++i) {
		if(indi2idx.find(indi_keep[i]) == indi2idx.end()) 
			throw("\n" + indi_keep[i] + " found in phenotype/genotype/subject files but not in the GRM file [" + indi_file + "].\nRemove it in phenotype/subject file and retry.\n");
		else vidx(i) = indi2idx.at(indi_keep[i]);
	}
	ArrayXi vidx0 = vidx;
	bool cont_idx = true;
	if(vidx(vidx.size()-1) - vidx(0) != vidx.size() - 1) {
		cont_idx = false;
		std::sort(vidx.data(), vidx.data()+vidx.size());
	}
	
	#pragma omp parallel
	{
		VectorXf buf(num);
		#pragma omp for
		for(int k=0; k<vidx.size(); ++k) {
			int i = vidx(k);
			std::uint64_t offset = static_cast<std::uint64_t>(2*num - i + 1) * i / 2 * sizeof(float) + sizeof(int) + sizeof(float);
			pread64(pfd, (void*)buf.tail(num-i).data(), (num-i)*sizeof(float), offset);
			if(cont_idx) grm.col(k).tail(grm_dim-k) = buf.segment(i,grm_dim-k);
			else grm.col(k).tail(grm_dim-k) = buf(vidx.tail(grm_dim-k));
		}
	}
	close(pfd);
	
	// Do not use selfadjointView to fill upper triangle. It is slow.
	// Better use grm.triangularView<StrictlyUpper>() = grm.transpose();
	
	if(cont_idx) return;
	else if((vidx-vidx0==0).all()) return;
	else {
		// Filling the upper triangular part is only necessary for reordering.
		#pragma omp parallel for
		for(int i=1; i<grm_dim; ++i) {
			grm.col(i).head(i) = grm.row(i).head(i);
		}
		
		std::map<int, int> midx;
		for(int i=0; i<vidx.size(); ++i) 
			midx[vidx(i)] = i;
		ArrayXi pidx(grm_dim);
		for(int i=0; i<vidx0.size(); ++i) 
			pidx(i) = midx.at(vidx0(i));
		
		#pragma omp parallel
		{
			VectorXf buf(grm_dim); // necessary for indexing
			#pragma omp for
			for(int i=0; i<grm_dim; ++i) {
				buf = grm.col(i)(pidx);
				grm.col(i) = buf;
			}
		}
		PermutationMatrix<Dynamic, Dynamic> tr;
		tr.indices() = pidx;
		grm = grm * tr;
	}
}

// TODO: modify GRM reading and GRM triangularView filling following the above function. 
void read_dGRM_from_file(const std::string binary_grm_file_prefix, const std::vector<std::string>& indi_keep, float& sum2pq, Ref<MatrixXd> grm)
{
	std::ifstream ifs;
	std::string line;
	int i=0, j;
	// read individual list file
	std::string indi_file = binary_grm_file_prefix + ".grm.iid";
	ifs.open(indi_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<indi_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading GRM IID file from ["+indi_file+"]."<<std::endl;
	
	std::map<std::string, int> indi2idx;
	std::vector<std::string> grm_iid;
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		indi2idx[line] = i++;
		grm_iid.push_back(line);
	}
	ifs.close();
	std::cout<<"There are "<<i<<" individuals in GRM file."<<std::endl;
	
	//read binary GRM file
	std::string grm_file = binary_grm_file_prefix + ".grm.bin";
	ifs.open(grm_file, std::ios::binary);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<grm_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading GRM BIN file from ["+grm_file+"]."<<std::endl;
	int num;
	ifs.read((char*)(&num), sizeof(int));
	ifs.read((char*)(&sum2pq), sizeof(float));
	if(num != i) throw("\nInconsistency between GRM IID file and GRM BIN file found regarding population size.\n");

	int grm_dim = indi_keep.size();
	for(i=0; i<grm_dim; ++i) {
		if(indi2idx.find(indi_keep[i]) == indi2idx.end()) 
			throw("\n" + indi_keep[i] + " found in phenotype/genotype/subject files but not in the GRM file [" + indi_file + "].\nRemove it in phenotype/subject file and retry.\n");
	}
	
	VectorXf colvec(num);
	for(i=0; i<num; ++i) {
		ifs.read((char*)colvec.tail(num-i).data(), (num-i)*sizeof(float));
		grm.col(i) = colvec.cast<double>();
	}
	ifs.close();
	
	grm.triangularView<StrictlyUpper>() = grm.transpose();

	for(i=0; i<grm_dim; ++i) {
		j = indi2idx[indi_keep[i]];
		if(i == j) continue;
		grm.col(i).swap(grm.col(j));
		grm.row(i).swap(grm.row(j));

		grm_iid[j] = grm_iid[i];
		indi2idx[grm_iid[j]] = j;
	}

}

void merge_grms(const std::string grm_list_file, const std::vector<std::string>& indi_keep, float& out_sum2pq, Ref<MatrixXf> out_grm)
{
	std::ifstream ifs;
	std::string line;
	// read individual list file
	ifs.open(grm_list_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<grm_list_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading the GRM list from ["+grm_list_file+"]."<<std::endl;
	
	out_grm.setZero();
	out_sum2pq = 0;
	
	MatrixXf grm (indi_keep.size(), indi_keep.size());
	float sum2pq;
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::stringstream ss(line);
		std::string vc_i;
		ss >> vc_i;
		read_grm_from_file(vc_i, indi_keep, sum2pq, grm);
		out_grm.triangularView<Lower>() += grm;
		out_sum2pq += sum2pq;
		std::cout<<"GRM ["<<vc_i<<"] has been merged."<<std::endl;
	}
	ifs.close();
}

void deduct_grms(const std::string grm_list_file, const std::vector<std::string>& indi_keep, float& out_sum2pq, Ref<MatrixXf> out_grm)
{
	std::ifstream ifs;
	std::string line;
	// read individual list file
	ifs.open(grm_list_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<grm_list_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading the GRM list from ["+grm_list_file+"]."<<std::endl;
	
	int line_num = 0;
	
	MatrixXf grm (indi_keep.size(), indi_keep.size());
	float sum2pq;
	while (std::getline(ifs, line)) {
		line_num ++;
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::stringstream ss(line);
		std::string vc_i;
		ss >> vc_i;
		read_grm_from_file(vc_i, indi_keep, sum2pq, grm);
		if(line_num == 1) {
			out_grm.triangularView<Lower>() = grm;
			out_sum2pq = sum2pq;
			std::cout<<"GRM ["<<vc_i<<"] has been read."<<std::endl;
		} else {
			out_grm.triangularView<Lower>() -= grm;
			out_sum2pq -= sum2pq;
			std::cout<<"GRM ["<<vc_i<<"] has been deducted."<<std::endl;
			
			if(out_sum2pq <=0) throw("\nError: --deduct_grms results in negative semidefinite GRM.\n");
		}
	}
	ifs.close();
}

void sum_grms(const std::string grm_list_file, const std::vector<std::string>& indi_keep, float& gvar, Ref<MatrixXf> out_grm)
{
	std::ifstream ifs;
	std::string line;
	// read individual list file
	ifs.open(grm_list_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<grm_list_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading the GRM list from ["+grm_list_file+"]."<<std::endl;
	
	gvar = 0;
	out_grm.setZero();
	
	MatrixXf grm (indi_keep.size(), indi_keep.size());
	float sum2pq;
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::stringstream ss(line);
		std::string vc_i;
		float wt = 1.0;
		ss >> vc_i >> wt;
		read_grm_from_file(vc_i, indi_keep, sum2pq, grm);
		out_grm.triangularView<Lower>() += grm * (wt/sum2pq);
		gvar += wt;
		std::cout<<"GRM ["<<vc_i<<"] has been added."<<std::endl;
	}
	ifs.close();
}

void make_core(const std::string grm_list_file, const std::vector<std::string>& indi_keep, const std::string output_file)
{
	std::vector<std::string> vc;
	std::vector<std::string> lab;
	std::ifstream ifs;
	std::string line;
	// read GRM list file
	ifs.open(grm_list_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<grm_list_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading the GRM list from ["+grm_list_file+"]."<<std::endl;
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::stringstream  ss(line);
		std::string vc_i;
		std::string lab_i = "";
		ss >> vc_i >> lab_i;
		vc.push_back(vc_i);
		lab.push_back(lab_i);
	}
	ifs.close();
	for(int i=0; i<lab.size(); ++i) {
		if(lab[i] == "") lab[i] = std::to_string(i+1);
	}
	std::cout<<"There are "<<vc.size()<<" matrices in the list.\n"<<std::endl;
	if(vc.size()<2) throw("\nError: --make_core requires at least 2 VCs in --grm_list.\n");
	
	const int n_packed = vc.size();
	std::vector<MatrixXf> gg(n_packed);
	MatrixXf grm(indi_keep.size(), indi_keep.size());
	float sumwt;
	for(int i=0; i<vc.size(); ++i) {
		line = vc[i];
		read_grm_from_file(line, indi_keep, sumwt, grm);
		std::cout<<"GRM ["<<line<<"] has been obtained.\n";
		
		std::cout<<"Computing Cholesky("<<line<<") = L * U."<<std::endl;
		grm.triangularView<Lower>() /= sumwt;
		grm.diagonal().array() += 1e-4; // to improve PD
		LLT<Ref<MatrixXf> > llt(grm);
		
		gg[i] = llt.matrixL();
		std::cout<<"Matrix L for ["<<line<<"] has been obtained."<<std::endl;
	}
	
	for(int i=0; i<vc.size(); ++i) {
		for(int j=i+1; j<vc.size(); ++j) {
			std::cout<<"\nMaking the CORE relationship matrix for "<<lab[i]<<"x"<<lab[j]<<std::endl;
			grm.noalias() = gg[i] * gg[j].transpose();
			grm.triangularView<StrictlyLower>() += grm.transpose();
			grm.diagonal() *= 2.;
			
			std::string this_output = output_file + ".core." + lab[i] + "x" + lab[j];
			write_grm_into_file(this_output, indi_keep, 1, grm);
			std::cout<<"Written the CORE relationship matrix to "<<this_output<<std::endl;
		}
	}
}

void make_fore(const std::string grm_list_file, const std::vector<std::string>& indi_keep, const std::string output_file)
{
	std::vector<std::string> vc;
	std::vector<std::string> lab;
	std::ifstream ifs;
	std::string line;
	// read GRM list file
	ifs.open(grm_list_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<grm_list_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading the GRM list from ["+grm_list_file+"]."<<std::endl;
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::stringstream  ss(line);
		std::string vc_i;
		std::string lab_i = "";
		ss >> vc_i >> lab_i;
		vc.push_back(vc_i);
		lab.push_back(lab_i);
	}
	ifs.close();
	for(int i=0; i<lab.size(); ++i) {
		if(lab[i] == "") lab[i] = std::to_string(i+1);
	}
	std::cout<<"There are "<<vc.size()<<" matrices in the list.\n"<<std::endl;
	
	const int n_packed = (vc.size() - 1)/2 + 1;
	MatrixXf dd (indi_keep.size(), vc.size());
	std::vector<MatrixXf> gg(n_packed);
	MatrixXf grm(indi_keep.size(), indi_keep.size());
	float sumwt;
	for(int i=0; i<vc.size(); ++i) {
		line = vc[i];
		read_grm_from_file(line, indi_keep, sumwt, grm);
		
		grm.triangularView<Lower>() /= sumwt;
		
		if(i % 2 == 0) gg[i/2] = grm.transpose();
		else gg[i/2].triangularView<Lower>() = grm;
		dd.col(i) = gg[i/2].diagonal();
		
		std::cout<<"GRM ["<<line<<"] has been obtained."<<std::endl;
	}
	
	for(int i=0; i<vc.size(); ++i) {
		for(int j=i; j<vc.size(); ++j) {
			std::cout<<"\nMaking the FORE relationship matrix for "<<lab[i]<<"x"<<lab[j]<<std::endl;
			if((i % 2 == 0) && (j % 2 == 0)) {
				grm.triangularView<StrictlyUpper>() = gg[i/2].cwiseProduct(gg[j/2]);
				grm.triangularView<StrictlyLower>() = grm.transpose();
			} else if((i % 2 == 0) && (j % 2 != 0)) {
				grm.triangularView<StrictlyLower>() = gg[i/2].transpose().cwiseProduct(gg[j/2]);
			} else if((i % 2 != 0) && (j % 2 == 0)) {
				grm.triangularView<StrictlyLower>() = gg[i/2].cwiseProduct(gg[j/2].transpose());
			} else {
				grm.triangularView<StrictlyLower>() = gg[i/2].cwiseProduct(gg[j/2]);
			}
			grm.diagonal() = dd.col(i).array() * dd.col(j).array();
			float scaling = grm.trace()/indi_keep.size();
			grm.triangularView<Lower>() /= scaling;
			
			std::string this_output = output_file + ".fore." + lab[i] + "x" + lab[j];
			write_grm_into_file(this_output, indi_keep, 1, grm);
			std::cout<<"Written the CORE relationship matrix to "<<this_output<<std::endl;
		}
	}
}
