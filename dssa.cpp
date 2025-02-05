#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector c_propensity_n(NumericVector n, NumericMatrix reactant_matrix){
	//n: should be a vector, not a 1-column matrix!
	NumericVector p(reactant_matrix.ncol());
	for(int j=0; j<reactant_matrix.ncol(); j++){
		p(j) = sum(lfactorial(n - reactant_matrix(_,j)));
	}
	return exp(sum(lfactorial(n)) - p);
}

// [[Rcpp::export]]
NumericVector c_fun_fr(NumericVector kn, NumericVector n, NumericMatrix reactant_matrix){
	return(kn*c_propensity_n(n,reactant_matrix));
}

// [[Rcpp::export]]
NumericVector c_k(Function k, NumericVector n){
	NumericVector kn = as<NumericVector>(k(n));
	return kn;
}

// [[Rcpp::export]]
int find_interval(double x, NumericVector v){
	//v: sorted non-decreasing
	std::vector vv = as<std::vector<double> >(v);
	auto loc = std::lower_bound(vv.begin(), vv.end(), x);
	int index = std::distance(vv.begin(), loc);
	return index; //0-based position index
}

// [[Rcpp::export]]
IntegerVector c_index_in_vector(IntegerVector a, int b){
    IntegerVector indices;
    for (int i = 0; i < a.size(); i++) {
        if (a[i] == b) {
            indices.push_back(i); 
		}
    }
	return indices;
}

// [[Rcpp::export]]
IntegerMatrix c_get_delay_effect_matrix(NumericMatrix S_matrix, NumericMatrix S_matrix_delay){
	IntegerMatrix index(2,0);
	IntegerVector tmp_vec;
	for(int j=0; j<S_matrix_delay.ncol(); j++){
		LogicalVector is_nonzero = (S_matrix_delay(_,j) != 0);
		if(sum(is_nonzero) > 0){
			for(int k=0; k<S_matrix.ncol(); k++){
				if(sum(S_matrix(_,k) == S_matrix_delay(_,j)) == S_matrix.nrow()){
					tmp_vec = NumericVector::create(k+1,j+1);
				}
			}
		}
	}
	return index;
}


// [[Rcpp::export]]
List Dmnr(int iterInit, int iterInc, double tmax, NumericVector n_initial, double t_initial, NumericMatrix S_matrix, NumericMatrix S_matrix_delay, Function k, IntegerVector delay_type, NumericVector delaytime_list, NumericMatrix reactant_matrix_delay, NumericMatrix reactant_matrix, IntegerMatrix delay_effect_matrix ){
	//iterInit: initial memory size for n_values/t_values, e.g. 1000; iterInc: if the initial memory size used up, such a size will be added each time, e.g. 1000
	NumericVector n0 = clone(n_initial);
	
	NumericMatrix n_values(n0.length(), iterInit+1); //n_initial
	NumericVector t_values(iterInit+1); //+1 because initial value also need 1 position
	n_values(_,0) = n0; //n_initial
	t_values(0) = t_initial;
	NumericVector n = n0; //n_initial
	double t = t_initial;
	double min_1,min_2,tau,add_tau;
	int r_1,r_2,r;
	NumericVector tau_vec;
	IntegerVector effect_r;
	int drop_index,index;
	IntegerVector drop_index_vec0;
	double u1;
	
	if(delay_effect_matrix.nrow() != 2){
		delay_effect_matrix = c_get_delay_effect_matrix(S_matrix, S_matrix_delay);
	}

	NumericVector t_vec(S_matrix.ncol());
	NumericVector Tstruct(0);
	IntegerVector Rstruct(0);

	NumericVector kn = c_k(k,n);
	NumericVector f_r = c_fun_fr(kn,n,reactant_matrix);
	NumericVector u2 = runif(S_matrix.ncol());
	NumericVector p_vec = log(1/u2);

	int i = 0;
	while(t < tmax){
		if(i >= t_values.length()-1){ //initial t already hold 1 position, so -1
			NumericMatrix n_values_old = n_values;
			if(FALSE){
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc, n_values_old.begin());
			}else{
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc);
				for(int j=0; j<n_values_old.ncol(); j++){
					n_values(_,j) = n_values_old(_,j);
				}
			}
			NumericVector t_values_old = t_values;
			t_values = NumericVector(t_values_old.length() + iterInc);
			for(int j=0; j<t_values_old.length(); j++){
				t_values[j] = t_values_old[j];
			}
		}
		tau_vec = (p_vec - t_vec)/f_r + t; //time waiting for each reaction
		
		// for (int j = 0; j < tau_vec.length(); j++) {
		// 	if (R_IsNA(tau_vec[j]) || R_finite(tau_vec[j]) == 0) {
		// 		tau_vec[j] = std::numeric_limits<double>::infinity();
		// 	}
		// }

		r_1 = which_min(tau_vec); //next reaction
		min_1 = tau_vec[r_1]; //time of next reaction
		if(Tstruct.length() > 0){
			r_2 = Rstruct(0); //next delay reaction
			min_2 = Tstruct(0); //time of next delay reaction
		}else{
			min_2 = std::numeric_limits<double>::infinity();
		}

		if((i == 0) || (Tstruct.length() == 0) || (min_1 < min_2)){ 
		//if next reaction is non-delay
			r = r_1;
			tau = min_1 - t; //time waiting
			t = min_1;
			if(delay_type(r) == 0){ 
				//reaction r has no delay
				n = n + S_matrix(_,r);
				LogicalVector has_r = delay_effect_matrix(0,_) == r + 1;
				if (sum(has_r) > 0) { 
					IntegerVector effect_r = delay_effect_matrix(1,_);
					effect_r = effect_r[has_r];
					effect_r = effect_r - 1;
					for (int ii = 0; ii < sum(has_r); ++ii) {
						IntegerVector drop_index_vec0 = c_index_in_vector(Rstruct,effect_r(ii));
						if(drop_index_vec0.length()>0){
							drop_index = sample(drop_index_vec0,1)[0];
							if(delay_type(Rstruct(drop_index))==1){
								n = n + reactant_matrix_delay(_,Rstruct(drop_index));	
							}
							Tstruct.erase(drop_index);
							Rstruct.erase(drop_index);

						}
					}
				}
			}else if((delay_type(r) == 1) || (delay_type(r) == 2)){
				//reaction r has delay
				n = n + S_matrix(_,r);
				if(delay_type(r) == 1){ //reaction r has delay type 1 (ICD)
					n = n - reactant_matrix_delay(_,r);
				}else{ //reaction r has delay type 2 (CD)
				}
				add_tau = delaytime_list(r) + t;
				index = find_interval(add_tau, Tstruct);
				Tstruct.insert(index, add_tau);
				Rstruct.insert(index, r);
			}else{
				stop("Error: delay_type is a vector of 0/1/2 only!");
			}
			u1 = runif(1)[0];
			p_vec(r) = p_vec(r) + log(1/u1);
		}else{
		//if next reaction is delay
			r = r_2;
			tau = min_2 - t;
			t = min_2;
			if(delay_type(r) == 0){
				//warning
			}else if(delay_type(r) == 1){
				n = n + S_matrix_delay(_,r) + reactant_matrix_delay(_,r);
			}else if(delay_type(r) == 2){
				n = n + S_matrix_delay(_,r);
			}else{
				stop("Error: delay_type is a vector of 0/1/2 only!");
			}
			Tstruct.erase(0);
			Rstruct.erase(0);
		}

		t_vec = t_vec + f_r * tau;
		kn = c_k(k,n);
		f_r = c_fun_fr(kn, n, reactant_matrix);
		if(t < tail(t_values,1)[0]){
			break;
		}
		t_values[i+1] = t;
		n_values(_,i+1) = n;

		i++;
	}
    return List::create(Named("t_values") = t_values, Named("n_values") = n_values);
}	

// [[Rcpp::export]]
List Drejection(int iterInit, int iterInc, double tmax, NumericVector n_initial, double t_initial, NumericMatrix S_matrix, NumericMatrix S_matrix_delay, Function k, IntegerVector delay_type, NumericVector delaytime_list, NumericMatrix reactant_matrix_delay, NumericMatrix reactant_matrix, IntegerMatrix delay_effect_matrix){
	//iterInit: initial memory size for n_values/t_values, e.g. 1000; iterInc: if the initial memory size used up, such a size will be added each time, e.g. 1000
	NumericVector n0 = clone(n_initial);
	
	NumericMatrix n_values(n0.length(), iterInit+1); //n_initial
	NumericVector t_values(iterInit+1); //+1 because initial value also need 1 position
	n_values(_,0) = n0; //n_initial
	t_values(0) = t_initial;
	NumericVector n = n0; //n_initial
	double t = t_initial;
	double min_2,tau,add_tau;
	int r_2,r;
	// NumericVector tau_vec;
	NumericVector tau_vec;
	IntegerVector effect_r;
	int drop_index,index;
	IntegerVector drop_index_vec0;
	double u1,u2;
	NumericVector kn;
	NumericVector f_r;
	if(delay_effect_matrix.nrow() != 2){
		delay_effect_matrix = c_get_delay_effect_matrix(S_matrix, S_matrix_delay);
	}
	NumericVector Tstruct(0);
	IntegerVector Rstruct(0);
	double lambda_sum;

	int i = 0;
	while(t < tmax){
		if(i >= t_values.length()-1){ //initial t already hold 1 position, so -1
			NumericMatrix n_values_old = n_values;
			if(FALSE){
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc, n_values_old.begin());
			}else{
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc);
				for(int j = 0; j < n_values_old.ncol(); j++){
					n_values(_,j) = n_values_old(_,j);
				}
			}
			NumericVector t_values_old = t_values;
			t_values = NumericVector(t_values_old.length() + iterInc);
			for(int j = 0; j < t_values_old.length(); j++){
				t_values[j] = t_values_old[j];
			}
		}
		u1 = runif(1)[0];
		kn = c_k(k,n);
		f_r = c_fun_fr(kn,n,reactant_matrix);
		lambda_sum = sum(f_r);
		tau = -log(u1) / lambda_sum; //time waiting for next reaction
		if(Tstruct.length() > 0){
			r_2 = Rstruct(0); //next delay reaction
			min_2 = Tstruct(0); //time of next delay reaction
		}else{
			min_2 = std::numeric_limits<double>::infinity();
		}
		if((i == 0) || (Tstruct.length() == 0) || (t+tau < min_2)){
			u2 = runif(1)[0];
			r = find_interval(u2 * lambda_sum, cumsum(f_r));
			t = t + tau;
			if(delay_type(r)==0){ 
				//reaction r has no delay
				n = n + S_matrix(_,r);
				LogicalVector has_r = delay_effect_matrix(0,_) == r + 1;
				if (sum(has_r) > 0) { 
					IntegerVector effect_r = delay_effect_matrix(1,_);
					effect_r = effect_r[has_r];
					effect_r = effect_r - 1;
					for (int ii = 0; ii < sum(has_r); ++ii) {
						IntegerVector drop_index_vec0 = c_index_in_vector(Rstruct,effect_r(ii));
						if(drop_index_vec0.length()>0){
							drop_index = sample(drop_index_vec0,1)[0];
							if(delay_type(Rstruct(drop_index))==1){
								n <- n + reactant_matrix_delay(_,Rstruct(drop_index));	
							}
							Tstruct.erase(drop_index);
							Rstruct.erase(drop_index);

						}
					}
				}
			}else if((delay_type(r) == 1) || (delay_type(r) == 2)){
				//reaction r has delay
				n = n + S_matrix(_,r);
				if(delay_type(r) == 1){ //reaction r has delay type 1 (ICD)
					n = n - reactant_matrix_delay(_,r);
				}else{ //reaction r has delay type 2 (CD)
				}
				add_tau = delaytime_list(r) + t;
				index = find_interval(add_tau, Tstruct);
				Tstruct.insert(index, add_tau);
				Rstruct.insert(index, r);
			}else{
				stop("Error: delay_type is a vector of 0/1/2 only!");
			}
		}else{
			r = r_2;
			tau = min_2 - t;
			t = min_2;
			if(delay_type(r) ==0 ){
				//warning
			}else if(delay_type(r) == 1){
				n = n + S_matrix_delay(_,r) + reactant_matrix_delay(_,r);
			}else if(delay_type(r) == 2){
				n = n + S_matrix_delay(_,r);
			}else{
				stop("Error: delay_type is a vector of 0/1/2 only!");
			}
			Tstruct.erase(0);
			Rstruct.erase(0);
		}

		if(t < tail(t_values,1)[0]){
			break;
		}
		t_values[i+1] = t;
		n_values(_,i+1) = n;

		i++;
	}
    return List::create(Named("t_values") = t_values, Named("n_values") = n_values);
}	

// [[Rcpp::export]]
List Ddirect(int iterInit, int iterInc, double tmax, NumericVector n_initial, double t_initial, NumericMatrix S_matrix, NumericMatrix S_matrix_delay, Function k, IntegerVector delay_type, NumericVector delaytime_list, NumericMatrix reactant_matrix_delay, NumericMatrix reactant_matrix, IntegerMatrix delay_effect_matrix){
	//iterInit: initial memory size for n_values/t_values, e.g. 1000; iterInc: if the initial memory size used up, such a size will be added each time, e.g. 1000
	NumericVector n0 = clone(n_initial);
	
	NumericMatrix n_values(n0.length(), iterInit+1); //n_initial
	NumericVector t_values(iterInit+1); //+1 because initial value also need 1 position
	n_values(_,0) = n0; //n_initial
	t_values(0) = t_initial;
	NumericVector n = n0; //n_initial
	double t = t_initial;
	double min_2,tau,add_tau;
	int r_2,r;
	// NumericVector tau_vec;
	NumericVector tau_vec;
	IntegerVector effect_r;
	int drop_index,index;
	IntegerVector drop_index_vec0;
	double u1,u2;
	if(delay_effect_matrix.nrow() != 2){
		delay_effect_matrix = c_get_delay_effect_matrix(S_matrix, S_matrix_delay);
	}

	NumericVector Tstruct(0);
	IntegerVector Rstruct(0);
	double lambda_sum;
	NumericVector kn;
	NumericVector f_r;

	int i = 0;
	while(t < tmax){
		if(i >= t_values.length()-1){ //initial t already hold 1 position, so -1
			NumericMatrix n_values_old = n_values;
			if(FALSE){
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc, n_values_old.begin());
			}else{
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc);
				for(int j=0; j<n_values_old.ncol(); j++){
					n_values(_,j) = n_values_old(_,j);
				}
			}
			NumericVector t_values_old = t_values;
			t_values = NumericVector(t_values_old.length() + iterInc);
			for(int j = 0; j < t_values_old.length(); j++){
				t_values[j] = t_values_old[j];
			}
		}
		u1 = runif(1)[0];
		if(Tstruct.length() > 0){
			r_2 = Rstruct(0); //next delay reaction
			min_2 = Tstruct(0); //time of next delay reaction
		}else{
			min_2 = std::numeric_limits<double>::infinity();
		}

		if((i == 0) || (Tstruct.length() == 0)){
			kn = c_k(k,n);
			f_r = c_fun_fr(kn,n,reactant_matrix);
			lambda_sum = sum(f_r);
			tau = -log(u1) / lambda_sum; 
		}else{
			int j = 0;
			kn = c_k(k,n);
			f_r = c_fun_fr(kn,n,reactant_matrix);
			lambda_sum = sum(f_r);
			double a0 = 0;
			double at = lambda_sum * (min_2 - t);
			double F = 1 - exp(- at);
			while (F < u1) {
				r = Rstruct(j);
				if(delay_type(r)==0){
					//warning
				}else if(delay_type(r) == 1){
					n = n + S_matrix_delay(_,r) + reactant_matrix_delay(_,r);
				}else if(delay_type(r) == 2){
					n = n + S_matrix_delay(_,r);
				}else{
					stop("Error: delay_type is a vector of 0/1/2 only!");
				}
				kn = c_k(k,n);
				f_r = c_fun_fr(kn,n,reactant_matrix);	
				lambda_sum = sum(f_r);
				if(j == (Tstruct.length() - 1)){
					a0 = at;	
					F = 1;
				}else{
					double add = lambda_sum * (Tstruct(j + 1) - Tstruct(j));
					a0 = at;	
					at = at + add;
					F = 1 - exp(- at);
				}
				t = Tstruct[j];
				j = j + 1;

				t_values[i+1] = t;
				n_values(_,i+1) = n;
				i++;
				if(i >= t_values.length()-1){ //initial t already hold 1 position, so -1
					NumericMatrix n_values_old = n_values;
					if(FALSE){
						n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc, n_values_old.begin());
					}else{
						n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc);
						for(int j=0; j<n_values_old.ncol(); j++){
							n_values(_,j) = n_values_old(_,j);
						}
					}
					NumericVector t_values_old = t_values;
					t_values = NumericVector(t_values_old.length() + iterInc);
					for(int j = 0; j < t_values_old.length(); j++){
						t_values[j] = t_values_old[j];
					}
				}
			}
			// j = j - 1;
			tau = - (log(1 - u1) + a0) / lambda_sum;
			Tstruct.erase(Tstruct.begin(), Tstruct.begin() + j);
			Rstruct.erase(Rstruct.begin(), Rstruct.begin() + j);
		}
		u2 = runif(1)[0];
		r = find_interval(u2 * lambda_sum, cumsum(f_r));
		t = t + tau;
		if(delay_type(r) == 0){ 
			//reaction r has no delay
			n = n + S_matrix(_,r);
			LogicalVector has_r = delay_effect_matrix(0,_) == r + 1;
			if (sum(has_r) > 0) { 
				IntegerVector effect_r = delay_effect_matrix(1,_);
				effect_r = effect_r[has_r];
				effect_r = effect_r - 1;
				for (int ii = 0; ii < sum(has_r); ++ii) {
					IntegerVector drop_index_vec0 = c_index_in_vector(Rstruct,effect_r(ii));
					if(drop_index_vec0.length()>0){
						drop_index = sample(drop_index_vec0,1)[0];
						if(delay_type(Rstruct(drop_index))==1){
							n <- n + reactant_matrix_delay(_,Rstruct(drop_index));	
						}
						Tstruct.erase(drop_index);
						Rstruct.erase(drop_index);

					}
				}
			}
		}else if((delay_type(r) == 1) || (delay_type(r) == 2)){
			//reaction r has delay
			n = n + S_matrix(_,r);
			if(delay_type(r) == 1){ //reaction r has delay type 1 (ICD)
				n = n - reactant_matrix_delay(_,r);
			}else{ //reaction r has delay type 2 (CD)
			}
			add_tau = delaytime_list(r) + t;
			index = find_interval(add_tau, Tstruct);
			Tstruct.insert(index, add_tau);
			Rstruct.insert(index, r);
		}else{
			stop("Error: delay_type is a vector of 0/1/2 only!");
		}

		if(t < tail(t_values,1)[0]){
			break;
		}
		t_values[i+1] = t;
		n_values(_,i+1) = n;
		i++;
	}
    return List::create(Named("t_values") = t_values, Named("n_values") = n_values);
}	

// [[Rcpp::export]]
List direct(int iterInit, int iterInc, double tmax, NumericVector n_initial, double t_initial, NumericMatrix S_matrix, Function k, NumericMatrix reactant_matrix){
	//iterInit: initial memory size for n_values/t_values, e.g. 1000; iterInc: if the initial memory size used up, such a size will be added each time, e.g. 1000
	NumericVector n0 = clone(n_initial);
	
	NumericMatrix n_values(n0.length(), iterInit+1); //n_initial
	NumericVector t_values(iterInit+1); //+1 because initial value also need 1 position
	n_values(_,0) = n0; //n_initial
	t_values(0) = t_initial;
	NumericVector n = n0; //n_initial
	double t = t_initial;
	double tau;
	int r;
	double u1,u2;
	double lambda_sum;
	NumericVector kn;
	NumericVector f_r;

	int i = 0;
	while(t < tmax){
		if(i >= t_values.length()-1){ //initial t already hold 1 position, so -1
			NumericMatrix n_values_old = n_values;
			if(FALSE){
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc, n_values_old.begin());
			}else{
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc);
				for(int j=0; j<n_values_old.ncol(); j++){
					n_values(_,j) = n_values_old(_,j);
				}
			}
			NumericVector t_values_old = t_values;
			t_values = NumericVector(t_values_old.length() + iterInc);
			for(int j = 0; j < t_values_old.length(); j++){
				t_values[j] = t_values_old[j];
			}
		}
		u1 = runif(1)[0];
		u2 = runif(1)[0];
		kn = c_k(k,n);
		f_r = c_fun_fr(kn,n,reactant_matrix);
		lambda_sum = sum(f_r);
		tau = - log(u1) / lambda_sum;
		r = find_interval(u2 * lambda_sum, cumsum(f_r));
		n = n + S_matrix(_,r);
		t = t + tau;
		if(t < tail(t_values,1)[0]){
			break;
		}
		t_values[i+1] = t;
		n_values(_,i+1) = n;
		i++;
	}
    return List::create(Named("t_values") = t_values, Named("n_values") = n_values);
}	

// [[Rcpp::export]]
List mnr(int iterInit, int iterInc, double tmax, NumericVector n_initial, double t_initial, NumericMatrix S_matrix, Function k, NumericMatrix reactant_matrix){
	//iterInit: initial memory size for n_values/t_values, e.g. 1000; iterInc: if the initial memory size used up, such a size will be added each time, e.g. 1000
	NumericVector n0 = clone(n_initial);
	
	NumericMatrix n_values(n0.length(), iterInit+1); //n_initial
	NumericVector t_values(iterInit+1); //+1 because initial value also need 1 position
	n_values(_,0) = n0; //n_initial
	t_values(0) = t_initial;
	NumericVector n = n0; //n_initial
	double t = t_initial;
	double tau;
	int r;
	double u1;
	double lambda_sum;
	NumericVector tau_vec;
	NumericVector t_vec(S_matrix.ncol());	
	NumericVector u2 = runif(S_matrix.ncol());
	NumericVector p_vec = log(1/u2);
	NumericVector kn = c_k(k,n);
	NumericVector f_r = c_fun_fr(kn,n,reactant_matrix);

	int i = 0;
	while(t < tmax){
		if(i >= t_values.length()-1){ //initial t already hold 1 position, so -1
			NumericMatrix n_values_old = n_values;
			if(FALSE){
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc, n_values_old.begin());
			}else{
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc);
				for(int j=0; j<n_values_old.ncol(); j++){
					n_values(_,j) = n_values_old(_,j);
				}
			}
			NumericVector t_values_old = t_values;
			t_values = NumericVector(t_values_old.length() + iterInc);
			for(int j = 0; j < t_values_old.length(); j++){
				t_values[j] = t_values_old[j];
			}
		}
		tau_vec = (p_vec - t_vec)/f_r;
		r = which_min(tau_vec);
		tau = tau_vec(r);
		n = n + S_matrix(_,r);
		t_vec = t_vec + f_r * tau;
		u1 = runif(1)[0];
		p_vec(r) = p_vec(r) + log(1/u1);
		kn = c_k(k,n);
		f_r = c_fun_fr(kn,n,reactant_matrix);

		t = t + tau;
		if(t < tail(t_values,1)[0]){
			break;
		}
		t_values[i+1] = t;
		n_values(_,i+1) = n;
		i++;
	}
    return List::create(Named("t_values") = t_values, Named("n_values") = n_values);
}	

// [[Rcpp::export]]
List nr(int iterInit, int iterInc, double tmax, NumericVector n_initial, double t_initial, NumericMatrix S_matrix, Function k, NumericMatrix reactant_matrix){
	//iterInit: initial memory size for n_values/t_values, e.g. 1000; iterInc: if the initial memory size used up, such a size will be added each time, e.g. 1000
	NumericVector n0 = clone(n_initial);
	
	NumericMatrix n_values(n0.length(), iterInit+1); //n_initial
	NumericVector t_values(iterInit+1); //+1 because initial value also need 1 position
	n_values(_,0) = n0; //n_initial
	t_values(0) = t_initial;
	NumericVector n = n0; //n_initial
	double t = t_initial;
	double tau;
	int r;
	double u1;
	double lambda_sum;
	NumericVector t_vec(S_matrix.ncol());	
	NumericVector u2 = runif(S_matrix.ncol());
	NumericVector kn = c_k(k,n);
	NumericVector f_r = c_fun_fr(kn,n,reactant_matrix);
	NumericVector f_r_update;
	NumericVector tau_vec = -log(u2)/f_r;

	int i = 0;
	while(t < tmax){
		if(i >= t_values.length()-1){ //initial t already hold 1 position, so -1
			NumericMatrix n_values_old = n_values;
			if(FALSE){
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc, n_values_old.begin());
			}else{
				n_values = NumericMatrix(n0.length(), n_values_old.ncol()+iterInc);
				for(int j=0; j<n_values_old.ncol(); j++){
					n_values(_,j) = n_values_old(_,j);
				}
			}
			NumericVector t_values_old = t_values;
			t_values = NumericVector(t_values_old.length() + iterInc);
			for(int j = 0; j < t_values_old.length(); j++){
				t_values[j] = t_values_old[j];
			}
		}
		r = which_min(tau_vec);
		tau = tau_vec(r);
		n = n + S_matrix(_,r);
		kn = c_k(k,n);
		f_r_update = c_fun_fr(kn,n,reactant_matrix);

		for (int j = 0; j < tau_vec.length(); j++) {
			if (j != r) {  
				tau_vec[j] = f_r[j] / f_r_update[j] * (tau_vec[j] - tau) + tau;
			}
		}
		for (int j = 0; j < tau_vec.length(); j++) {
			if (R_IsNA(tau_vec[j]) || R_finite(tau_vec[j]) == 0) {
				tau_vec[j] = -log(runif(1)[0]) / f_r_update[j] + tau;
			}
		}

		u1 = runif(1)[0];
		tau_vec(r) = 1/f_r_update(r)*log(1/u1)+tau;
		f_r = f_r_update;

		t = tau;
		if(t < tail(t_values,1)[0]){
			break;
		}
		t_values[i+1] = t;
		n_values(_,i+1) = n;
		i++;
	}
    return List::create(Named("t_values") = t_values, Named("n_values") = n_values);
}	
