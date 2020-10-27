# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <fstream>

using namespace std;
using namespace arma;

//[[Rcpp::export()]]
arma::rowvec EM(int n1,int n2,int num_snp,int num_per,arma::mat sigma,double alpha,double sigma_beta,
	double sigma_gamma, arma::colvec ldscore, arma::mat zscore) {

	double sqrt_n = sqrt(n1)*sqrt(n2);

	arma::mat Al_Ga(2, 2);
	Al_Ga(0, 0) = n1 * alpha*alpha*sigma_beta + n1 * sigma_gamma;
	Al_Ga(0, 1) = sqrt_n * alpha*sigma_beta;
	Al_Ga(1, 0) = sqrt_n * alpha*sigma_beta;
	Al_Ga(1, 1) = n2 * sigma_beta;

	arma::mat sigma_inv = inv(sigma);

	arma::colvec D_alpha(2);
	D_alpha(0) = alpha * sqrt(n1);
	D_alpha(1) = sqrt(n2);

	arma::colvec D_gamma(2);
	D_gamma(0) = sqrt(n1);
	D_gamma(1) = 0;

	arma::colvec D(2);
	D(0) = sqrt(n1);
	D(1) = 0;

	arma::rowvec par_vec(6);
	par_vec.zeros();

	for (int i = 1; i <= num_snp; i++) {

		double lj = ldscore(i - 1);
		double rj = lj + num_snp / num_per;
		arma::colvec zj(2);
		zj(0) = zscore(i - 1, 0);
		zj(1) = zscore(i - 1, 1);

		arma::mat var_zj_inverse = inv(rj*Al_Ga + sigma);
		arma::mat E_uj_zj = rj * sigma_gamma*D_gamma.t()*var_zj_inverse*zj;
		arma::mat V_uj_zj = rj * sigma_gamma - pow(rj*sigma_gamma, 2)*D_gamma.t()*var_zj_inverse*D_gamma;
		double E_uj_uj = pow(E_uj_zj(0, 0), 2) + V_uj_zj(0, 0);

		arma::mat E_bj_zj = rj * sigma_beta*D_alpha.t()*var_zj_inverse*zj;
		arma::mat V_bj_zj = rj * sigma_beta - pow(rj*sigma_beta, 2)*D_alpha.t()*var_zj_inverse*D_alpha;
		double E_bj_bj = pow(E_bj_zj(0, 0), 2) + V_bj_zj(0, 0);

		arma::mat alpha_num = E_bj_zj * D.t()*sigma_inv*zj / lj - E_bj_zj * D.t()*sigma_inv*D_gamma*E_uj_zj / lj - sqrt_n * sigma_inv(0, 1)*E_bj_bj / lj;

		double alpha_den = n1 * sigma_inv(0, 0)*E_bj_bj / lj;


		arma::mat mar_like = -zj.t()*var_zj_inverse*zj / lj + log(det(var_zj_inverse)) / lj;

		par_vec(0) += alpha_num(0, 0);
		par_vec(1) += alpha_den;
		par_vec(2) += E_bj_bj / rj / lj;
		par_vec(3) += E_uj_uj / rj / lj;
		par_vec(4) += mar_like(0, 0) / 2;
		par_vec(5) += 1 / lj;
	}
	double alpha_update = par_vec(0) / par_vec(1);
	double sigma_beta_update = par_vec(2) / par_vec(5);
	double sigma_gamma_update = par_vec(3) / par_vec(5);
	double logL = par_vec(4);

	arma::rowvec out(4);
	out(0) = alpha_update;
	out(1) = sigma_beta_update;
	out(2) = sigma_gamma_update;
	out(3) = logL;
	return(out);
}

//[[Rcpp::export()]]
arma::colvec NR(int n1,int n2,int num_snp,int num_per,arma::mat sigma,double alpha,double sigma_beta,
	double sigma_gamma, arma::colvec ldscore, arma::mat zscore) {

	double sqrt_n = sqrt(n1)*sqrt(n2);

	arma::colvec initial_est(3);
	initial_est(0) = alpha;
	initial_est(1) = sigma_beta;
	initial_est(2) = sigma_gamma;

	arma::mat Al_Ga(2, 2);
	Al_Ga(0, 0) = n1 * alpha*alpha*sigma_beta + n1 * sigma_gamma;
	Al_Ga(0, 1) = sqrt_n * alpha*sigma_beta;
	Al_Ga(1, 0) = sqrt_n * alpha*sigma_beta;
	Al_Ga(1, 1) = n2 * sigma_beta;

	arma::rowvec par_out(10);
	par_out.zeros();

	for (int i = 1; i <= num_snp; i++) {

		double lj = ldscore(i - 1);
		double rj = lj + num_snp / num_per;
		arma::colvec zj(2);
		zj(0) = zscore(i - 1, 0);
		zj(1) = zscore(i - 1, 1);

		arma::mat H_inv = inv(sigma + rj * Al_Ga);

		/************************

		First Order Derivatives


		************************/

		arma::mat sigma_beta_mat(2, 2);
		sigma_beta_mat(0, 0) = n1 * alpha*alpha;
		sigma_beta_mat(0, 1) = sqrt_n * alpha;
		sigma_beta_mat(1, 0) = sqrt_n * alpha;
		sigma_beta_mat(1, 1) = n2;
		arma::mat H_sigma_beta = rj * sigma_beta_mat; //Sigma_beta
		arma::mat H_inv_sigma_beta = H_inv * H_sigma_beta;
		arma::mat H_inv_sigma_beta_inv = zj.t()*H_inv_sigma_beta*H_inv*zj;

		double sigma_beta_der = -0.5*trace(H_inv_sigma_beta) / lj + 0.5*H_inv_sigma_beta_inv(0, 0) / lj;

		arma::mat alpha_mat(2, 2); //Alpha
		alpha_mat(0, 0) = 2 * n1*alpha;
		alpha_mat(0, 1) = sqrt_n;
		alpha_mat(1, 0) = sqrt_n;
		alpha_mat(1, 1) = 0;
		arma::mat H_alpha = rj * sigma_beta*alpha_mat;
		arma::mat H_inv_alpha = H_inv * H_alpha;
		arma::mat H_inv_alpha_inv = zj.t()*H_inv_alpha*H_inv*zj;

		double alpha_der = -0.5*trace(H_inv_alpha) / lj + 0.5*H_inv_alpha_inv(0, 0) / lj;

		arma::mat sigma_gamma_mat(2, 2);
		sigma_gamma_mat.zeros();
		sigma_gamma_mat(0, 0) = n1;
		arma::mat H_sigma_gamma = rj * sigma_gamma_mat;
		arma::mat H_inv_sigma_gamma = H_inv * H_sigma_gamma;
		arma::mat H_inv_sigma_gamma_inv = zj.t()*H_inv_sigma_gamma*H_inv*zj;

		double sigma_gamma_der = -0.5*trace(H_inv_sigma_gamma) / lj + 0.5*H_inv_sigma_gamma_inv(0, 0) / lj;


		/***************************

		Second Order Derivatives


		*****************************/
		arma::mat alpha_alpha_mat(2, 2);
		alpha_alpha_mat.zeros();
		alpha_alpha_mat(0, 0) = 2 * n1;
		arma::mat H_alpha_alpha = rj * sigma_beta*alpha_alpha_mat;

		arma::mat H_alpha_sigma_beta = rj * alpha_mat;

		arma::cube H_inverse_paramter(2, 2, 3, arma::fill::zeros);
		H_inverse_paramter.slice(0) = H_inv_alpha;
		H_inverse_paramter.slice(1) = H_inv_sigma_beta;
		H_inverse_paramter.slice(2) = H_inv_sigma_gamma;

		par_out(0) += alpha_der;
		par_out(1) += sigma_beta_der;
		par_out(2) += sigma_gamma_der;

		arma::rowvec par_j(6);
		for (int i = 0; i<3; i++) {
			for (int j = 0; j <= i; j++) {
				arma::mat H_inv_par_inv_par = H_inverse_paramter.slice(i)*H_inverse_paramter.slice(j);
				arma::mat H_inv_par_inv_par_inv = zj.t()*H_inv_par_inv_par*H_inv*zj;
				int index = i * (i + 1) / 2 + j;
				par_j(index) = 0.5*trace(H_inv_par_inv_par) / lj - H_inv_par_inv_par_inv(0, 0) / lj;
			}
		}

		//For alpha_alpha
		arma::mat H_inv_alpha_inv_alpha = H_inv_alpha * H_inv_alpha;
		arma::mat H_inv_alpha_inv_alpha_inv = zj.t()*H_inv_alpha_inv_alpha*H_inv*zj;
		arma::mat H_inv_alpha_alpha = H_inv * H_alpha_alpha;
		arma::mat H_inv_alpha_alpha_inv = zj.t() * H_inv_alpha_alpha * H_inv *zj;
		par_j(0) = 0.5*trace(H_inv_alpha_inv_alpha - H_inv_alpha_alpha) / lj - H_inv_alpha_inv_alpha_inv(0, 0) / lj + 0.5*H_inv_alpha_alpha_inv(0, 0) / lj;

		//For alpha_sigma_beta
		arma::mat H_inv_alpha_inv_sigma = H_inv_alpha * H_inv_sigma_beta;
		arma::mat H_inv_alpha_inv_sigma_inv = zj.t()*H_inv_alpha_inv_sigma*H_inv*zj;
		arma::mat H_inv_alpha_sigma_beta = H_inv * H_alpha_sigma_beta;
		arma::mat H_inv_alpha_sigma_inv = zj.t() * H_inv_alpha_sigma_beta * H_inv *zj;
		par_j(1) = 0.5*trace(H_inv_alpha_inv_sigma - H_inv_alpha_sigma_beta) / lj - H_inv_alpha_inv_sigma_inv(0, 0) / lj + 0.5*H_inv_alpha_sigma_inv(0, 0) / lj;

		for (int i = 0; i<6; ++i) {
			par_out(i + 3) += par_j(i);
		}
		arma::mat LogL = zj.t()*H_inv*zj;
		par_out(9) += -0.5*log(det(sigma + rj * Al_Ga)) / lj - 0.5*LogL(0, 0) / lj;
	}
	arma::colvec sum_first_der(3);
	for (int i = 0; i<3; i++) {
		sum_first_der(i) = par_out(i);
	}
	arma::mat sum_sec_der(3, 3);
	for (int i = 0; i<3; i++) {
		for (int j = 0; j<3; j++) {
			sum_sec_der(i, j) = par_out(i*(i + 1) / 2 + j + 3);
		}
	}
	arma::mat tot_sum_sec_der = symmatl(sum_sec_der);
	arma::colvec update = initial_est - inv(tot_sum_sec_der)*sum_first_der;
	//Rcout<<update;
	arma::colvec res(4);
	for (int i = 0; i<3; i++) {
		res(i) = update(i);
	}
	res(3) = par_out(9);
	return(res);

}