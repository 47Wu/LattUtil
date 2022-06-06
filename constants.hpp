#pragma once
#include <complex>
#include <string>
#include <iostream>
#include <armadillo>


#ifndef MAXNATOM
#define MAXNATOM 1024
#endif // !MAXNATOM

const std::string LUversion = "0.0";

typedef std::complex<double> dcomplex;
#define MKL_Complex16 dcomplex

enum errcode {
	ec_norm=0,
	ec_io,
	ec_math,
	ec_mpi,
	ec_alloc
};

enum latttype {
	frac, 
	cart
};

const dcomplex cmplx_0(0.0, 0.0);
const dcomplex cmplx_1(1.0, 0.0);
const dcomplex cmplx_i(0.0, 1.0);

const double eps2 = 1E-2;
const double eps3 = 1E-3;
const double eps4 = 1E-4;
const double eps5 = 1E-5;
const double eps6 = 1E-6;
const double eps7 = 1E-7;
const double eps8 = 1E-8;

const double infinit_small = eps8;

const double thr_2 = eps2;
const double thr_3 = eps3;
const double thr_4 = eps4;
const double thr_5 = eps5;
const double thr_6 = eps6;
const double thr_7 = eps7;
const double thr_8 = eps8;

// thresholds for some processes
const double prec_default = infinit_small;
const double prec_lattice = thr_5;
const double thr_angle = thr_2;
const double thr_lattice = thr_5;

const double pi = 3.14159265358979323846;
const double twopi = 2.0 * pi;

const double sqrt2 = 1.414213562373095048801688;
const double sqrt2_2 = sqrt2 / 2.0;
const double sqrt3 = 1.732050807568877293527446;
const double sqrt5 = 2.236067977499789696409173;
const double sqrt6 = sqrt2 * sqrt3;
const double sqrt7 = 2.6457513110645907;
const double sqrt11 = 3.3166247903554;

const double pi2_3 = pi * 2.0 / 3.0;
const double pi4_3 = pi * 4.0 / 3.0;
const double pi_2 = pi / 2.0;
const double pi3_2 = pi * 3.0 / 2.0;
const double pi_3 = pi / 3.0;
const double pi5_3 = pi * 5.0 / 3.0;

const double pi_4 = pi / 4.0;
const double pi3_4 = pi * 3.0 / 4.0;
const double pi_6 = pi / 6.0;
const double pi5_6 = pi * 5.0 / 6.0;

namespace paulimat {
	const arma::cx_mat22 sigma_0 = { {cmplx_1, cmplx_0}, {cmplx_0, cmplx_1} };
	const arma::cx_mat22 sigma_x = { {cmplx_0, cmplx_1}, {cmplx_1, cmplx_0} };
	const arma::cx_mat22 sigma_y = { {cmplx_0, -cmplx_i}, {cmplx_i, cmplx_0} };
	const arma::cx_mat22 sigma_z = { {cmplx_1, cmplx_0}, {cmplx_0, -cmplx_1} };
	const arma::cx_mat22 isigma_y = { {cmplx_0, cmplx_1}, {-cmplx_1, cmplx_0} };
	const arma::cx_mat22 sigma_d0 = isigma_y;
	const arma::cx_mat22 sigma_dx = sigma_x * isigma_y;
	const arma::cx_mat22 sigma_dy = sigma_y * isigma_y;
	const arma::cx_mat22 sigma_dz = sigma_z * isigma_y;
	const arma::sp_cx_mat sp_sigma_d0 = arma::conv_to<arma::sp_cx_mat>::from(sigma_d0);
	const arma::sp_cx_mat sp_sigma_dx = arma::conv_to<arma::sp_cx_mat>::from(sigma_dx);
	const arma::sp_cx_mat sp_sigma_dy = arma::conv_to<arma::sp_cx_mat>::from(sigma_dy);
	const arma::sp_cx_mat sp_sigma_dz = arma::conv_to<arma::sp_cx_mat>::from(sigma_dz);
}

#ifdef __DEBUG
// employ some debuggers
static int ccxu_count = 0;
static int gxzhi_count = 0;
#endif