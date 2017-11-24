// firfilter.cpp : Defines the entry point for the console application.
// Reference : Digital Signal Processing by Proakis

#include "stdafx.h"
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <conio.h>

using namespace std;

//Matrix multiplication

void multiply(double **x, double **y, double **z)
{
	int N;
	N = sizeof(x) / sizeof(*x);
	int i, j, k;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			z[i][j] = 0;
			for (k = 0; k < N; k++)
				z[i][j] += x[i][k] * y[k][j];
		}
	}

}

// Transpose 
void *transpose(double **x, double **y)
{
	int row, col;
	y = new double *[col]; //creates a new array of pointers to int objects
	for (int i = 0; i < col; i++)
		y[i] = new double[row];

	// transposing
	for (int i = 0; i<row; i++)
		for (int j = 0; j<col; j++)
			y[j][i] = x[i][j];
	return y;
}


void mod(int n, int N, int m) //Computes m = (n mod N) index; m = mod(n, N)
{
	m = remainder(n, N);
	m = m + N;
	m = remainder(m, N);

}

void *circshift(double **x, int m, int N)
{
	//z = output sequence containing the circular shift
	//x = input sequence of length <= N
	//m = sample shift, N = size of circular buffer

	double **a, **y; double **z;

	int len, len_x;
	len_x = sizeof(x) / sizeof(*x);
	len = N - len_x;

	//Check for length of x
	if (len_x > N)
	{
		std::cout << "error('N must be >= the length of x')" << std::endl;
	}
	for (int j = 0; j <len; j++)
	{
		a[1][j] = { 0 };
	}
	for (int i = 0; i <len; i++)
	{
		y[i] = y[i];
		y[i + len] = a[i];
	}
	for (int n = 0; n <= (N - 1); n++)
	{

		mod(n - m, N, n);
		*z = y[n + 1];
	}
	return *z;
}

//convolution with succesive blocks, N-point circular convolution between x1 and x2.

void *circonvt(double **x1, double **x2, int N)
{

	//y = output sequence containing the circular convolution
	// x1 = input sequence of length N1 <= N, x2 = input sequence of length N2 <= N
	//N = size of circular buffer
	//Method: y(n) = sum(x1(m)*x2((n - m) mod N))

	double **a, **b, **x_1, **x_2, **y;
	int len_1, len_2;
	len_1 = N - (sizeof(x1) / sizeof(*x1));
	len_2 = N - (sizeof(x2) / sizeof(*x2));

	//Check for length of x1
	if (len_1 > N)
	{
		std::cout << "error (N must be >= the length of x1)" << std::endl;
	}

	// Check for length of x2
	if (len_2 > N)
	{
		std::cout << "error(N must be <=the length of x2)" << std::endl;
	}

	for (int j = 0; j <len_1; j++)
	{
		a[1][j] = { 0 };
	}

	for (int j = 0; j <len_2; j++)
	{
		b[1][j] = { 0 };
	}

	for (int i = 0; i <len_1; i++)
	{
		x_1[i] = x1[i];
		x_1[i + len_1] = a[i];
	}
	for (int i = 0; i <len_2; i++)
	{
		x_2[i] = x2[i];
		x_2[i + len_2] = b[i];
	}
	for (int m = 0; m <= (N - 1); m++)
	{
		int u;
		mod(-m, N, u);
		*x_2 = x_2[u + 1];
	}

	double **H;
	for (int i = 0; i <= N; i++)
	{
		H[i][i] = { 0 };
	}

	int H_row, H_col;
	H_row = (sizeof(H) / sizeof(H[0]));
	H_col = (sizeof(H) / sizeof(H[0][0])) / H_row;

	for (int n = 0; n <= N; n++)
	{
		for (int i = 0; i < H_col; i++)
		{
			double H = *(double *)circshift;
			H[n][i] = circshift(x_2, n - 1, N);
		}
	}

	transpose(H, H);

	multiply(x_1, H, y);

	return *y;
}


int main()
{
	double delta_f, w_p, w_s, w_c, w_tr;
	int  N_length;
	const double f_pass = 2000;
	const double f_stop = 10000;
	const double Fs = 48000;
	const double dB = 80;
	const double eps = 2.2204e-16;
	const double pi = 3.1415926535897931;

	//Transition Region
	delta_f = f_stop - f_pass;
	w_p = 3.14*(f_pass / (Fs / 2));  //Passband
	w_s = 3.14*(f_stop / (Fs / 2)); //Stopband
									//cutoff frequency in radians
	w_c = 0.5*(w_p + w_s);
	w_tr = w_s - w_p; //Transition

					  // length of the filter
	N_length = round(dB * Fs / (22.0 * (f_stop - f_pass)));

	//blackman window for stopband attenuation 80dB

	int M1, alpha; double **w_black, **m, **hd, **h;
	M1 = N_length - 1;

	for (int n = 0; n <= M1; n++)
	{
		*w_black[n] = abs(0.42 - 0.5*cos(2 * 3.14* n / (M1)) + 0.08*cos(4 * 3.14*n / (M1)));
		//std::cout << w_black[n] << std::endl;

		// Ideal low pass filter

		alpha = M1 / 2;
		*m[n] = n - alpha + eps;
		//std::cout << m[n] << std::endl;
		*hd[n] = sin(w_c* *m[n]) / (pi* *m[n]); //ideal impulse response between 0 to M1
												//std::cout << hd[n] << std::endl;

												//Impulse response of the filter
		*h[n] = *hd[n] * *w_black[n];
		std::cout << (sizeof(h) / sizeof(*h)) << std::endl;
	}

	// Read Input audio wav file (data.txt)

	double **x;
	std::ifstream in("data.txt");
	if (!in) {
		std::cout << "Cannot open file.\n";
	}
	for (int j = 0; j < 1; j++) {
		for (int i = 0; i < 48000; i++) {
			in >> x[i][j];
		}
	}
	in.close();

	// overlap save method

	transpose(x, x);
	int len_x, len_h; double **a, **h_1, **z_1, **z_2, **x_1;
	len_x = sizeof(x) / sizeof(*x);
	len_h = sizeof(h) / sizeof(*h);
	int M_O, L_O, Z_O;
	M_O = len_h - 1;
	L_O = N_length - M_O;
	Z_O = N_length - L_O;

	for (int j = 0; j <= Z_O; j++)
	{
		a[j] = { 0 };
	}
	for (int i = 0; i <= Z_O; i++)
	{
		h_1[i] = h[i];
		h_1[i + Z_O] = a[i];
	}

	// Preappend (M-1)zeros

	for (int j = 0; j <= L_O; j++)
	{
		z_1[j] = { 0 };
	}
	for (int j = 0; j <= (N_length - 1); j++)
	{
		z_2[j] = { 0 };
	}
	for (int i = 0; i <= L_O; i++)
	{
		x_1[i] = h[i];
		x_1[i + L_O] = z_1[i];
		x_1[i + L_O + L_O] = z_2[i];
	}

	double K; double  **Y, **x_k, **y;
	K = floor((len_x + L_O - 1) / (Z_O));  //no.of blocks, round
	for (int i = 0; i <= (K + 1); i++)
		for (int j = 0; j <= N_length; j++)
		{
			{
				Y[i][j] = { 0 };
			}
		}

	int Y_row, Y_col;
	Y_row = (sizeof(Y) / sizeof(Y[0]));
	Y_col = (sizeof(Y) / sizeof(Y[0][0])) / Y_row;

	for (int k = 0; k <= K; k++)
	{
		for (int k = (k*Z_O + 1); k <= (k*Z_O + N_length); k++)
		{
			*x_k = x_1[k];
		}
		for (int i = 0; i < Y_col; i++)
		{
			//circular convolution
			double Y = *(double *)circonvt;
			Y[k + 1][i] = circonvt(x_k, h_1, N_length);
		}
	}

	for (int i = M_O; i <= N_length; i++)
	{
		for (int j = 0; j<Y_row; j++)
			**Y = Y[j][i];
	}
	transpose(Y, Y); // discard the first(M - 1) samples

	for (int q = 0; q < Y_row; q++)
	{
		for (int t = 0; t < Y_col; t++)
		{
			*Y[q * Y_col + t] = Y[q][t];
		}
	}
	transpose(Y, y);     // Output sequence
	std::cout << y << std::endl;

	const int size = 48021;
	ofstream out("output.txt");
	if (out.is_open())
	{
		for (int count = 0; count < size; count++) {
			out << y[count] << " ";
		}
		out.close();
	}
	else cout << "Unable to open file";

	return 0;
}


