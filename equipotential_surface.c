/**
 *
 * Copyright (c) 2024, Jason Koci
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 * Author: Jason Koci <iasonaskotsis@hotmail.com>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define SCOEFF 2401333 /* Number of the the fully normalized dimensionless spherical harmonic coefficients (Cnm, Snm) */
#define GM 0.3986004415e15 /* Geocentric gravitational constant (s^-2) */
#define G 6.673e-11 /* Gravitational constant (kg^-1 s^-2) */
#define R 0.63781363e7 /* Reference radius of the GGM (m) */
#define OMEGA 7.29211500000e-5 /* Angular velocity (1 / s) */
#define DENSITY 2670.0 /* Density (kg / m^3) */
#define CONVTOL 1e-5 /* Convergence tolerance */
#define MYABS(x) ((x) > 0) ? (x):-(x)

const double pi = 4.0l * atan(1.0l);
const double pih = pi / 2.0l;
const double d2r = pi / 180.0l; /* Conversion from deg to rad */

/* Coefficients for the computation of the ALFs */
double *root, *rooti, *norm1, *norm2;
double *PNN, *PN0, *Tni_odd, *Tni_even, sqr2, sqr3;
double *anm, *bnm, *cnm;

void geomalf(int, double, double *);
void cholesky_tridiagonal(double *, double *, double *, int);

int main(void)
{
	register int i, j, t, q, n, m;
	char ch, line[BUFSIZ], str[BUFSIZ];
	int r, p, N, STnieven, STniodd, S, header1, header2, temp;
	int n2, r2, np1, nm1, nm3, r2m1, i_1, points, iterations;
	double coeff, coefn, z, A, B, C, pnn, pn0, sign;
	double *pnm, *Snm, *Cnm, _snm, _cnm, test;
	double cosnm1x, cosnm2x, cosnx, sinnm1x, sinnm2x, sinnx, *cosn, *sinn;
	double Wa, DiffWa, deltaW, DiffdeltaW, lambda, theta, rho;
	double *colatitude, *longitude, *radius, *wpi;
	double *bij, *Di, *si, *mij, *wi, *di;
	double *Wp, *DiffW;
	double sigma01, sigma0;
	double bii, bi_ip1, Fi0, ridi, summation1, resi, resii;
	double frac, var, sintheta, sin2theta;
	FILE *fp1, *fp2;

	N = 2190; /* The maximum degree of the ALFs */
	char file1[BUFSIZ] = "EGM2008.gfc"; /* The name of the global gravity model file */
	header1 = 22; /* The header of the global gravity model file */
	char file2[BUFSIZ] = "spherical_coords_geoid_step_10.txt"; /* The name of the measurements file */
	header2 = 1; /* The header of the measurements file */

	STnieven = (N % 2) ? ((N - 1) * (N - 3) / 8):(N * (N - 2) / 8); /* The number of the Tni ratios for n even */
	STniodd = (N % 2) ? ((N * N - 1) / 8):(N * (N - 2) / 8); /* The number of the Tni ratios for n odd */
	S = (N + 1) * (N + 2) / 2; /* The number of the ALFs */

	if((fp1 = fopen(file1, "r")) == NULL)
	{
		printf("\nCant open the global gravity model file");
		exit(1000);
	}

	if((fp2 = fopen(file2, "r")) == NULL)
	{
		printf("\nCant open the measurements file (spherical coordinates)");
		exit(2000);
	}

	/* Counting the measurements */
	i = 0;
	while ((ch = getc(fp2)) != EOF)
		if (ch == '\n')
			i++;
	rewind(fp2);

	points = i - header2; /* The number of points */
	r = points - 1; /* Degrees of freedom */
	r2 = 2 * r;
	r2m1 = r2 - 1;
	printf("\nN = %d\nS = %d ALFs\nPoints = %d\nDegrees of freedom = %d", N, S, points, r);

	/* Allocating memory */
	if ((root = malloc((2 * N + 5) * sizeof(double))) == NULL)
		exit(1);
	else if ((rooti = malloc((2 * N + 4) * sizeof(double))) == NULL)
		exit(2);
	else if ((norm1 = malloc((N + 1) * sizeof(double))) == NULL)
		exit(3);
	else if ((norm2 = malloc(N * sizeof(double))) == NULL)
		exit(4);
	else if ((PNN = malloc((N + 1) * sizeof(double))) == NULL)
		exit(5);
	else if ((PN0 = malloc((N + 1) * sizeof(double))) == NULL)
		exit(6);
	else if ((Tni_even = malloc(STnieven * sizeof(double))) == NULL)
		exit(7);
	else if ((Tni_odd = malloc(STniodd * sizeof(double))) == NULL)
		exit(8);
	else if ((pnm = malloc(S * sizeof(double))) == NULL)
		exit(9);
	else if ((Snm = malloc(SCOEFF * sizeof(double))) == NULL)
		exit(10);
	else if ((Cnm = malloc(SCOEFF * sizeof(double))) == NULL)
		exit(11);
	else if ((cosn = malloc((N + 1) * sizeof(double))) == NULL)
		exit(12);
	else if ((sinn = malloc((N + 1) * sizeof(double))) == NULL)
		exit(13);
	else if ((colatitude = malloc(points * sizeof(double))) == NULL)
		exit(14);
	else if ((longitude = malloc(points * sizeof(double))) == NULL)
		exit(15);
	else if ((radius = malloc(points * sizeof(double))) == NULL)
		exit(16);
	else if ((wpi = malloc(points * sizeof(double))) == NULL)
		exit(17);
	else if ((anm = malloc(S * sizeof(double))) == NULL)
		exit(18);
	else if ((bnm = malloc(S * sizeof(double))) == NULL)
		exit(19);
	else if ((cnm = malloc(S * sizeof(double))) == NULL)
		exit(20);
	else if ((Di = malloc(r * sizeof(double))) == NULL)
		exit(21);
	else if ((si = malloc(r * sizeof(double))) == NULL)
		exit(22);
	else if ((di = malloc(points * sizeof(double))) == NULL)
		exit(23);
	else if ((bij = malloc((2 * r) * sizeof(double))) == NULL)
		exit(24);
	else if ((wi = malloc(r * sizeof(double))) == NULL)
		exit(25);
	else if ((Wp = malloc(points * sizeof(double))) == NULL)
		exit(26);
	else if ((DiffW = malloc(points * sizeof(double))) == NULL)
		exit(27);
	else if ((mij = malloc(((r * (r + 1)) / 2) * sizeof(double))) == NULL)
		exit(28);

	/*===============================GGM MODEL=================================*/
	for (i = 0; i < header1; i++)
	{
		if(fgets(line, BUFSIZ, fp1) == NULL)
				printf("\nerror in function fgets\n");
	}

	Cnm[0] = 1.0l;
	Snm[0] = 0.0l;
	Cnm[1] = 0.0l;
	Snm[1] = 0.0l;
	Cnm[2] = 0.0l;
	Snm[2] = 0.0l;
	i = 3;
	while(fgets(line, BUFSIZ, fp1))
	{
		sscanf(line, "%s%d%d%lf%lf%lf%lf", str, &temp, &temp, &_cnm, &_snm, &test, &test);
		Cnm[i] = _cnm;
		Snm[i] = _snm;
		i++;
	}
	/*=========================================================================*/

	/*============================Measurements=================================*/
	for (i = 0; i < header2; i++)
	{
		if(fgets(line, BUFSIZ, fp2) == NULL)
				printf("\nerror in function fgets\n");
	}

	i = 0;
	while(fgets(line, BUFSIZ, fp2))
	{
		sscanf(line, "%lf%lf%lf", &lambda, &theta, &rho);
		colatitude[i] = theta;
		longitude[i] = lambda;
		radius[i] = rho;
		i++;
	}
	/*=========================================================================*/

	/*============================GLOBAL VARIABLES ============================*/
	/* Square roots look up table */
	for (t = 0; t <= 2 * N + 4; t++)
		root[t] = sqrt(1.0l * t);
	p = sqrt(2.0l * N + 4) + 1;
	for (t = 0; t < p; t++)
	{
		q = t * t;
		root[q] = 1.0l * t;
	}

	sqr2 = root[2];
	sqr3 = root[3];
	coeff = 0.0l;
	for (t = 0; t <= 2 * N + 3; t++)
	{
		z = root[t + 1];
		rooti[t] = coeff * z;
		coeff = z;
	}

	/* Computation of the normalization coefficients */
	norm1[0] = 1.0l;
	for (n = 1; n <= N; n++)
	{
		coeff = root[2 * n + 1];
		norm1[n] = coeff;
		norm2[n] = sqr2 * coeff / root[n] / root[n + 1];
	}

	/* Computation of the odd pnn coefficients */
	pnn = 1.0l;
	PNN[1] = pnn;
	for (n = 3; n <= N; n += 2)
	{
		p = n - 1;
		A = 1.0l / p;
		B = 1.0l / n;
		coeff = 1.0l - 0.75L * B;
		pnn *= 1.0l - A * coeff;
		PNN[n] = pnn;
	}

	/* Computation of the even pnn coefficients */
	pnn = 2.0l;
	for (n = 2; n <= N; n += 2)
	{
		p = n - 1;
		A = 1.0l / p;
		B = 1.0l / n;
		coeff = 1.0l - 0.75L * B;
		pnn *= 1.0l - A * coeff;
		PNN[n] = pnn;
	}

	/* Computation of the pn0 coefficients */
	pn0 = 1.0l;
	PN0[0] = pn0;
	for (n = 2; n <= N; n += 2)
	{
		B = 1.0l / n;
		A = 1.0l - B;
		pn0 *= A * A;
		PN0[n] = pn0;
	}

	/* Computation of the odd Tni ratios */
	i = 0;
	for (n = 3; n <= N; n += 2)
	{
		p = n - 1;
		q = n + 2;
		while (p > 0)
		{
			/* Calculates the Tni ratio */
			B = 1.0l - 1.0l / p;
			C = 1.0l + 1.0l / q;
			Tni_odd[i] = B * C;
			p -= 2;
			q += 2;
			i++;
		}
	}

	/* Computation of the even Tni ratios */
	i = 0;
	for (n = 2; n <= N; n += 2)
	{
		p = n - 2;
		q = n + 3;
		while (p > 0)
		{
			/* Calculates the Tni ratio */
			B = 1.0l - 1.0l / p;
			C = 1.0l + 1.0l / q;
			Tni_even[i] = B * C;
			p -= 2;
			q += 2;
			i++;
		}
	}

	/* Computation of the anm, bnm and cnm coefficients */
	for (n = 4; n <= N; n++)
	{
		n2 = 2 * n;
		np1 = n + 1;
		nm1 = n - 1;
		nm3 = n - 3;
		coefn = root[n2 + 1] / root[n2 - 3];
		for (m = 2; m <= n; m++)
		{
			sign = (m > 2) ? 1.0l:sqr2;
			q = m + (n * np1) / 2;
			//printf("\nq = %d", q);
			z = rooti[nm1 + m];
			anm[q] = coefn * rooti[nm1 - m] / z;
			bnm[q] = sign * coefn * rooti[nm3 + m] / z;
			cnm[q] = sign * rooti[np1 - m] / z;
		}
	}
	/*=========================================================================*/

	/*===================Least-squares adjustment procedure====================*/

	/* Setting initial values for the radial distances and weights for each point */
	for (i = 0; i < points; i++)
	{
		di[i] = 0.0l;
		wpi[i] = 1.0l; /* Weight */
	}

	sigma0 = 1.0l;
	sigma01 = 2.0l;
	iterations = 0;
	do {
		sigma0 = sigma01;
		/* Computation of the gravity potential and its derivatives for each point */
		printf("\n\nComputation of the gravity potential and its derivatives for each point: Start (%lds)", clock() / CLOCKS_PER_SEC);
		for (i = 0; i < points; i++)
		{
			theta = colatitude[i];
			lambda = longitude[i] * d2r;
			ridi = radius[i] - di[i];

			/* Computation of the multiple cosines by using the Chebyshev's method */
			cosnm1x = cos(lambda);
			coeff = 2.0l * cosnm1x;
			cosnm2x = 1.0l;
			cosn[0] = 1.0l;
			cosn[1] = cosnm1x;

			for (n = 2; n <= N; n++)
			{
				cosnx = coeff * cosnm1x - cosnm2x;
				cosn[n] = cosnx;
				cosnm2x = cosnm1x;
				cosnm1x = cosnx;
			}

			/* Computation of the multiple sines by using the Chebyshev's method */
			coeff = 2.0l * cosn[1];
			sinnm2x = 0.0l;
			sinnm1x = sin(lambda);
			sinn[0] = 0.0l;
			sinn[1] = sinnm1x;

			for (n = 2; n <= N; n++)
			{
				sinnx = coeff * sinnm1x - sinnm2x;
				sinn[n] = sinnx;
				sinnm2x = sinnm1x;
				sinnm1x = sinnx;
			}
			/*==================================================================*/

			Wa = DiffWa = 0.0l;
			geomalf(N, theta, pnm); /* Calling the function geomalf(), to compute the ALFs */
			frac = R / ridi;
			coeff = 1.0l;
			j = 0;
			for(n = 0; n <= N; n++)
			{
				summation1 = 0.0l;
				for(m = 0; m <= n; m++)
				{
					summation1 += pnm[j] * (Cnm[j] * cosn[m] + Snm[j] * sinn[m]);
					j++;
				}
				var = coeff * summation1;
				Wa += var;
				DiffWa += (n + 1) * var;
				coeff *= frac;
			}

			sintheta = sin(theta * d2r);
			sin2theta = sintheta * sintheta;

			Wa *= GM / ridi;
			Wa += 0.5l * ridi * ridi * OMEGA * OMEGA * sin2theta;
			DiffWa *= GM / ridi / ridi;
			DiffWa -= ridi * OMEGA * OMEGA * sin2theta;

			/* Indirect effect of the topographic masses */
			if(di[i] <= 0)
			{
				Wp[i] = Wa;
				DiffW[i] = DiffWa;
			}
			else
			{
				var = 2.0l * pi * G * DENSITY * di[i];
				deltaW = var * di[i] * (1.0l + 2.0l * radius[i] / ridi) / 3.0l;
				Wp[i] = Wa - deltaW;
				DiffdeltaW = 2.0l * var * (1.0l + (radius[i] / ridi) * (1.0l + radius[i] / ridi)) / 3.0l;
				DiffW[i] = DiffWa - DiffdeltaW;
			}
		}
		printf("\nComputation of the gravity potential and its derivatives for each point: End (%lds)", clock() / CLOCKS_PER_SEC);

		printf("\nComputation of the elements of the matrix B");

		/* Computation of the elements of the design matrix B */
		for (i = 0; i < r - 1; i++)
		{
			bii = DiffW[i];
			bi_ip1 = -DiffW[i + 1];
			bij[2 * i] = bii;
			bij[2 * i + 1] = bi_ip1;
			Di[i] = bii * bii / wpi[i] + bi_ip1 * bi_ip1 / wpi[i + 1]; /* Computation of the diagonal elements of the matrix M */
		}
		bii = DiffW[r - 1];
		bi_ip1 = -DiffW[r];
		bij[2 * r - 2] = bii;
		bij[2 * r - 1] = bi_ip1;
		Di[i] = bii * bii / wpi[r - 1] + bi_ip1 * bi_ip1 / wpi[r];

		/* Computation of the elements of the matrix M */
		for (i = 0; i < r - 1; i++)
			si[i] = bij[2 * i + 1] * bij[2 * i + 2] / wpi[i + 1];

		printf("\nInversion of the matrix M");

		/* Inversion of the matrix M */
		cholesky_tridiagonal(Di, si, mij, r);

		printf("\nComputation of the vector w");

		/* Computation of the vector w */
		for (i = 0; i < r; i++)
		{
			Fi0 = Wp[i] - Wp[i + 1];
			wi[i] = bij[2 * i] * di[i] + bij[2 * i + 1] * di[i + 1] - Fi0;
		}

		printf("\nComputation of the radial distances di");

		/*========== Computation of the radial distances di ==========*/

		/* Computation of the radial distance d1 */
		resii = coeff = 0.0l;
		for (j = 0; j < r; j++)
				resii += mij[j] * wi[j];
		di[0] = bij[0] * resii / wpi[0];

		/* Computation of the radial distances di for 1 < i < n */
		for (i = 1; i < r; i++)
		{
			i_1 = i - 1;
			resi = resii = 0.0l;

			for (j = 0; j < i_1; j++)
				resi += mij[j * (r2m1 - j) / 2 + i_1] * wi[j];
			for (j = i_1; j < r; j++)
				resi += mij[i_1 * (r2 - i) / 2 + j] * wi[j];
			resi *= bij[2 * i - 1];

			for (j = 0; j < i; j++)
				resii += mij[j * (r2m1 - j) / 2 + i] * wi[j];
			for (j = i; j < r; j++)
				resii += mij[i * (r2m1 - i) / 2 + j] * wi[j];
			resii *= bij[2 * i];
			resi += resii;
			di[i] = resi / wpi[i];
		}

		/* Computation of the radial distance dn */
		resii = 0.0l;
		for (j = 0; j < r; j++)
				resii += mij[(j + 1) * (2 * r - j) / 2 - 1] * wi[j];
		di[r] = bij[2 * r - 1] * resii / wpi[r];

		for (i = 0; i < points; i++)
		{
			resi = di[i];
			coeff += resi * resi * wpi[i];
		}

		sigma01 = sqrt(coeff / r);
		printf("\niteration %d completed", iterations);
		printf("\nsigma0 = +/- %-14.8lf", sigma01);

		/* Correcting the values of the radius for each point */
		for (i = 0; i < points; i++)
			radius[i] -= di[i];

		iterations++;
		if (iterations > 5)
			break;

	}while(MYABS(sigma0 - sigma01) > CONVTOL);

	/* Printing the results */
	printf("\nn = %d measurements\nr = %d degrees of freedom\niterations = %d\nsigma0 = +/- %-10.4lf\n", points, r, iterations, sigma01);
	printf("\n%-15s%-15s%-20s%-20s%-20s", "Longitude", "Co-latitude", "Radius", "di", "W0\n\n");
	for (i = 0; i < points; i++)
	{
		ridi = radius[i] - di[i];
		printf("%-15.4lf%-15.4lf%-20.3lf%-20.3lf%-20.3lf\n", longitude[i], colatitude[i], ridi, di[i], Wp[i]);
	}

	return 0;
}

/**
 * \brief           This functions computes the fully normalized associated Legendre functions (ALFs)
 *
 * \param[in]       N: The maximum degree of the ALFs
 * \param[in]       theta: Colatitude
 * \param[in]       pnm: The vector for the ALFs
 *
 */
void geomalf(int N, double theta, double *pnm)
{
	register int p, t, i;
	int n, m, q, np1, nm1, nm2, nm3, nnp12;
	double pihalf, raddeg, cosnm1x, cosnm2x, cosnx, sinnm1x, sinnm2x, sinnx;
	double coeff, coss, sins, Pi0, Pi1, Tni;
	double *cosn, *sinn, *Pn0, *Pn1, *Pn_even, *Pn_odd, *Pn;
	double Pn_2_m, Pn_2_m_2, Pn_m_2, Pnm;

	pihalf = 2.0l * atan(1.0l);
	raddeg = pihalf / 90.0l;
	theta *= raddeg;

	if ((cosn = malloc((N + 1) * sizeof(double))) == NULL)
		exit(100);
	else if ((sinn = malloc((N + 1) * sizeof(double))) == NULL)
		exit(101);
	else if ((Pn0 = malloc((N + 1) * sizeof(double))) == NULL)
		exit(102);
	else if ((Pn1 = malloc((N + 1) * sizeof(double))) == NULL)
		exit(103);
	else if ((Pn_even = malloc((N + 3) * sizeof(double))) == NULL)
		exit(104);
	else if ((Pn_odd = malloc((N + 3) * sizeof(double))) == NULL)
		exit(105);
	else if ((Pn = malloc((N + 3) * sizeof(double))) == NULL)
		exit(106);

	/* Calculates the multiple angle cosines with Chebyshev's algorithm
	https://en.wikipedia.org/wiki/List_of_trigonometric_identities */

	cosnm1x = cos(theta);
	coeff = 2.0l * cosnm1x;
	cosnm2x = 1.0l;
	cosn[0] = 1.0l;
	cosn[1] = cosnm1x;
	for (n = 2; n <= N; n++)
	{
		cosnx = coeff * cosnm1x - cosnm2x;
		cosn[n] = cosnx;
		cosnm2x = cosnm1x;
		cosnm1x = cosnx;
	}

	/* Calculates the multiple angle sines with Chebyshev's algorithm
	https://en.wikipedia.org/wiki/List_of_trigonometric_identities */

	coeff = 2.0l * cosn[1];
	sinnm2x = 0.0l;
	sinnm1x = sin(theta);
	sinn[0] = 0.0l;
	sinn[1] = sinnm1x;
	for (n = 2; n <= N; n++)
	{
		sinnx = coeff * sinnm1x - sinnm2x;
		sinn[n] = sinnx;
		sinnm2x = sinnm1x;
		sinnm1x = sinnx;
	}

	/* Calculates the Legendre Polynomials (m = 0) and Functions (m = 1) */
	Pn0[0] = 1.0l;
	Pn0[1] = root[3] * cosn[1];
	Pn1[0] = 0.0l;
	Pn1[1] = root[3] * sinn[1];

	/* Calculates the odd Legendre polynomials (m = 0) and Functions (m = 1) */
	i = 0;
	coss = cosn[1];
	sins = sinn[1];
	for (n = 3; n <= N; n += 2)
	{
		p = n - 1;
		Pi0 = coss;
		Pi1 = sins;
		t = 3;
		while (p > 0)
		{
			Tni = Tni_odd[i];
			Pi0 *= Tni;
			Pi1 *= Tni;
			Pi0 += cosn[t];
			Pi1 += t * sinn[t];
			t += 2;
			p -= 2;
			i++;
		}
		Pn0[n] = Pi0 * PNN[n] * norm1[n];
		Pn1[n] = Pi1 * PNN[n] * norm2[n];
	}

	/* Calculates the even Legendre polynomials (m = 0) and Functions (m = 1) */
	i = 0;
	coss = cosn[2];
	sins = sinn[2];
	for (n = 2; n <= N; n += 2)
	{
		p = n - 2;
		Pi0 = coss;
		Pi1 = 2.0l * sins;
		t = 4;
		while (p > 0)
		{
			Tni = Tni_even[i];
			Pi0 *= Tni;
			Pi1 *= Tni;
			Pi0 += cosn[t];
			Pi1 += t * sinn[t];
			t += 2;
			p -= 2;
			i++;
		}
		coeff = root[2 * n + 1];
		Pn0[n] = (Pi0 * PNN[n] + PN0[n]) * norm1[n];
		Pn1[n] = Pi1 * PNN[n] * norm2[n];
	}

	n = 0;
	pnm[0] = Pn0[0];

	n = 1;
	Pn_odd[0] = Pn0[n];
	Pn_odd[1] = Pn1[n];
	pnm[1] = Pn_odd[0];
	pnm[2] = Pn_odd[1];

	n = 2;
	Pn_even[0] = Pn0[n];
	Pn_even[1] = Pn1[n];
	Pn_even[2] = sqrt(3.0l) * (sqrt(5.0l) - Pn0[2]) / 3.0l;
	pnm[3] = Pn_even[0];
	pnm[4] = Pn_even[1];
	pnm[5] = Pn_even[2];

	n = 3;
	Pn_odd[0] = Pn0[n];
	Pn_odd[1] = Pn1[n];
	Pn_odd[2] = sqrt(15.0l) * (sqrt(21.0l) * Pn0[1] / 3.0l - Pn0[3]) / 5.0l;
	Pn_odd[3] = sqrt(15.0l) * (sqrt(14.0l) * Pn1[1] - Pn1[3]) / 15.0l;
	pnm[6] = Pn_odd[0];
	pnm[7] = Pn_odd[1];
	pnm[8] = Pn_odd[2];
	pnm[9] = Pn_odd[3];

	i = 10;
	for (n = 4; n <= N; n++)
	{
		np1 = n + 1;
		nm1 = n - 1;
		nm2 = n - 2;
		nm3 = n - 3;
		nnp12 = (n * np1) / 2;

		if(n % 2) /*========== n is odd ==========*/
		{
			/*<========= m is odd =========>*/
			Pn[1] = Pn1[n];

			for (m = 3; m <= nm2; m += 2)
			{
				q = m + nnp12;
				Pn_m_2 = Pn[m - 2];
				Pn_2_m_2 = Pn_odd[m - 2];
				Pn_2_m = Pn_odd[m];
				Pnm = anm[q] * Pn_2_m + bnm[q] * Pn_2_m_2 - cnm[q] * Pn_m_2;
				Pn[m] = Pnm;
			}

			q = m + nnp12;
			Pn_2_m_2 = Pn_odd[nm2];
			Pn_m_2 = Pn[nm2];
			Pnm = bnm[q] * Pn_2_m_2 - cnm[q] * Pn_m_2; /* calculates the Pnn Legendre function */
			Pn[n] = Pnm;

			/*<========= m is even =========>*/

			Pn[0] = Pn0[n]; /* polynomial */
			for (m = 2; m <= nm2; m += 2)
			{
				q = m + nnp12;
				Pn_m_2 = Pn[m - 2];
				Pn_2_m_2 = Pn_odd[m - 2];
				Pn_2_m = Pn_odd[m];
				Pnm = anm[q] * Pn_2_m + bnm[q] * Pn_2_m_2 - cnm[q] * Pn_m_2;
				Pn[m] = Pnm;
			}

			q = m + nnp12;
			Pn_2_m_2 = Pn_odd[nm3];
			Pn_m_2 = Pn[nm3];
			Pnm = bnm[q] * Pn_2_m_2 - cnm[q] * Pn_m_2; /* calculates the Pn,n-1 Legendre function */
			Pn[nm1] = Pnm;

			for (m = 0; m <= n; m++)
			{
				Pnm = Pn[m];
				Pn_odd[m] = Pnm;
				pnm[i++] = Pnm;
			}
		}
		else /*========== n is even ==========*/
		{
			/*<========= m is even =========>*/

			Pn[0] = Pn0[n];
			for (m = 2; m <= nm2; m += 2)
			{
				q = m + nnp12;
				Pn_m_2 = Pn[m - 2];
				Pn_2_m_2 = Pn_even[m - 2];
				Pn_2_m = Pn_even[m];
				Pnm = anm[q] * Pn_2_m + bnm[q] * Pn_2_m_2 - cnm[q] * Pn_m_2;
				Pn[m] = Pnm;
			}

			q = m + nnp12;
			Pn_2_m_2 = Pn_even[nm2];
			Pn_m_2 = Pn[nm2];
			Pnm = bnm[q] * Pn_2_m_2 - cnm[q] * Pn_m_2; /* calculates the Pnn Legendre function */
			Pn[n] = Pnm;

			/*<========= m is odd =========>*/

			Pn[1] = Pn1[n];
			for (m = 3; m <= nm2; m += 2)
			{
				q = m + nnp12;
				Pn_m_2 = Pn[m - 2];
				Pn_2_m_2 = Pn_even[m - 2];
				Pn_2_m = Pn_even[m];
				Pnm = anm[q] * Pn_2_m + bnm[q] * Pn_2_m_2 - cnm[q] * Pn_m_2;
				Pn[m] = Pnm;
			}

			q = m + nnp12;
			Pn_2_m_2 = Pn_even[nm3];
			Pn_m_2 = Pn[nm3];
			Pnm = bnm[q] * Pn_2_m_2 - cnm[q] * Pn_m_2; /* calculates the Pn,n-1 Legendre function */
			Pn[nm1] = Pnm;

			for (m = 0; m <= n; m++)
			{
				Pnm = Pn[m];
				Pn_even[m] = Pnm;
				pnm[i++] = Pnm;
			}
		}
	}
	free(cosn);
	free(sinn);
	free(Pn0);
	free(Pn1);
	free(Pn_odd);
	free(Pn_even);
	free(Pn);
}

/**
 * \brief           This function computes the inverse of a symmetric tridiagonal matrix
 *
 * \param[in]       D: The vector of the diagonal elements of the matrix
 * \param[in]       s: The vector of the other elements of the matrix
 * \param[in]       b: The vector of the elements of the inverse matrix
 * \param[in]       n: The dimension of the matrix
 *
 */
void cholesky_tridiagonal(double *D, double *s, double *b, int n)
{
	register int i, j, k;
	int uts, n2p1, n2m1, nm1, index, indexi, indexj, ind;
	double *c, *d, sumb, coeff;
	double ci_ip1, cii;

	uts = (n * (n + 1)) / 2; /* The number of elements of the upper triangular part */
	n2p1 = 2 * n + 1;
	n2m1 = 2 * n - 1;
	nm1 = n - 1;

	if ((c = malloc(n2m1 * sizeof(double))) == NULL)
		exit(200);
	else if ((d = malloc(uts * sizeof(double))) == NULL)
		exit(201);

	/* calculation of the matrix c (bidiagonal matrix) */
	c[0] = cii = sqrt(D[0]);
	d[0] = 1.0l / cii;
	for (i = 1; i < n; i++)
	{
		ci_ip1 = s[i - 1] / cii;
		cii = sqrt(D[i] - ci_ip1 * ci_ip1);
		d[i * (n2p1 - i) / 2] = 1.0l / cii;
		c[2 * i - 1] = ci_ip1;
		c[2 * i] = cii;
	}

	/* calculation of the matrix d (upper triangular matrix) */
	for (i = n - 1; i >= 0; i--)
	{
		ind = nm1 - i;
		index = i * (n2p1 - i) / 2;
		indexi = i * (n2m1 - i) / 2;
		coeff = -d[index] * c[2 * i + 1];
		for (j = i + 1; j < n; j++)
		{
			indexj = indexi + j;
			d[indexj] = coeff * d[indexj + ind];
		}
	}

	/* calculation of the inverse matrix (symmetric matrix) */
	for (i = 0; i < n; i++)
	{
		index = i * (n2m1 - i) / 2;
		for (j = i; j < n; j++)
		{
			indexj = j * (n2m1 - j) / 2;
			sumb = 0.0l;
			for (k = j; k < n; k++)
				sumb += d[index + k] * d[indexj + k];
			b[index + j] = sumb;
		}
	}

	free(c);
	free(d);
}
