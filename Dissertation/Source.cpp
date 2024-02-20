// Diplom.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// REWENIE SKALYARNOY ZADA4i DIFRAKCII NA SYSTEME AKUSTI4ESKY MYAGKIH EKRANOV METODOM GALERKINA
// HOSPADE ZACHEM VSE ETO PROISHODIT...


//начать переводить в трехмерную задачу, параметры, все такое, начинать дипломную записку.
#include <fstream>
//#include "Complex.h"
#include <mpi.h>
#include "Complex (2).h"
#include <iostream>
#include <cmath>

using namespace std;

complex t;
double pi = 4.0 * atan(1.0), k0 = 2 * pi;
const int n = 4, N = 2 * n * n;  // n - число разбиения
//const double lymda = 0.5;
const double GranA1 = 0.0, GranA2 = 2.0;
const double GranB1 = 1.0, GranB2 = 3.0;
const double GranC1 = 0.0, GranC2 = 2.0;
const double GranD1 = 1.0, GranD2 = 3.0;
const double H11 = (GranB1 - GranA1) / n, H12 = (GranD1 - GranC1) / n;
const double H21 = (GranB2 - GranA2) / n, H22 = (GranD2 - GranC2) / n;
//double X11[n + 1], X12[n + 1];
//double X21[n + 1], X22[n + 1];
//complex A1[n * n][n * n + 1], C1[n * n];
//complex A2[n * n][n * n + 1], C2[n * n];

complex A[N][N + 1], C[N]; // !!!
double U1[n + 1], V1[n + 1];
double U2[n + 1], V2[n + 1];


void printcomplex(complex z) {
	printf("%6.3f,%6.3f|", real(z), imag(z));
	//cout <<"(" << z.real << "+" << z.imag<< ""<< ") ";
}


double X_Param(double u, double v, int num) {
	if (num)
	{
		//return cos(u)*cos(v);
		return u;
	}
	else
	{
		//return cos(u) * cos(v);
		return u; //пока u потом для другого типа экрана, например сфера, или кусок любой другой фигуры
	}

};
double Y_Param(double u, double v, int num) {
	if (num)
	{
		//return cos(u)*sin(v);
		return v;
	}
	else
	{
		//return cos(u) * sin(v);
		return v; //пока u потом для другого типа экрана, например сфера, или кусок любой другой фигуры
	}
}
double Z_Param(double u, double v, int num) {
	if (num)
	{
		//return sin(u);
		return 1;
	}
	else
	{
		//return cos(v);
		return 1; //пока u потом для другого типа экрана, например сфера, или кусок любой другой фигуры
	}
}


//тоже самое
double DXu_Param(double u, double v, int num) {
	return (X_Param(u + 0.0000001, v, num) - X_Param(u - 0.0000001, v, num)) / (2 * 0.0000001);
}
double DYu_Param(double u, double v, int num) {
	return (Y_Param(u + 0.0000001, v, num) - Y_Param(u - 0.0000001, v, num)) / (2 * 0.0000001);
}
double DZu_Param(double u, double v, int num) {
	return (Z_Param(u + 0.0000001, v, num) - Z_Param(u - 0.0000001, v, num)) / (2 * 0.0000001);
}
double DXv_Param(double u, double v, int num) {
	return (X_Param(u, v + 0.0000001, num) - X_Param(u, v - 0.0000001, num)) / (2 * 0.0000001);
}
double DYv_Param(double u, double v, int num) {
	return (Y_Param(u, v + 0.0000001, num) - Y_Param(u, v - 0.0000001, num)) / (2 * 0.0000001);
}
double DZv_Param(double u, double v, int num) {
	return (Z_Param(u, v + 0.0000001, num) - Z_Param(u, v - 0.0000001, num)) / (2 * 0.0000001);
}

// так же иф есле от номера экрана.
double sqrtEGF2(double u, double v,/* double u2, double v2,*/ int num) {

	double E, G, F;
	E = DXu_Param(u, v, num) * DXu_Param(u, v, num) + DYu_Param(u, v, num) * DYu_Param(u, v, num) + DZu_Param(u, v, num) * DZu_Param(u, v, num);
	G = DYv_Param(u, v, num) * DYv_Param(u, v, num) + DYv_Param(u, v, num) * DYv_Param(u, v, num) + DZv_Param(u, v, num) * DZv_Param(u, v, num);
	F = DXu_Param(u, v, num) * DXv_Param(u, v, num) + DYu_Param(u, v, num) * DYv_Param(u, v, num) + DZu_Param(u, v, num) * DZv_Param(u, v, num);
	//cout « E « " " « G « " " « F « " " « endl;
	return sqrt(E * G - F * F);
}


//ядро
complex Ker(double u1, double v1, double u2, double v2, int num1, int num2) {   //добавить z1, z2, тк теперь будет объемное тело
	//return(_i * (x1 - y2));
	//double rho = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	double x1 = X_Param(u1, v1, num1);
	double y1 = Y_Param(u1, v1, num1);
	double z1 = Z_Param(u1, v1, num1);
	double x2 = X_Param(u2, v2, num2);
	double y2 = Y_Param(u2, v2, num2);
	double z2 = Z_Param(u2, v2, num2);
	double rho = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
	if (rho > 1e-5)
	{
		return  exp(_i * k0 * rho) / (4.0 * pi * rho);
	}
	else
	{
		return 0.0 * _i;
	}

};
// ядро это функция грина. 

// это ядро для вычисления поля вне тела
complex KerVneEc(double x1, double y1, double z1, double u2, double v2, int num2) {
	//return(_i * (x1 - y2));
	//double rho = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

	double x2 = X_Param(u2, v2, num2);
	double y2 = Y_Param(u2, v2, num2);
	double z2 = Z_Param(u2, v2, num2);
	double rho = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
	if (rho > 1e-5)
	{
		return  exp(_i * k0 * rho) / (4.0 * pi * rho);
	}
	else
	{
		return 0.0 * _i;
	}

}
// правая часть
complex U0(double u1, double v1, int num1) {
	//return (x1 * x2) - (_i * lymda * (3.0 * x1 - 2.0)) / 12.0;
	//return 1 - lymda * _i * (x1 - 0.5);
	//double x1 = X_Param(u1, v1, num1);
	//return exp(_i * k0 * x1); //плоская волна
	return _i; //вначале с красивым этим
}
complex Ux(double x1, double x2) {
	complex ux(x1 * x2, 0);

	return ux;
}
complex Uy(double y1, double y2) {
	complex uy(y1 * y2, 0);

	return uy;
}




//pravilno?
//double phi(double xi, int i) { //poka odnomernoe potom peredelat nado na 2 mernoe
//	return((xi >= X11[i]) && (xi < X11[i + 1]));
//}
//

double phi2(double xi1, double xi2, int i, int j, int num) {
	if (num)
	{
		return((xi1 >= U1[i]) && (xi1 < U1[i + 1]) && (xi2 >= V1[j]) && (xi2 < V1[j + 1]));
	}
	else {
		return((xi1 >= U2[i]) && (xi1 < U2[i + 1]) && (xi2 >= V2[j]) && (xi2 < V2[j + 1]));
	}
}

// num = 1 если первая матрица, если нет то вторая
complex Integral_voln(int i1, int j1, int num) {
	double aa1, bb1, cc1, dd1;
	if (num)
	{
		aa1 = GranA1 + i1 * H11; bb1 = aa1 + H11; cc1 = GranC1 + j1 * H12; dd1 = cc1 + H12;
	}
	else {
		aa1 = GranA2 + i1 * H21; bb1 = aa1 + H21; cc1 = GranC2 + j1 * H22; dd1 = cc1 + H22;
	};

	double nn = 4, h1, h2, t1, t2;
	complex in(0.0, 0.0);
	h1 = (bb1 - aa1) / nn;
	h2 = (dd1 - cc1) / nn;

	for (int i1 = 0; i1 < nn; i1++) {
		for (int i2 = 0; i2 < nn; i2++) {
			t1 = aa1 + (i1 + 0.5) * h1;
			t2 = cc1 + (i2 + 0.5) * h2;

			in = in + U0(t1, t2, num) * sqrtEGF2(t1, t2, num);
		}
	}
	// in = in * h1 * h2;
	return in * h1 * h2;
}

//complex del(/*double xi,*/ int I, int J) { //vmesto etogo integral ot * itoi and jtoi basisnoy function
//	return complex(I == J, 0.0);
//}
//
////правильно
//complex del2(int I1, int J1, int I2, int J2) {
//	return complex(I1 == I2 && J1 == J2, 0.0);
//}


complex un(double x1, double x2, int num) {
	complex s(0.0, 0.0);
	int i = 0;
	for (int i1 = 0; i1 < n; i1++)
	{
		for (int j1 = 0; j1 < n; j1++) {
			i = i1 + n * j1;
			if (num)
			{
				s = s + C[i] * phi2(x1, x2, i1, j1, num); // or vice versa ... i%n, i/n)
			}
			else {
				s = s + C[i + n * n] * phi2(x1, x2, i1, j1, num); // or vice versa ... i%n, i/n)
			}
		}
	}
	return(s);
}


// num = 1 если первая матрица, если нет то вторая
complex Integral_ecran(int i1, int j1, int i2, int j2, int num1, int num2) {
	double aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2;
	if (num1)
	{
		aa1 = GranA1 + i1 * H11; bb1 = aa1 + H11; cc1 = GranC1 + j1 * H12; dd1 = cc1 + H12;
	}
	else {
		aa1 = GranA2 + i1 * H21; bb1 = aa1 + H21; cc1 = GranC2 + j1 * H22; dd1 = cc1 + H22;
	};
	if (num2)
	{
		aa2 = GranA1 + i2 * H11; bb2 = aa2 + H11; cc2 = GranC1 + j2 * H12; dd2 = cc2 + H12;  //тут
	}
	else {
		aa2 = GranA2 + i2 * H21; bb2 = aa2 + H21; cc2 = GranC2 + j2 * H22; dd2 = cc2 + H22;  // и тут было bb2 = aa1 + H11  и bb2 = aa1 + H21 соответственно
	};																						// aa2 bb2 получались равны => шаг h21 = 0 ну и пошло поехало

	int nn = 8;
	double h11, h12, h21, h22, t11, t12, t21, t22, rho;
	complex in(0.0, 0.0);
	h11 = (bb1 - aa1) / nn;
	h12 = (dd1 - cc1) / nn;
	h21 = (bb2 - aa2) / nn;
	h22 = (dd2 - cc2) / nn;

	for (int kk = 0; kk < nn; kk++)
	{
		for (int ll = 0; ll < nn; ll++)
		{
			for (int ii = 0; ii < nn; ii++)
			{
				for (int jj = 0; jj < nn; jj++) {
					//!!!! phi2(double xi1, double xi2, int i, int j) * phi2(double xi1, double xi2, int i, int j)
					t11 = aa1 + (jj + 0.5) * h11;
					t12 = cc1 + (ll + 0.5) * h12;

					t21 = aa2 + (ii + 0.5) * h21;
					t22 = cc2 + (kk + 0.5) * h22;
					rho = sqrt((t11 - t21) * (t11 - t21) + (t12 - t22) * (t12 - t22));
					//cout << sqrtEGF2(t11, t12, num1) << endl;
					if (rho > 1e-7) in = in + Ker(t11, t12, t21, t22, num1, num2) * sqrtEGF2(t11, t12, num1) * sqrtEGF2(t21, t22, num2);
				}
			}
		}
	}

	return in * h11 * h12 * h21 * h22;
}

//нужен новый мидлпрям с двойным интегралом для вычисления поля вне экрана.
complex Integral_ecran_VNE(double x, double y, double z, int num) { //double x,y,z, int num
	double aa, aa2, bb, bb2, cc, cc2, dd, dd2;
	if (num == 1)
	{
		aa = GranA1; bb = GranB1; cc = GranC1; dd = GranD1;
	}
	else {
		aa = GranA2; bb = GranB2; cc = GranC2; dd = GranD2;
	};

	int nn = 8;
	double h1, h2, t1, t2, rho;
	complex in(0.0, 0.0);
	h1 = (bb - aa) / nn;
	h2 = (dd - cc) / nn;

	for (int kk = 0; kk < nn; kk++)
	{
		//for (int ll = 0; ll < nn; ll++)
		//{
		for (int ii = 0; ii < nn; ii++)
		{

			t1 = aa + (ii + 0.5) * h1;
			t2 = cc + (kk + 0.5) * h2;
			in = in + KerVneEc(x, y, z, t1, t2, num) * sqrtEGF2(t1, t2, num) * un(t1, t2, num);
			//}  // need !!!!  K
		}
		//}
	}

	return in /** h11 * h12*/ * h1 * h2;
}



void Gauss(int k, complex Matrix[N][N + 1]) {
	if (Matrix[k][k] != complex(1.0, 0.0)) {
		complex T = Matrix[k][k];
		for (int j = k; j < N + 1; j++) {
			Matrix[k][j] = Matrix[k][j] / T;
		}
	}
	for (int i = 0; i < N; i++) {
		if ((Matrix[i][k] != complex(0.0, 0.0)) && (i != k)) {
			complex T = Matrix[i][k];
			Matrix[i][k] = complex(0.0, 0.0);
			for (int j = k + 1; j < N + 1; j++) {
				Matrix[i][j] -= Matrix[k][j] * T;
			}
		}
	}
	if (k < N - 1) {
		Gauss(k + 1, Matrix);
	}
}


void print_un(int pn, int num) {
	// double localA1, localB1, localC1, localD1;
	double t1, t2;
	if (num) {
		printf("\n");
		for (int i1 = 0; i1 < pn; i1++) {
			for (int i2 = 0; i2 < pn; i2++) {
				t1 = GranA1 + (GranB1 - GranA1) / pn * i1;
				t2 = GranC1 + (GranD1 - GranC1) / pn * i2;
				//t2 = GranC2 + (GranD2 - GranC2) / pn * i2;
				printf("%5.5f ", abs(un(t1, t2, num)));
				/// File1 << abs(un(t1, t2)) << " ";
			}
			printf("\n");
			// File1 << std::endl;
		}
	}
	else {
		printf("\n");
		for (int i1 = 0; i1 < pn; i1++) {
			for (int i2 = 0; i2 < pn; i2++) {
				t1 = GranA2 + (GranB2 - GranA2) / pn * i1;
				t2 = GranC2 + (GranD2 - GranC2) / pn * i2;
				printf("%5.5f ", abs(un(t1, t2, num)));
				/// File1 << abs(un(t1, t2)) << " ";
			}
			printf("\n");
			// File1 << std::endl;
		}
	}

}

//void Zapis_v_File(int pn, int f) {
//	//double localA1, localB1, localC1, localD1;
//	if (f) {
//		std::ofstream File1("./Matrix_1.txt");
//		FILE* tab_file;
//		fopen_s(&tab_file, "result1.xls", "w");
//		double t1, t2;
//		// printf("\n");
//		for (int i1 = 0; i1 < pn; i1++) {
//			for (int i2 = 0; i2 < pn; i2++) {
//				t1 = GranA1 + (GranB1 - GranA1) / pn * i1;
//				t2 = GranC1 + (GranD1 - GranC1) / pn * i2;
//				// printf("%6.3f ", abs(un(t1, t2)));
//				File1 << XP[kk] << " " << YP[kk] << " " << ZP[kk] << " " << abs(un(t1, t2, f)) << "\n";
//				kk--;
//				fprintf(tab_file, "%5.5f\t", abs(un(t1, t2, f)));
//			}
//			//printf("\n");
//			File1 << "\n";
//			fprintf(tab_file, "\n");
//		}
//		File1.close();
//		fclose(tab_file);
//	}
//	else {
//		std::ofstream File1("./Matrix_2.txt");
//		FILE* tab_file;
//		fopen_s(&tab_file, "result2.xls", "w");
//		double t1, t2;
//		// printf("\n");
//		for (int i1 = 0; i1 < pn; i1++) {
//			for (int i2 = 0; i2 < pn; i2++) {
//				t1 = GranA2 + (GranB2 - GranA2) / pn * i1;
//				t2 = GranC2 + (GranD2 - GranC2) / pn * i2;
//				// printf("%6.3f ", abs(un(t1, t2)));
//				File1 << abs(un(t1, t2, f)) << "\t";
//				fprintf(tab_file, "%5.5f\t", abs(un(t1, t2, f)));
//			}
//			//printf("\n");
//			File1 << "\n";
//			fprintf(tab_file, "\n");
//		}
//		File1.close();
//		fclose(tab_file);
//	}
//}

//запись в файл для графика в Visit
void Zapis_v_File_Visit(int pn) {
	//std::ofstream File3("../Matrix_3.txt");
	FILE* tab_file1;
	fopen_s(&tab_file1, "resultVIZIT.txt", "w");
	//int pn = 50;
	double t1, t2;
	// printf("\n");
	for (int i1 = 0; i1 < pn; i1++) {
		for (int i2 = 0; i2 < pn; i2++) {
			t1 = GranA2 + (GranB2 - GranA2) / pn * i1;
			t2 = GranC2 + (GranD2 - GranC2) / pn * i2;
			// printf("%6.3f ", abs(un(t1, t2)));
			//File3 << X_Param(t1, t2, 0) << " " << Y_Param(t1, t2, 0) << " " << Z_Param(t1, t2, 0) << " " << abs(un(t1, t2, 0)) << "\n";
			fprintf(tab_file1, "%5.5f\t%5.5f\t%5.5f\t%5.5f\n", X_Param(t1, t2, 0), Y_Param(t1, t2, 0), Z_Param(t1, t2, 0), abs(un(t1, t2, 0)));
			t1 = GranA1 + (GranB1 - GranA1) / pn * i1;
			t2 = GranC1 + (GranD1 - GranC1) / pn * i2;
			// printf("%6.3f ", abs(un(t1, t2)));
			//File3 << X_Param(t1, t2, 1) << " " << Y_Param(t1, t2, 1) << " " << Z_Param(t1, t2, 1) << " " << abs(un(t1, t2, 1)) << "\n";
			fprintf(tab_file1, "%5.5f\t%5.5f\t%5.5f\t%5.5f\n", X_Param(t1, t2, 1), Y_Param(t1, t2, 1), Z_Param(t1, t2, 1), abs(un(t1, t2, 1)));

		}

	}
	fclose(tab_file1);
	//File3.close();

}
void Zapis_v_File_Visit_VNE() {
	FILE* tab_file2;
	fopen_s(&tab_file2, "resultVIZIT_VNE.txt", "w");
	double x_vne = 0.0, y_vne = 0.0, z_vne = 1.0;
	double field_vne;
	for (int i = -150; i < 150; i++)
	{
		x_vne = i * 0.02;
		for (int j = -150; j < 150; j++) {
			y_vne = j * 0.02;

			// в плоскости x,y
			//field_vne = abs(Integral_ecran_VNE(x_vne, y_vne, z_vne, 0) + Integral_ecran_VNE(x_vne, y_vne, z_vne, 1)); // !!!!!!
			//fprintf(tab_file2, "%5.5f\t%5.5f\t%5.5f\t%5.5f\t\n", x_vne, y_vne, z_vne, field_vne); //!!!!!!
			// в плоскости y,z
			field_vne = abs(Integral_ecran_VNE(x_vne, x_vne, y_vne, 0) + Integral_ecran_VNE(x_vne, x_vne, y_vne, 1)); // !!!!!!
			//найти максимум field_vne и отсеять 10% верхних значений
			fprintf(tab_file2, "%5.5f\t%5.5f\t%5.5f\t%5.5f\t\n", x_vne, x_vne, y_vne, field_vne); //!!!!!!

		}
	}
	fclose(tab_file2);
}


complex buff1[N][N + 1];
double buff2[N][N + 1];
int main() {
	//cout << " AAAAAAAAAAAAAA";

	//double xi2[n + 1], xi1[n + 1];
	MPI_Init(NULL, NULL);
	std::cout << " ya der'mo1" << std::endl;

	double starttime, endtime;
	starttime = MPI_Wtime();
	double starttimeZ1 = MPI_Wtime();
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	printf_s("rank = %i, world_size = %i", world_rank, world_size);
	std::fflush(stdout);


	for (int j = 0; j < n + 1; j++) {
		U1[j] = GranA1 + j * H11;
		V1[j] = GranC1 + j * H12;
		//cout << U1[j] << "   " << V1[j] << endl;
	}
	for (int j = 0; j < n + 1; j++) {
		U2[j] = GranA2 + j * H21;
		V2[j] = GranC2 + j * H22;
		//cout << U2[j] << "   " << V2[j] << endl;
	}
	//if (world_rank == 0)
	//{
	//	cout << endl;
	//	complex o = Integral_ecran(0, 0, 0, 0, 1, 1);
	//	printf("%g+%gi \n", real(o), imag(o));
	//	complex o1 = Integral_ecran(0, 0, 0, 2, 1, 1); //только его
	//	printf("%g+%gi \n", real(o1), imag(o1));
	//	complex o2 = Integral_ecran(0, 0, 0, 2, 1, 0);
	//	printf("%g+%gi \n", real(o2), imag(o2));
	//	complex p = Integral_voln(0, 0, 1);
	//	printf("%g+%gi \n", real(p), imag(p));
	//	//complex p1 = Integral_voln(3, 3, 0);
	//	//printf("re = %g, im = %g \n", real(p1), imag(p1));
	//	//printcomplex(o);
	//}
	//std::cout << " ----------------------------------------------------- " << std::endl;
	int i, j;
	std::fflush(stdout);
	for (int i1 = 0; i1 < n; i1++)
	{
		for (int j1 = 0; j1 < n; j1++)
		{
			i = j1 + n * i1;
			for (int i2 = 0; i2 < n; i2++)
			{
				for (int j2 = 0; j2 < n; j2++)
				{

					j = j2 + n * i2;
					if (world_rank == 0)
					{
						A[i][j] = Integral_ecran(i1, j1, i2, j2, 1, 1);
					}
					if (world_rank == 1)
					{
						A[i][j + n * n] = Integral_ecran(i1, j1, i2, j2, 1, 0);
					}
					if (world_rank == 2)
					{
						A[i + n * n][j] = Integral_ecran(i1, j1, i2, j2, 0, 1);
					}
					if (world_rank == 3)
					{
						A[i + n * n][j + n * n] = Integral_ecran(i1, j1, i2, j2, 0, 0);
					}



				}


			}
			if (world_rank == 0)
			{
				A[i][N] = Integral_voln(i1, j1, 1);
			}
			if (world_rank == 1)
			{
				A[i + n * n][N] = Integral_voln(i1, j1, 0);
			}

		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(A, buff1, 2 * N * (N + 1) /*/ row_size*/, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (world_rank == 0)
	{
		double endtime1 = MPI_Wtime();
		std::printf("samo rasparallel zanyalo %f seconds\n", endtime1 - starttimeZ1);
	}
	std::fflush(stdout);
	if (world_rank == 0)
	{
		double starttimeG = MPI_Wtime();
		Gauss(0, buff1);
		double endtimeG = MPI_Wtime();
		std::printf("Gauss zanyal %f seconds\n", endtimeG - starttimeG);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//Gauss(0, A2);
	//cout << " ----------------------------------------------------- " << endl;
	if (world_rank == 0) {

		double starttimeZ = MPI_Wtime();
		for (int i = 0; i < N; i++)
		{
			C[i] = buff1[i][N];

		}


		std::cout << " ----------------------------------------------------- " << std::endl;
		//print_un(n, 1);
		//Zapis_v_File(n, 1);
		//print_un(n, 0);
		Zapis_v_File_Visit(50);
		/*Zapis_v_File_Visit_VNE();*/

		endtime = MPI_Wtime();
		std::printf("Zapis zanyala %f seconds\n", endtime - starttimeZ);
		std::printf("Obwee vipolnenie zanyalo %f seconds\n", endtime - starttime);

	}
	//MPI_Barrier(MPI_COMM_WORLD);

	//MPI_Comm_free(&row_comm);
	MPI_Finalize();


	return 0;
}