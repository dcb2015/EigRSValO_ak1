// EigRSValO_ak1 - Program for calculating the Eigenvalues ONLY of a N X N real, symmetric matrix.
//
// This program is a severely edited translation of the EISPACK routine RS.F.
//
// References:
//
// Smith, B.T.; J.M. Boyle; J.J. Dongarra; B.S. Garbow; Y. Ikebe; V.C. Klema; and C.B. Moler.
//      "Matrix Eigensystem Routines--(EISPACK) Guide"
//      Springer - Verlag, Berlin.
//      1976
//
// Garbow, B.S.; J.M. Boyle; J.J. Dongarra; and C.B. Moler.
//      "Matrix Eigensystem Routines--(EISPACK) Guide Extension"
//      Springer - Verlag, Berlin.
//      1977
//
// To distinguish this version of RS from other translations,
// an '_ak1' suffix has been appended to its name.
//
// The main program accepts a full matrix and tests for symmetry. However, the basic algorithm uses
// only the lower triangular part of the matrix. The strictly upper triangular part of the matrix
// is not used. In other words, the upper triangular part of the matrix could be filled with random
// numbers and the test for symmetry omitted. Also, the matrix is changed by this program so, 
// if it needs to be re-used elsewhere, it should be copied and a copy passed into the routine.
//
// 5 June 2018
//
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
//

#include <iostream>
#include <fstream>
#include <cctype>
#include <cmath>
#include <vector>
#include <cfloat>

using namespace std;

typedef vector<double> C1DArray;
typedef vector<vector<double>> C2DArray;

double pythag(double a, double b);
void Echo2DArray(const C2DArray a2);

double pythag(double a, double b)
{ // Returns the square root of (a*a + b*b)
  // without overflow or destructive underflow

	double p, r, s, t, u;

	p = t = fabs(a);
	r = u = fabs(b);
	if (u > t) {
		p = u;
		r = t;
	}

	if (p > 0){
		r /= p;
		r *= r;
		t = 4.0 + r;

		while (t > 4.0){
			s = r / t;
			u = 1.0 + 2.0*s;
			p *= u;
			t = s / u;
			r *= t*t;
			t = 4.0 + r;
		} // while (t > 4.0);
	} // End if (p > 0)

	return p;
} // End pythag

void eigRSvalo_ak1(const int N, C2DArray& A_Matrix, C1DArray& fv1_vec, C1DArray& wr_vec, int *info);

void eigRSvalo_ak1(const int N, C2DArray& A_Matrix, C1DArray& fv1_vec, C1DArray& wr_vec, int *info){

	int i = N - 1, ii = i, j, k, l = i, l1;
	double c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, scale, tst1, tst2;

	// ======BEGINNING OF TRED1 ===================================

	wr_vec[i] = A_Matrix[i][i];

	while (i > 0){
		--i;
		wr_vec[i] = A_Matrix[ii][i];
		A_Matrix[ii][i] = A_Matrix[i][i];
	} // End while loop

	for (i = l; i > 0; --i){

		--l;
		scale = h = 0.0;

		for (j = l; j >= 0; --j)  scale += fabs(wr_vec[j]);

		if (scale == 0.0){

			for (j = 0; j <= l; ++j){
				wr_vec[j] = A_Matrix[l][j];
				A_Matrix[l][j] = A_Matrix[i][j];
				A_Matrix[i][j] = 0.0;
			}//End for j

			fv1_vec[i] = 0.0;
			continue;

		} // End if (scale == 0.0)

		for (j = l; j >= 0; --j){
			wr_vec[j] /= scale;
			h += wr_vec[j] * wr_vec[j];
		}//End for j

		f = wr_vec[l];
		g = ((f >= 0) ? -sqrt(h) : sqrt(h));
		fv1_vec[i] = g*scale;
		h -= f*g;
		wr_vec[l] = f - g;

		if (l != 0){

			for (j = l; j >= 0; --j) fv1_vec[j] = 0.0;

			for (j = 0; j <= l; ++j){
				f = wr_vec[j];
				g = fv1_vec[j] + f*A_Matrix[j][j];
				for (ii = (j + 1); ii <= l; ++ii){
					g += wr_vec[ii] * A_Matrix[ii][j];
					fv1_vec[ii] += f*A_Matrix[ii][j];
				} // End for ii
				fv1_vec[j] = g;
			}//End for j

			// Form p

			f = 0.0;
			for (j = l; j >= 0; --j){
				fv1_vec[j] /= h;
				f += fv1_vec[j] * wr_vec[j];
			}//End for j

			h = f / h * 0.5;

			// Form q

			for (j = l; j >= 0; --j)  fv1_vec[j] -= h*wr_vec[j];

			// Form reduced A

			for (j = 0; j <= l; ++j){

				f = wr_vec[j];
				g = fv1_vec[j];

				for (ii = j; ii <= l; ++ii)	 A_Matrix[ii][j] -= f*fv1_vec[ii] + g*wr_vec[ii];

			}//End for j

		} // End if (l != 0)

		for (j = 0; j <= l; ++j){
			f = wr_vec[j];
			wr_vec[j] = A_Matrix[l][j];
			A_Matrix[l][j] = A_Matrix[i][j];
			A_Matrix[i][j] = f*scale;
		}//End for j

	}//End for i

	fv1_vec[0] = 0.0;

	// ======END OF TRED1 =========================================

	// ======BEGINNING OF TQL1 ===================================

	*info = -1;

	if (N == 1)  return;

	for (i = 1; i < N; ++i)  fv1_vec[i - 1] = fv1_vec[i];

	fv1_vec[N - 1] = tst1 = f = 0.0;
	l1 = 0;

	for (l = 0; l < N; ++l){
		++l1;   // l1 = l + 1 in this for l loop
		j = 0;
		h = fabs(wr_vec[l]) + fabs(fv1_vec[l]);

		if (tst1 < h) tst1 = h;

		// Look for small sub-diagonal element

		for (k = l; k < N; ++k){
			tst2 = tst1 + fabs(fv1_vec[k]);
			if (tst2 == tst1) break; // fv1[N-1] is always 0, so there is no exit through the bottom of the loop
		}//End for k

		if (k != l){

			do {

				if (j == 30){
					*info = l;
					return;
				} // End if (j == 30)

				++j;

				// Form shift

				g = wr_vec[l];
				p = (wr_vec[l1] - g) / (2.0*fv1_vec[l]);
				scale = r = pythag(p, 1.0);	// Use scale as a dummy variable
				if (p < 0) scale = -scale;
				scale += p;
				wr_vec[l] = fv1_vec[l] / scale;
				dl1 = wr_vec[l1] = fv1_vec[l] * scale;
				h = g - wr_vec[l];

				for (i = l1 + 1; i < N; ++i)  wr_vec[i] -= h;

				f += h;

				// q1 transformation

				p = wr_vec[k];
				c2 = c = 1.0;
				el1 = fv1_vec[l1];
				s = 0.0;

				// Look for i = k - 1 until l in steps of -1

				for (i = (k - 1); i >= l; --i){
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c*fv1_vec[i];
					h = c*p;
					r = pythag(p, fv1_vec[i]);
					fv1_vec[i + 1] = s*r;
					s = fv1_vec[i] / r;
					c = p / r;
					p = c*wr_vec[i] - s*g;
					wr_vec[i + 1] = h + s*(c*g + s*wr_vec[i]);
				}//End for i

				p = -s*s2*c3*el1*fv1_vec[l] / dl1;
				fv1_vec[l] = s*p;
				wr_vec[l] = c*p;
				tst2 = tst1 + fabs(fv1_vec[l]);
			} while (tst2 > tst1);  // End do-while loop

		} // End if (k != l)

		p = wr_vec[l] + f;

		// Order eigenvalues

		// For i = l to 1, in steps of -1
		for (i = l; i >= 1; --i){
			if (p >= wr_vec[i - 1])  break;
			wr_vec[i] = wr_vec[i - 1];
		}//End for i

		wr_vec[i] = p;

	}//End for l

	// ======END OF TQL1 =========================================

	return;
} // End eigRSvalo_ak1

void Echo2DArray(const C2DArray a2) {
	//Routine to output a 2D Array of type double to the console
	cout << "\n";
	for (unsigned int i = 0; i < a2.size(); ++i) {
		for (unsigned int j = 0; j < a2[0].size(); ++j) {
			cout << a2[i][j] << " ";
		} // End for j
		cout << "\n";
	} // End for i
	cout << "\n";
	return;
} // End Echo2DArray

int main()
{
	char rflag = 0;	//Readiness flag

	cout << "                     EigRSValO_ak1   (5 June 2018)\n";
	cout << "=========================================================================== \n";
	cout << "This program calculates the eigenvalues ONLY of a N X N, real,\n";
	cout << "symmetric, Matrix, A.\n";
	cout << "\nThe A Matrix to be input should have been saved beforehand in a file named\n";
	cout << "EigSysRS.txt, which should be in the same folder as the EigRSValO executable.\n";
	cout << "The first entry in this file should be N, the size of the N X N matrix.\n";
	cout << "The entries for A should follow, with data for row 1 first, then row 2,\n";
	cout << "then row 3, etc.\n";
	cout << "\nThe data is assumed to be of type double. Variables used within this program\n";
	cout << "are type double.\n";
	cout << "\nOutput--eigenvalues--is written to the file EigOutRS.txt.\n";

	cout << "Note the Error Code output.\n";
	cout << "If normal return, ierr = 0. If Error Code > 0,\n";
	cout << "it indicates that more than 30 iterations of a subroutine were required to\n";
	cout << "determine an eigenvalue. In this case, the subroutine terminated.\n";
	cout << "Error Code gives the index of the eigenvalue for which the failure occurred.\n";
	cout << "Eigenvalue[1], Eigenvalue[2], . . . Eigenvalue[ErCode - 1] should be correct.\n";
	cout << "\nIs everything ready (are you ready to continue?)? If yes, Enter y. \n";
	cout << "Otherwise Enter any other key. \n";
	cin >> rflag;

	if (toupper(rflag) == 'Y') {
		C2DArray A;				// A is a 2D matrix of type double
		C1DArray fv1;			// Temporary vector
		C1DArray wr;			// Vector for eigenvalues
		int i, ierr, j, mDim;
		bool erFlag = false;  // Error flag

		ifstream in("EigSysRS.txt", ios::in);

		if (!in) {
			cout << "Cannot open the input file." << endl;
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;
		}

		in >> mDim;  //Input the Matrix dimension from the file
		if (mDim < 1)  {
			in.close(); //Close the input file before terminating
			cout << "\nInvalid dimension entered. Program terminated." << endl;
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;
		}

		ofstream out("EigOutRS.txt", ios::out);

		if (!out) {
			in.close(); //Close the input file before terminating
			cout << "\nCannot open the output file. Program terminated." << endl;
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;
		}

		// Beginning of try block, if vector re-sizing unsuccessful
		try {

			// Resize the arrays to the appropriate sizes

			A.resize(mDim); // N rows
			fv1.resize(mDim); // N rows
			wr.resize(mDim); // N rows

			for (i = 0; i < mDim; ++i)  A[i].resize(mDim); // N columns

		} // End of try block

		catch (bad_alloc& xa) { // Catch block, for exceptions

			in.close();
			out.close();
			cerr << "\nIn catch block, so an exception occurred: " << xa.what() << endl;
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;

		} // End of catch block

		for (i = 0; i < mDim; ++i){ //Input the A Matrix from the file
			for (j = 0; j < mDim; ++j) in >> A[i][j];
		}//End for i

		in.close();  //Close the input file

		Echo2DArray(A);

		// Confirm the symmetry of the matrix
		i = mDim - 1;
		while ((i > 0) && (!erFlag)){

			j = i;
			do {
				--j;
				if (fabs(A[i][j] - A[j][i]) > DBL_EPSILON){
					erFlag = true;
					break;
				} // End if
			} while (j > 0); // End do-while

			--i;

		} // End while (i > 0)

		if (erFlag) {
			cout << "Non-symmetry in matrix detected. No further action taken. Program terminated." << endl << endl;
			out.close();
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;
		} // End if erFlag

		cout << "Matrix seems to be symmetric." << endl << endl;

		eigRSvalo_ak1(mDim, A, fv1, wr, &ierr);

		out.precision(DBL_DIG);

		out << "ierr = " << (ierr + 1) << "\n";
		out << "The eigenvalues are:\n\n";

		for (i = 0; i < mDim; ++i)  out << wr[i] << " \n";

		out << "\n";

		out.close();
		cout << "\nDone! The solution is in the text file EigOutRS.txt \n";

	} //End if rflag = 'Y'
	else cout << "\nNot ready. Try again when ready with information. \n";
	cout << "\nEnter any key to continue. \n";
	cin >> rflag;
	return 0;
}                      // End main program.
