#include <vector>
#include <math.h>
#include <iostream>
#include "bresenham.h"

template<typename T>
class Vector2D {
		T* vec;
		int l;
		int c;
	public:
		Vector2D(int lin, int col) {
			vec = new T[lin * col];
			this->l = lin;
			this->c = col;
		}
		virtual ~Vector2D() {
			delete[] vec;
		}

		T& operator()(int lin, int col) {
			return vec[lin * c + col];
		}
};

class ControlPoint {
	public:
		int x;
		int y;

		ControlPoint(int a, int b) {
			x = a;
			y = b;
		}

		ControlPoint() {
		}
};

/** contains the points to be drawn - T is a class representing a point */
template<typename T>
class DrawContainer {
	private:
		std::vector<T> points;
	public:
		DrawContainer() {

		}

		void operator()(int x, int y) {
			points.push_back(T(x, y));
		}

		T& operator[](int idx) {
			return points[idx];
		}

		int size() {
			return points.size();
		}

		std::vector<T>& getPoints() {
			return points;
		}
};

class CurvePanel {
	private:
		int i;
		ConicPlotter<DrawContainer<ControlPoint> > plotter;
		std::vector<ControlPoint> points;
	public:
		CurvePanel() {

		}

		void clearPoints() {
			points.clear();
		}

		void addPoint(int x, int y) {
			points.push_back(ControlPoint(x, y));
		}

	private:
		void ROTATE(Vector2D<double>& a, int i, int j, int k, int l,
				double tau, double s) {
			double g, h;
			g = a(i, j);
			h = a(k, l);
			a(i, j) = g - s * (h + g * tau);
			a(k, l) = h + s * (g - h * tau);
		}

		void jacobi(Vector2D<double>& a, int n, double d[],
				Vector2D<double>& v, int nrot) {
			int j, iq, ip, i;
			double tresh, theta, tau, t, sm, s, h, g, c;

			double b[n];
			double z[n];

			for (ip = 0; ip < n; ip++) {
				for (iq = 0; iq < n; iq++)
					v(ip, iq) = 0.0;
				v(ip, ip) = 1.0;
			}
			for (ip = 0; ip < n; ip++) {
				b[ip] = d[ip] = a(ip, ip);
				z[ip] = 0.0;
			}
			nrot = 0;
			for (i = 1; i <= 50; i++) {
				sm = 0.0;
				for (ip = 0; ip < n - 1; ip++) {
					for (iq = ip + 1; iq < n; iq++)
						sm += fabs(a(ip, iq));
				}
				if (sm == 0.0) {
					return;
				}
				if (i < 4)
					tresh = 0.2 * sm / (n * n);
				else
					tresh = 0.0;
				for (ip = 0; ip < n - 1; ip++) {
					for (iq = ip + 1; iq < n; iq++) {
						g = 100.0 * fabs(a(ip, iq));
						if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) && fabs(
								d[iq]) + g == fabs(d[iq]))
							a(ip, iq) = 0.0;
						else if (fabs(a(ip, iq)) > tresh) {
							h = d[iq] - d[ip];
							if (fabs(h) + g == fabs(h))
								t = (a(ip, iq)) / h;
							else {
								theta = 0.5 * h / (a(ip, iq));
								t = 1.0 / (fabs(theta) + sqrt(1.0 + theta
										* theta));
								if (theta < 0.0)
									t = -t;
							}
							c = 1.0 / sqrt(1 + t * t);
							s = t * c;
							tau = s / (1.0 + c);
							h = t * a(ip, iq);
							z[ip] -= h;
							z[iq] += h;
							d[ip] -= h;
							d[iq] += h;
							a(ip, iq) = 0.0;
							for (j = 0; j < ip - 1; j++) {
								ROTATE(a, j, ip, j, iq, tau, s);
							}
							for (j = ip; j < iq - 1; j++) {
								ROTATE(a, ip, j, j, iq, tau, s);
							}
							for (j = iq; j < n; j++) {
								ROTATE(a, ip, j, iq, j, tau, s);
							}
							for (j = 0; j < n; j++) {
								ROTATE(v, j, ip, j, iq, tau, s);
							}
							++nrot;
						}
					}
				}
				for (ip = 0; ip < n; ip++) {
					b[ip] += z[ip];
					d[ip] = b[ip];
					z[ip] = 0.0;
				}
			}
		}

		/* Perform the Cholesky decomposition
		 * Return the lower triangular L such that L*L'=A
		 */
		void choldc(Vector2D<double>& a, int n, Vector2D<double>& l) {
			int i, j, k;
			double sum;
			double p[n];

			for (i = 0; i < n; i++) {
				for (j = i; j < n; j++) {
					for (sum = a(i, j), k = i - 1; k >= 0; k--)
						sum -= a(i, k) * a(j, k);
					if (i == j) {
						if (sum > 0.0)
							p[i] = sqrt(sum);
					} else {
						a(j, i) = sum / p[i];
					}
				}
			}
			for (i = 0; i < n; i++)
				for (j = i; j < n; j++)
					if (i == j)
						l(i, i) = p[i];
					else {
						l(j, i) = a(j, i);
						l(i, j) = 0.0;
					}
		}

		int inverse(Vector2D<double>& TB, Vector2D<double>& InvB, int N) {
			int k, i, j, p, q;
			double mult;
			double D, temp;
			double maxpivot;
			int npivot;
			Vector2D<double> B(N, N + 1);
			Vector2D<double> A(N, 2 * N + 1);
			Vector2D<double> C(N, N);
			double eps = 10e-20;

			for (k = 0; k < N; k++)
				for (j = 0; j < N; j++)
					B(k, j) = TB(k, j);

			for (k = 0; k < N; k++) {
				for (j = 0; j < N + 1; j++)
					A(k, j) = B(k, j);
				for (j = N + 1; j < 2 * N + 1; j++)
					A(k, j) = (float) 0;
				A(k, k - 1 + N + 2) = (float) 1;
			}
			for (k = 0; k < N; k++) {
				maxpivot = fabs((double) A(k, k));
				npivot = k;
				for (i = k; i < N; i++)
					if (maxpivot < fabs((double) A(i, k))) {
						maxpivot = fabs((double) A(i, k));
						npivot = i;
					}
				if (maxpivot >= eps) {
					if (npivot != k)
						for (j = k; j < 2 * N + 1; j++) {
							temp = A(npivot, j);
							A(npivot, j) = A(k, j);
							A(k, j) = temp;
						};
					D = A(k, k);
					for (j = 2 * N; j > k; j--)
						A(k, j) = A(k, j) / D;
					for (i = 0; i < N; i++) {
						if (i != k) {
							mult = A(i, k);
							for (j = 2 * N; j > k; j--)
								A(i, j) = A(i, j) - mult * A(k, j);
						}
					}
				} else {
					return (-1);
				};
			}
			for (k = 0, p = 0; k < N; k++, p++)
				for (j = N + 1, q = 0; j < 2 * N + 1; j++, q++)
					InvB(p, q) = A(k, j);
			return (0);
		}

		void mulAB(Vector2D<double>& _A, Vector2D<double>& _B,
				Vector2D<double>& _res, int _righA, int _colA, int _righB,
				int _colB) {
			int p, q, l;
			for (p = 0; p < _righA; p++)
				for (q = 0; q < _colB; q++) {
					_res(p, q) = 0.0;
					for (l = 0; l < _colA; l++)
						_res(p, q) = _res(p, q) + _A(p, l) * _B(l, q);
				}
		}

		void mulA_TB(Vector2D<double>& _A, Vector2D<double>& _B, Vector2D<
				double>& _res, int _righA, int _colA, int _righB, int _colB) {
			int p, q, l;
			for (p = 0; p < _colA; p++)
				for (q = 0; q < _colB; q++) {
					_res(p, q) = 0.0;
					for (l = 0; l < _righA; l++)
						_res(p, q) = _res(p, q) + _A(l, p) * _B(l, q);
				}
		}

		void mulAB_T(Vector2D<double>& _A, Vector2D<double>& _B, Vector2D<
				double>& _res, int _righA, int _colA, int _righB, int _colB) {
			int p, q, l;
			for (p = 0; p < _colA; p++)
				for (q = 0; q < _colB; q++) {
					_res(p, q) = 0.0;
					for (l = 0; l < _righA; l++)
						_res(p, q) = _res(p, q) + _A(p, l) * _B(q, l);
				}
		}
	public:
		/** fits the points with an ellipse that has the axes alligned with Ox Oy */
		void fitStraightEllipse(double& el_a, double& el_b, double& el_c,
				double& el_d, double& el_e, double& el_f) {
			int np = points.size(); // number of points
			Vector2D<double> D(np, 5);
			Vector2D<double> S(5, 5);
			Vector2D<double> Const(5, 5);
			Vector2D<double> temp(5, 5);
			Vector2D<double> L(5, 5);
			Vector2D<double> C(5, 5);

			Vector2D<double> invL(5, 5);
			double d[5];
			Vector2D<double> V(5, 5);
			Vector2D<double> sol(5, 5);
			double tx, ty;
			int nrot = 0;

			double pvec[6];

			Const(0, 1) = -2;
			Const(1, 0) = -2;

			if (np < 6)
				return;

			/* generate the design matrix */
			for (int i = 0; i < np; i++) {
				tx = points[i].x;
				ty = points[i].y;

				D(i, 0) = tx * tx;
				D(i, 1) = ty * ty;
				D(i, 2) = tx;
				D(i, 3) = ty;
				D(i, 4) = 1.0;
			}

			/* compute scatter matrix */
			mulA_TB(D, D, S, np, 5, np, 5);

			/* get L lower triangular s.t. S = L*L' (conjugate transpose) */
			choldc(S, 5, L);

			/* inverse of L */
			inverse(L, invL, 5);

			mulAB_T(Const, invL, temp, 5, 5, 5, 5);
			mulAB(invL, temp, C, 5, 5, 5, 5);

			/* get d the eigenvalues and V the eigenvectors for C */
			jacobi(C, 5, d, V, nrot);

			/* get the solution */
			mulA_TB(invL, V, sol, 5, 5, 5, 5);

			/* normalize the solution */
			for (int j = 0; j < 5; j++) /* Scan columns */
			{
				double mod = 0.0;
				for (int i = 0; i < 5; i++)
					mod += sol(i, j) * sol(i, j);
				for (int i = 0; i < 5; i++)
					sol(i, j) /= sqrt(mod);
			}

			int solind = 0;

			/* get the nonnull eigenvalue */
			for (int i = 0; i < 5; i++)
				if (d[i] != 0)
					solind = i;

			/* fetch the solution corresponding to the nonzero eigenvalue */
			for (int j = 0; j < 5; j++)
				pvec[j] = sol(j, solind);

			el_a = pvec[0];
			el_b = 0;
			el_c = pvec[1];
			el_d = pvec[2];
			el_e = pvec[3];
			el_f = pvec[4];

			/*double delta = pvec[2] * pvec[2] + pvec[3] * pvec[3] * pvec[0]
			 / pvec[1] - 4 * pvec[0] * pvec[4];
			 xc = (int) (-pvec[2] / (2 * pvec[0]));
			 yc = (int) (-pvec[3] / (2 * pvec[1]));
			 el_a = (int) (sqrt(delta) / (2 * fabs(pvec[0])));
			 el_b = (int) (sqrt(delta) / (2 * sqrt(pvec[0] * pvec[1])));
			 */
		}

		/** fits the points with an ellipse that might be tilted at a certain angle */
		void fitTiltedEllipse(double& el_a, double& el_b, double& el_c,
				double& el_d, double& el_e, double& el_f) {
			int np = points.size(); // number of points
			Vector2D<double> D(np, 6);
			Vector2D<double> S(6, 6);
			Vector2D<double> Const(6, 6);
			Vector2D<double> temp(6, 6);
			Vector2D<double> L(6, 6);
			Vector2D<double> C(6, 6);

			Vector2D<double> invL(6, 6);
			double d[6];
			Vector2D<double> V(6, 6);
			Vector2D<double> sol(6, 6);
			double tx, ty;
			int nrot = 0;

			double pvec[6];

			Const(0, 2) = -2;
			Const(1, 1) = 1;
			Const(2, 0) = -2;

			if (np < 6)
				return;

			/* generate the design matrix */
			for (int i = 0; i < np; i++) {
				tx = points[i].x;
				ty = points[i].y;

				D(i, 0) = tx * tx;
				D(i, 1) = tx * ty;
				D(i, 2) = ty * ty;
				D(i, 3) = tx;
				D(i, 4) = ty;
				D(i, 5) = 1.0;
			}

			/* compute scatter matrix */
			mulA_TB(D, D, S, np, 6, np, 6);

			/* get L lower triangular s.t. S = L*L' (conjugate transpose) */
			choldc(S, 6, L);

			/* inverse of L */
			inverse(L, invL, 6);

			mulAB_T(Const, invL, temp, 6, 6, 6, 6);
			mulAB(invL, temp, C, 6, 6, 6, 6);

			/* get d the eigenvalues and V the eigenvectors for C */
			jacobi(C, 6, d, V, nrot);

			/* get the solution */
			mulA_TB(invL, V, sol, 6, 6, 6, 6);

			/* normalize the solution */
			for (int j = 0; j < 6; j++) /* Scan columns */
			{
				double mod = 0.0;
				for (int i = 0; i < 6; i++)
					mod += sol(i, j) * sol(i, j);
				for (int i = 0; i < 6; i++)
					sol(i, j) /= sqrt(mod);
			}

			int solind = 0;

			/* get the nonnull eigenvalue */
			for (int i = 0; i < 6; i++)
				if (d[i] != 0)
					solind = i;

			/* fetch the solution corresponding to the nonzero eigenvalue */
			for (int j = 0; j < 6; j++)
				pvec[j] = sol(j, solind);

			el_a = pvec[0];
			el_b = pvec[1];
			el_c = pvec[2];
			el_d = pvec[3];
			el_e = pvec[4];
			el_f = pvec[5];
		}

		/** finds the points the ellipse passes through
		 *	ellipse given by F(x,y)=Ax^2+Bxy+Cy^2+Dx+Ey+F = 0
		 */
		void getDrawingPoints(const double& A, const double& B,
				const double& C, const double& D, const double& E,
				const double& F, std::vector<ControlPoint>& pts) {

			/* find a point on the ellipse xs is the x at the center*/
			/* xc = (2CD-BE)/(B^2-4AC)     */
			double xs = (int) ((2 * C * D - B * E) / (B * B - 4 * A * C));
			double a = C;
			double b = B * xs + E;
			double c = A * xs * xs + D * xs + F;
			double ys = ((-b + sqrt(b * b - 4 * a * c)) / (2 * a));

			plotter(A, B, C, D, E, F, xs, ys);
			pts = plotter.getPlottingObject().getPoints();
		}

};

