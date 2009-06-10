#include <vector>
#include <math.h>
#include <iostream>

template<typename T>
class Vector2D {
		T** vec;
		int l;
		int c;
	public:
		Vector2D(int lin, int col) {
			vec = new T*[lin];
			for (int i = 0; i < lin; i++) {
				vec[i] = new T[col];
			}
			this->l = lin;
			this->c = col;
		}
		virtual ~Vector2D() {
			for (int i = 0; i < l; i++) {
				delete[] vec[i];
			}
			delete[] vec;
		}

		T& operator()(int lin, int col) {
			return vec[lin][col];
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

class CurvePanel {
	private:
		int i;

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

			double b[n + 1];
			double z[n + 1];

			for (ip = 1; ip <= n; ip++) {
				for (iq = 1; iq <= n; iq++)
					v(ip, iq) = 0.0;
				v(ip, ip) = 1.0;
			}
			for (ip = 1; ip <= n; ip++) {
				b[ip] = d[ip] = a(ip, ip);
				z[ip] = 0.0;
			}
			nrot = 0;
			for (i = 1; i <= 50; i++) {
				sm = 0.0;
				for (ip = 1; ip <= n - 1; ip++) {
					for (iq = ip + 1; iq <= n; iq++)
						sm += fabs(a(ip, iq));
				}
				if (sm == 0.0) {
					/*
					 * free_vector(z,1,n); free_vector(b,1,n);
					 */
					return;
				}
				if (i < 4)
					tresh = 0.2 * sm / (n * n);
				else
					tresh = 0.0;
				for (ip = 1; ip <= n - 1; ip++) {
					for (iq = ip + 1; iq <= n; iq++) {
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
							for (j = 1; j <= ip - 1; j++) {
								ROTATE(a, j, ip, j, iq, tau, s);
							}
							for (j = ip + 1; j <= iq - 1; j++) {
								ROTATE(a, ip, j, j, iq, tau, s);
							}
							for (j = iq + 1; j <= n; j++) {
								ROTATE(a, ip, j, iq, j, tau, s);
							}
							for (j = 1; j <= n; j++) {
								ROTATE(v, j, ip, j, iq, tau, s);
							}
							++nrot;
						}
					}
				}
				for (ip = 1; ip <= n; ip++) {
					b[ip] += z[ip];
					d[ip] = b[ip];
					z[ip] = 0.0;
				}
			}
			// printf("Too many iterations in routine JACOBI");
		}

		// Perform the Cholesky decomposition
		// Return the lower triangular L such that L*L'=A
		void choldc(Vector2D<double>& a, int n, Vector2D<double>& l) {
			int i, j, k;
			double sum;
			double p[n + 1];

			for (i = 1; i <= n; i++) {
				for (j = i; j <= n; j++) {
					for (sum = a(i, j), k = i - 1; k >= 1; k--)
						sum -= a(i, k) * a(j, k);
					if (i == j) {
						if (sum <= 0.0)
						// printf("\nA is not poitive definite!");
						{
						} else
							p[i] = sqrt(sum);
					} else {
						a(j, i) = sum / p[i];
					}
				}
			}
			for (i = 1; i <= n; i++)
				for (j = i; j <= n; j++)
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
			Vector2D<double> B(N + 1, N + 2);
			Vector2D<double> A(N + 1, 2 * N + 2);
			Vector2D<double> C(N + 1, N + 1);
			double eps = 10e-20;

			for (k = 1; k <= N; k++)
				for (j = 1; j <= N; j++)
					B(k, j) = TB(k, j);

			for (k = 1; k <= N; k++) {
				for (j = 1; j <= N + 1; j++)
					A(k, j) = B(k, j);
				for (j = N + 2; j <= 2 * N + 1; j++)
					A(k, j) = (float) 0;
				A(k, k - 1 + N + 2) = (float) 1;
			}
			for (k = 1; k <= N; k++) {
				maxpivot = fabs((double) A(k, k));
				npivot = k;
				for (i = k; i <= N; i++)
					if (maxpivot < fabs((double) A(i, k))) {
						maxpivot = fabs((double) A(i, k));
						npivot = i;
					}
				if (maxpivot >= eps) {
					if (npivot != k)
						for (j = k; j <= 2 * N + 1; j++) {
							temp = A(npivot, j);
							A(npivot, j) = A(k, j);
							A(k, j) = temp;
						};
					D = A(k, k);
					for (j = 2 * N + 1; j >= k; j--)
						A(k, j) = A(k, j) / D;
					for (i = 1; i <= N; i++) {
						if (i != k) {
							mult = A(i, k);
							for (j = 2 * N + 1; j >= k; j--)
								A(i, j) = A(i, j) - mult * A(k, j);
						}
					}
				} else { // printf("\n The matrix may be singular !!") ;
					return (-1);
				};
			}
			/** Copia il risultato nella matrice InvB ***/
			for (k = 1, p = 1; k <= N; k++, p++)
				for (j = N + 2, q = 1; j <= 2 * N + 1; j++, q++)
					InvB(p, q) = A(k, j);
			return (0);
		} /* End of INVERSE */

		void AperB(Vector2D<double>& _A, Vector2D<double>& _B,
				Vector2D<double>& _res, int _righA, int _colA, int _righB,
				int _colB) {
			int p, q, l;
			for (p = 1; p <= _righA; p++)
				for (q = 1; q <= _colB; q++) {
					_res(p, q) = 0.0;
					for (l = 1; l <= _colA; l++)
						_res(p, q) = _res(p, q) + _A(p, l) * _B(l, q);
				}
		}

		void A_TperB(Vector2D<double>& _A, Vector2D<double>& _B, Vector2D<
				double>& _res, int _righA, int _colA, int _righB, int _colB) {
			int p, q, l;
			for (p = 1; p <= _colA; p++)
				for (q = 1; q <= _colB; q++) {
					_res(p, q) = 0.0;
					for (l = 1; l <= _righA; l++)
						_res(p, q) = _res(p, q) + _A(l, p) * _B(l, q);
				}
		}

		void AperB_T(Vector2D<double>& _A, Vector2D<double>& _B, Vector2D<
				double>& _res, int _righA, int _colA, int _righB, int _colB) {
			int p, q, l;
			for (p = 1; p <= _colA; p++)
				for (q = 1; q <= _colB; q++) {
					_res(p, q) = 0.0;
					for (l = 1; l <= _righA; l++)
						_res(p, q) = _res(p, q) + _A(p, l) * _B(q, l);
				}
		}
	public:
		void operator()(double& xc, double& yc, double& el_a, double& el_b) {
			int np = points.size(); // number of points
			Vector2D<double> D(np + 1, 6);
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
			int npts = 50;

			Vector2D<double> XY(3, npts + 1);
			double pvec[6];

			Const(1, 2) = -2;
			Const(2, 1) = -2;

			if (np < 6)
				return;

			// Now first fill design matrix
			for (int i = 1; i <= np; i++) {
				tx = points[i - 1].x;
				ty = points[i - 1].y;

				D(i, 1) = tx * tx;
				D(i, 2) = ty * ty;
				D(i, 3) = tx;
				D(i, 4) = ty;
				D(i, 5) = 1.0;
			}

			A_TperB(D, D, S, np, 5, np, 5);

			choldc(S, 5, L);

			inverse(L, invL, 5);

			AperB_T(Const, invL, temp, 5, 5, 5, 5);
			AperB(invL, temp, C, 5, 5, 5, 5);

			jacobi(C, 5, d, V, nrot);

			A_TperB(invL, V, sol, 5, 5, 5, 5);

			for (int j = 1; j <= 5; j++) /* Scan columns */
			{
				double mod = 0.0;
				for (int i = 1; i <= 5; i++)
					mod += sol(i, j) * sol(i, j);
				for (int i = 1; i <= 5; i++)
					sol(i, j) /= sqrt(mod);
			}

			double zero = 10e-20;
			int solind = 0;

			for (int i = 1; i <= 5; i++)
				if (d[i] < 0 && fabs(d[i]) > zero)
					solind = i;

			// Now fetch the right solution
			for (int j = 1; j <= 5; j++)
				pvec[j] = sol(j, solind);

			double delta = pvec[3] * pvec[3] + pvec[4] * pvec[4] * pvec[1]
					/ pvec[2] - 4 * pvec[1] * pvec[5];
			xc = (-pvec[3] / (2 * pvec[1]));
			yc = (-pvec[4] / (2 * pvec[2]));
			el_a = (sqrt(delta) / (2 * fabs(pvec[1])));
			el_b = (sqrt(delta) / (2 * sqrt(pvec[1] * pvec[2])));

		}
};

