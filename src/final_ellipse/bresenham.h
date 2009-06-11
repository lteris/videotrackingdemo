/*
 * bresenham.h
 *
 *  Created on: Jun 11, 2009
 *      Author: liviu
 */

#ifndef BRESENHAM_H_
#define BRESENHAM_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/** Foley and Van Dam Computer Graphics: Principles and Practice in C p 959 */

template<typename PLOT>
class ConicPlotter {
	private:
		/* A xx + B xy + C yy + D x + E y + F = 0 */
		double A;
		double B;
		double C;
		double D;
		double E;
		double F;

		PLOT plot;

		int odd(int n) {
			return n & 1;
		}

		int getoctant(int gx, int gy) {
			/* Use gradient to identify octant */
			int upper = abs(gx) > abs(gy);
			if (gx >= 0) // Right-pointing
				if (gy >= 0) //    Up
					return 4 - upper;
				else
					//    Down
					return 1 + upper;
			else // Left
			if (gy > 0) //    Up
				return 5 + upper;
			else
				//    Down
				return 8 - upper;
		}

		/** routine that detects the points - calls plot(x,y) for each detected point */
		void draw(int xs, int ys, int xe, int ye) {
			static int DIAGx[] = { 999, 1, 1, -1, -1, -1, -1, 1, 1 };
			static int DIAGy[] = { 999, 1, 1, 1, 1, -1, -1, -1, -1 };
			static int SIDEx[] = { 999, 1, 0, 0, -1, -1, 0, 0, 1 };
			static int SIDEy[] = { 999, 0, 1, 1, 0, 0, -1, -1, 0 };

			A *= 4;
			B *= 4;
			C *= 4;
			D *= 4;
			E *= 4;
			F *= 4;

			/* Translate start point to origin... */
			F = A * xs * xs + B * xs * ys + C * ys * ys + D * xs + E * ys + F;
			D = D + 2 * A * xs + B * ys;
			E = E + B * xs + 2 * C * ys;

			/* Work out starting octant */
			int octant = getoctant(D, E);

			int dxS = SIDEx[octant];
			int dyS = SIDEy[octant];
			int dxD = DIAGx[octant];
			int dyD = DIAGy[octant];

			double d, u, v;
			switch (octant) {
			case 1:
				d = A + B / 2 + C / 4 + D + E / 2 + F;
				u = A + B / 2 + D;
				v = u + E;
				break;
			case 2:
				d = A / 4 + B / 2 + C + D / 2 + E + F;
				u = B / 2 + C + E;
				v = u + D;
				break;
			case 3:
				d = A / 4 - B / 2 + C - D / 2 + E + F;
				u = -B / 2 + C + E;
				v = u - D;
				break;
			case 4:
				d = A - B / 2 + C / 4 - D + E / 2 + F;
				u = A - B / 2 - D;
				v = u + E;
				break;
			case 5:
				d = A + B / 2 + C / 4 - D - E / 2 + F;
				u = A + B / 2 - D;
				v = u - E;
				break;
			case 6:
				d = A / 4 + B / 2 + C - D / 2 - E + F;
				u = B / 2 + C - E;
				v = u - D;
				break;
			case 7:
				d = A / 4 - B / 2 + C + D / 2 - E + F;
				u = -B / 2 + C - E;
				v = u + D;
				break;
			case 8:
				d = A - B / 2 + C / 4 + D - E / 2 + F;
				u = A - B / 2 + D;
				v = u - E;
				break;
			default:
				return;
			}

			double k1sign = dyS * dyD;
			double k1 = 2 * (A + k1sign * (C - A));
			double Bsign = dxD * dyD;
			double k2 = k1 + Bsign * B;
			double k3 = 2 * (A + C + Bsign * B);

			int octantcount = 8;

			int x = xs;
			int y = ys;

			while (octantcount > 0) {
				if (odd(octant)) {
					while (2 * v <= k2) {
						plot(x, y);
						if (d < 0) { /* Inside the contour */
							x = x + dxS;
							y = y + dyS;
							u = u + k1;
							v = v + k2;
							d = d + u;
						} else { /* outside the contour */
							x = x + dxD;
							y = y + dyD;
							u = u + k2;
							v = v + k3;
							d = d + v;
						}
					}

					d = d - u + v / 2 - k2 / 2 + 3 * k3 / 8;
					/* error (^) in Foley and van Dam p 959, "2nd ed, revised 5th printing" */
					u = -u + v - k2 / 2 + k3 / 2;
					v = v - k2 + k3 / 2;
					k1 = k1 - 2 * k2 + k3;
					k2 = k3 - k2;
					int tmp = dxS;
					dxS = -dyS;
					dyS = tmp;
				} else { /* Octant is even */
					while (2 * u < k2) {
						plot(x, y);
						if (d > 0) { /* Outside */
							x = x + dxS;
							y = y + dyS;
							u = u + k1;
							v = v + k2;
							d = d + u;
						} else { /* Inside */
							x = x + dxD;
							y = y + dyD;
							u = u + k2;
							v = v + k3;
							d = d + v;
						}
					}
					int tmpdk = k1 - k2;
					d = d + u - v + tmpdk;
					v = 2 * u - v + tmpdk;
					u = u + tmpdk;
					k3 = k3 + 4 * tmpdk;
					k2 = k1 + tmpdk;

					int tmp = dxD;
					dxD = -dyD;
					dyD = tmp;
				}

				octant = (octant & 7) + 1;
				octantcount--;
			}

			/* Draw final octant until we reach the endpoint */
			if (odd(octant)) {
				while (2 * v <= k2) {
					plot(x, y);
					if (x == xe && y == ye)
						break;
					if (d < 0) { /* Inside */
						x = x + dxS;
						y = y + dyS;
						u = u + k1;
						v = v + k2;
						d = d + u;
					} else { /* outside */
						x = x + dxD;
						y = y + dyD;
						u = u + k2;
						v = v + k3;
						d = d + v;
					}
				}
			} else { /* Octant is even */
				while ((2 * u < k2)) {
					plot(x, y);
					if (x == xe && y == ye)
						break;
					if (d > 0) { /* Outside*/
						x = x + dxS;
						y = y + dyS;
						u = u + k1;
						v = v + k2;
						d = d + u;
					} else { /* Inside */
						x = x + dxD;
						y = y + dyD;
						u = u + k2;
						v = v + k3;
						d = d + v;
					}
				}
			}

		}
	public:
		ConicPlotter() {
			A = B = C = D = E = F = 0;
		}
		ConicPlotter(const double& A, const double& B, const double& C,
				const double& D, const double& E, const double& F) {
			this->A = A;
			this->B = B;
			this->C = C;
			this->D = D;
			this->E = E;
			this->F = F;
		}

		/** arguments - ellipse parameters and a points on the ellipse */
		void operator()(const double& A, const double& B, const double& C, const double& D,
				const double& E, const double& F, int xs, int ys) {
			this->A = A;
			this->B = B;
			this->C = C;
			this->D = D;
			this->E = E;
			this->F = F;

			draw(xs, ys, xs, ys);
		}

		PLOT& getPlottingObject() {
			return plot;
		}
};

#endif /* BRESENHAM_H_ */
