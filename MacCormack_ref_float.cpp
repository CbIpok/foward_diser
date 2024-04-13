#include <cmath>
#include <iostream>

#define GA 9.81
#define H_MIN 0.0

// Quite good
// #define SMOOTH      0.25
// #define MU          0.8

// Unstable
// #define SMOOTH      0.5
// #define MU          0.8

#define N_SMOOTH 1
#define SMOOTH 0.5
#define MU 0.8

static float* eh = 0;
static float* vh = 0;
static float* uh = 0;

static float* ef = 0;
static float* vf = 0;
static float* uf = 0;

void mac_cormack_float(float* u, float* v, float* e, float* D, float* dx, float* dy, float* cosv,
    float* tg, float& dT, float grnd, float* nu, float* nv, float* ne, int nx, int ny)
{
    if ((eh == 0) || (vh == 0) || (uh == 0)) {
        eh = new float[nx * ny];
        vh = new float[nx * ny];
        uh = new float[nx * ny];

        ef = new float[nx * ny];
        vf = new float[nx * ny];
        uf = new float[nx * ny];
    }

    static int step = 0;

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            int imj = j * nx + i - 1;
            int ijm = (j - 1) * nx + i;

            float dxi = dx[i - 1] * cosv[j];
            float dyi = dy[j - 1];

            float cG = dT * sqrt(GA * D[ij]);

            float eij = e[ij], eimj = e[imj], eijm = e[ijm];
            float uij = u[ij], uimj = u[imj], uijm = u[ijm];
            float vij = v[ij], vimj = v[imj], vijm = v[ijm];
            float Dij = D[ij], Dimj = D[imj], Dijm = D[ijm];

            if (j == 0 || i == 0) {
                eh[ij] = eij;
                uh[ij] = uij;
                vh[ij] = vij;
                continue;
            }

            if (D[ij] < grnd) {
                eh[ij] = 0.0f;
                uh[ij] = 0.0f;
                vh[ij] = 0.0f;
                continue;
            }

            if ((Dimj < grnd) && (Dijm < grnd)) {
                eimj = eij;
                eijm = eij;
                uimj = 0.0f;
                uijm = 0.0f;
                vimj = 0.0f;
                vijm = 0.0f;
                Dimj = Dij;
                Dijm = Dij;
            } else if (Dimj < grnd) {
                eimj = eij;
                uimj = 0.0f;
                vimj = vij;
                Dimj = Dij;
            } else if (Dijm < grnd) {
                eijm = eij;
                uijm = uij;
                vijm = 0.0f;
                Dijm = Dij;
            }

            if (j == ny - 1) {
                eh[ij] = eij - cG * dyi * (eij - eijm);
                uh[ij] = uij - cG * dyi * (uij - uijm);
                vh[ij] = vij - cG * dyi * (vij - vijm);
            } else if (i == nx - 1) {
                eh[ij] = eij - cG * dxi * (eij - eimj);
                uh[ij] = uij - cG * dxi * (uij - uimj);
                vh[ij] = vij - cG * dxi * (vij - vimj);
            } else {
                float Hij = eij + Dij, Himj = eimj + Dimj, Hijm = eijm + Dijm;

                eh[ij] = eij
                    - dT * (dxi * (Hij * uij - Himj * uimj) + dyi * (Hij * vij - Hijm * vijm));

                uh[ij] = uij
                    - dT
                        * (0.5f * dxi * (uij * uij - uimj * uimj) + dyi * vij * (uij - uijm)
                            + dxi * GA * (eij - eimj) - 0 /* Coriolis force here */);
                vh[ij] = vij
                    - dT
                        * (dxi * uij * (vij - vimj) + 0.5f * dyi * (vij * vij - vijm * vijm)
                            + dyi * GA * (eij - eijm) - 0 /* Coriolis force here */);
            }

            if (Dij + eh[ij] < H_MIN) {
                eh[ij] = -Dij + 0.01f;
                uh[ij] = 0.0f;
                vh[ij] = 0.0f;
            }
        }
    }

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            int ipj = j * nx + i + 1;
            int ijp = (j + 1) * nx + i;

            float dxi = dx[i] * cosv[j];
            float dyi = dy[j];

            float cG = dT * sqrt(GA * D[ij]);

            float ehij = eh[ij], ehipj = eh[ipj], ehijp = eh[ijp];
            float uhij = uh[ij], uhipj = uh[ipj], uhijp = uh[ijp];
            float vhij = vh[ij], vhipj = vh[ipj], vhijp = vh[ijp];
            float Dij = D[ij], Dipj = D[ipj], Dijp = D[ijp];

            if ((j == ny - 1) || (i == nx - 1)) {
                ef[ij] = ehij;
                uf[ij] = uhij;
                vf[ij] = vhij;
                continue;
            }

            if (D[ij] < grnd) {
                ef[ij] = 0.0f;
                uf[ij] = 0.0f;
                vf[ij] = 0.0f;
                continue;
            }

            if ((Dipj < grnd) && (Dijp < grnd)) {
                ehipj = ehij;
                ehijp = ehij;
                uhipj = 0.0f;
                uhijp = 0.0f;
                vhipj = 0.0f;
                vhijp = 0.0f;
                Dipj = D[ij];
                Dijp = D[ij];
            } else if (Dipj < grnd) {
                ehipj = ehij;
                uhipj = 0.0f;
                vhipj = vhij;
                Dipj = D[ij];
            } else if (Dijp < grnd) {
                ehijp = ehij;
                uhijp = uhij;
                vhijp = 0.0f;
                Dijp = D[ij];
            }

            else if (j == 0) {
                ef[ij] = ehij + cG * dyi * (ehijp - ehij);
                uf[ij] = uhij + cG * dyi * (uhijp - uhij);
                vf[ij] = vhij + cG * dyi * (vhijp - vhij);
            } else if (i == 0) {
                ef[ij] = ehij + cG * dxi * (ehipj - ehij);
                uf[ij] = uhij + cG * dxi * (uhipj - uhij);
                vf[ij] = vhij + cG * dxi * (vhipj - vhij);
            } else {
                float hhij = ehij + Dij, hhipj = ehipj + Dipj, hhijp = ehijp + Dijp;

                ef[ij] = 0.5f * (e[ij] + ehij)
                    - 0.5f * dT
                        * (dxi * (hhipj * uhipj - hhij * uhij)
                            + dyi * (hhijp * vhijp - hhij * vhij));

                uf[ij] = 0.5f * (u[ij] + uhij)
                    - 0.5f * dT
                        * (0.5f * dxi * (uhipj * uhipj - uhij * uhij) + dyi * vhij * (uhijp - uhij)
                            + dxi * GA * (ehipj - ehij) - 0 /* Coriolis force here */);
                vf[ij] = 0.5f * (v[ij] + vhij)
                    - 0.5f * dT
                        * (dxi * uhij * (vhipj - vhij) + 0.5f * dyi * (vhijp * vhijp - vhij * vhij)
                            + dyi * GA * (ehijp - ehij) - 0 /* Coriolis force here */);
            }

            if (Dij + ef[ij] < H_MIN) {
                ef[ij] = -Dij + 0.01f;
                uf[ij] = 0.0f;
                vf[ij] = 0.0f;
            }
        }
    }

    step++;

    float dxi0 = dx[0] * cosv[0];
    float dxin = dx[nx - 2] * cosv[nx - 2];

    float dyi0 = dy[0];
    float dyin = dy[ny - 2];

    float wx00 = dxi0 / (dxi0 + dyi0);
    float wy00 = dyi0 / (dxi0 + dyi0);
    float wxn0 = dxin / (dxin + dyi0);
    float wyn0 = dyi0 / (dxin + dyi0);
    float wx0n = dxi0 / (dxi0 + dyin);
    float wy0n = dyin / (dxi0 + dyin);
    float wxnn = dxin / (dxin + dyin);
    float wynn = dyin / (dxin + dyin);

    ef[0] = wx00 * e[nx] + wy00 * e[1];
    uf[0] = wx00 * u[nx] + wy00 * u[1];
    vf[0] = wx00 * v[nx] + wy00 * v[1];

    ef[nx - 1] = wxn0 * e[nx + nx - 1] + wyn0 * e[nx - 2];
    uf[nx - 1] = wxn0 * u[nx + nx - 1] + wyn0 * u[nx - 2];
    vf[nx - 1] = wxn0 * v[nx + nx - 1] + wyn0 * v[nx - 2];

    ef[(ny - 1) * nx] = wx0n * e[(ny - 2) * nx] + wy0n * e[(ny - 1) * nx + 1];
    uf[(ny - 1) * nx] = wx0n * u[(ny - 2) * nx] + wy0n * u[(ny - 1) * nx + 1];
    vf[(ny - 1) * nx] = wx0n * v[(ny - 2) * nx] + wy0n * v[(ny - 1) * nx + 1];

    ef[(ny - 1) * nx + nx - 1]
        = wynn * e[(ny - 1) * nx + nx - 2] + wxnn * e[(ny - 2) * nx + nx - 1];
    uf[(ny - 1) * nx + nx - 1]
        = wynn * u[(ny - 1) * nx + nx - 2] + wxnn * u[(ny - 2) * nx + nx - 1];
    vf[(ny - 1) * nx + nx - 1]
        = wynn * v[(ny - 1) * nx + nx - 2] + wxnn * v[(ny - 2) * nx + nx - 1];

    int sc = 0;

    /*for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            int ij = j * nx + i;
            ne[ij] = ef[ij];
            nu[ij] = uf[ij];
            nv[ij] = vf[ij];
        }
    }*/

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            int ipj = j * nx + i + 1;
            int ijp = (j + 1) * nx + i;
            int imj = j * nx + i - 1;
            int ijm = (j - 1) * nx + i;

            if ((i == 0) || (j == 0) || (i == nx - 1) || (j == ny - 1)) {
                ne[ij] = ef[ij];
                nu[ij] = uf[ij];
                nv[ij] = vf[ij];
                continue;
            }

            float Hij = ef[ij] + D[ij];

            if ((D[ij] < grnd) || (Hij <= 0.0) || (fabs(ef[ij] / Hij) <= SMOOTH)
                || ((step % N_SMOOTH) != 0)) {
                ne[ij] = ef[ij];
                nu[ij] = uf[ij];
                nv[ij] = vf[ij];
                continue;
            }
            sc++;

            float Himj = ef[imj] + D[imj];
            float Hijm = ef[ijm] + D[ijm];
            float Hipj = ef[ipj] + D[ipj];
            float Hijp = ef[ijp] + D[ijp];

            float Ae = 0.0f, Au = 0.0f, Av = 0.0f;
            int C = 0;

            if ((i > 0) && (D[imj] >= grnd) && (Himj > 0.0)) {
                Ae += ef[imj];
                Au += uf[imj];
                Av += vf[imj];
                C++;
            }
            if ((i < nx - 1) && (D[ipj] >= grnd) && (Hipj > 0.0)) {
                Ae += ef[ipj];
                Au += uf[ipj];
                Av += vf[ipj];
                C++;
            }
            if ((j > 0) && (D[ijm] >= grnd) && (Hijm > 0.0)) {
                Ae += ef[ijm];
                Au += uf[ijm];
                Av += vf[ijm];
                C++;
            }
            if ((j < ny - 1) && (D[ijp] >= grnd) && (Hijp > 0.0)) {
                Ae += ef[ijp];
                Au += uf[ijp];
                Av += vf[ijp];
                C++;
            }

            if (!C) {
                ne[ij] = ef[ij];
                nu[ij] = uf[ij];
                nv[ij] = vf[ij];
                continue;
            }

            ne[ij] = (1.0 - MU) * ef[ij] + MU * Ae / C;
            nu[ij] = (1.0 - MU) * uf[ij] + MU * Au / C;
            nv[ij] = (1.0 - MU) * vf[ij] + MU * Av / C;

            if (D[ij] + ne[ij] < H_MIN) {
                ne[ij] = -D[ij] + 0.01f;
                nu[ij] = 0.0f;
                nv[ij] = 0.0f;
            }
        }
    }
    std::cout << "Smoothed " << sc << " points" << std::endl;

    /*#pragma omp parallel for
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            int ij = j * nx + i;
            int ipj = j * nx + i + 1;
            int ijp = (j + 1) * nx + i;
            int imj = j * nx + i - 1;
            int ijm = (j - 1) * nx + i;

            if ((i == 0) || (j == 0) || (i == nx - 1) || (j == ny - 1))
            {
                ne[ij] = ef[ij];
                nu[ij] = uf[ij];
                nv[ij] = vf[ij];
                continue;
            }

            float Hij = ef[ij] + D[ij];

            float c = M_RUNUP * sqrt(GA * (ef[ij] + D[ij]));

            if ((D[ij] < grnd) || (D[ij] > D_DEEP) || ((uf[ij] <= c) && (vf[ij] <= c)))
            {
                ne[ij] = ef[ij];
                nu[ij] = uf[ij];
                nv[ij] = vf[ij];
                continue;
            }

            int ijmin = ij;

            if ((fabs(uf[ipj]) < fabs(uf[ijmin])) && (D[ipj] >= grnd))
                ijmin = ipj;
            if ((fabs(uf[ijp]) < fabs(uf[ijmin])) && (D[ijp] >= grnd))
                ijmin = ijp;
            if ((fabs(uf[imj]) < fabs(uf[ijmin])) && (D[imj] >= grnd))
                ijmin = imj;
            if ((fabs(uf[ijm]) < fabs(uf[ijmin])) && (D[ijm] >= grnd))
                ijmin = ijm;

            uf[ij] = uf[ijmin];
            ef[ij] = ef[ijmin];

            if ((fabs(vf[ipj]) < fabs(vf[ijmin])) && (D[ipj] >= grnd))
                ijmin = ipj;
            if ((fabs(vf[ijp]) < fabs(vf[ijmin])) && (D[ijp] >= grnd))
                ijmin = ijp;
            if ((fabs(vf[imj]) < fabs(vf[ijmin])) && (D[imj] >= grnd))
                ijmin = imj;
            if ((fabs(vf[ijm]) < fabs(vf[ijmin])) && (D[ijm] >= grnd))
                ijmin = ijm;

            vf[ij] = vf[ijmin];
            ef[ij] = ef[ijmin];
        }
    }*/
}
