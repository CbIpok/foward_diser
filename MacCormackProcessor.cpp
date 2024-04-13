#include "MacCormackProcessor.hpp"
#include <cmath>

#define GA 9.81

namespace tsunami {
template <class F> void MacCormackProcessor<F>::prepare(DataArrayPtr<F> bathymetry)
{
    TsunamiProcessor<F>::prepare(bathymetry);

    nx = bathymetry->width();
    ny = bathymetry->height();

    D = new F[nx * ny];
    u = new F[nx * ny];
    v = new F[nx * ny];
    e = new F[nx * ny];
    uh = new F[nx * ny];
    vh = new F[nx * ny];
    eh = new F[nx * ny];
    uf = new F[nx * ny];
    vf = new F[nx * ny];
    ef = new F[nx * ny];
    cosv = new F[ny];
    dx = new F[nx];
    dy = new F[ny];

    for (int j = 0; j < ny; j++) {
        bathymetry->data(0, j) = bathymetry->data(1, j);
        bathymetry->data(nx - 1, j) = bathymetry->data(nx - 2, j);
    }

    for (int j = 0; j < nx; j++) {
        bathymetry->data(j, 0) = bathymetry->data(j, 1);
        bathymetry->data(j, ny - 1) = bathymetry->data(j, ny - 2);
    }

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            D[ij] = bathymetry->data(i, j);
            u[ij] = 0.0;
            v[ij] = 0.0;
            e[ij] = 0.0;
        }
    }

    if (bathymetry->type() == DataArray<F>::Spheric) {
        for (int j = 0; j < nx - 1; j++)
            dx[j]
                = 1.0 / (fabs((F)bathymetry->xCoord(j + 1) - (F)bathymetry->xCoord(j)) * 111320.0);

        for (int j = 0; j < ny - 1; j++)
            dy[j]
                = 1.0 / (fabs((F)bathymetry->yCoord(j + 1) - (F)bathymetry->yCoord(j)) * 111320.0);

        for (int j = 0; j < ny; j++)
            cosv[j] = 1.0 / cos(M_PI / 180.0 * (F)bathymetry->yCoord(j));
    } else {
        for (int j = 0; j < nx - 1; j++)
            dx[j] = 1.0 / bathymetry->xStep();

        for (int j = 0; j < ny - 1; j++)
            dy[j] = 1.0 / bathymetry->yStep();

        for (int j = 0; j < ny; j++)
            cosv[j] = 1.0;
    }
    m_outHeight = DataArray<F>::create(
        nx, ny, bathymetry->xStep(), bathymetry->yStep(), bathymetry->type());
    m_outDataValid = false;
}

template <class F> void MacCormackProcessor<F>::addDeformation(DataArrayPtr<F> deformation)
{
    TsunamiProcessor<F>::addDeformation(deformation);
}

template <class F>
void MacCormackProcessor<F>::init(DataArrayPtr<F> h, DataArrayPtr<F> vx, DataArrayPtr<F> vy)
{
    TsunamiProcessor<F>::init(h, vx, vy);

    memset(u, 0, nx * ny * sizeof(F));
    memset(v, 0, nx * ny * sizeof(F));
    memset(e, 0, nx * ny * sizeof(F));

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            if (vx)
                u[ij] = vx->data(i, j);
            if (vy)
                v[ij] = vy->data(i, j);
            if (h)
                e[ij] = h->data(i, j);
        }
    }
}

template <class F> void MacCormackProcessor<F>::parameters(const json& params)
{
    grnd = params["ground"].get<F>();
}

template <class F> void MacCormackProcessor<F>::step(F dT, int& number)
{
    m_outDataValid = false;
    F time = number * dT;

    for (auto proc : MacCormackProcessor<F>::m_preProcessors) {
        (*proc)(nx, ny, number, dT, e, u, v, D, ef, uf, vf);
        if (!proc->inPlace()) {
            std::swap(e, ef);
            std::swap(u, uf);
            std::swap(v, vf);
        }
    }

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            int imj = j * nx + i - 1;
            int ijm = (j - 1) * nx + i;

            F dxi = dx[i - 1] * cosv[j];
            F dyi = dy[j - 1];

            F cG = dT * sqrt(GA * D[ij]);

            F eij = e[ij], eimj = e[imj], eijm = e[ijm];
            F uij = u[ij], uimj = u[imj], uijm = u[ijm];
            F vij = v[ij], vimj = v[imj], vijm = v[ijm];
            F Dij = D[ij], Dimj = D[imj], Dijm = D[ijm];

            if (j == 0 || i == 0) {
                eh[ij] = eij;
                uh[ij] = uij;
                vh[ij] = vij;
                continue;
            }

            if (D[ij] < grnd) {
                eh[ij] = 0.0;
                uh[ij] = 0.0;
                vh[ij] = 0.0;
                continue;
            }

            if ((Dimj < grnd) && (Dijm < grnd)) {
                eimj = eij;
                eijm = eij;
                uimj = 0.0;
                uijm = 0.0;
                vimj = 0.0;
                vijm = 0.0;
                Dimj = Dij;
                Dijm = Dij;
            } else if (Dimj < grnd) {
                eimj = eij;
                uimj = 0.0;
                vimj = vij;
                Dimj = Dij;
            } else if (Dijm < grnd) {
                eijm = eij;
                uijm = uij;
                vijm = 0.0;
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
                F Hij = eij + Dij, Himj = eimj + Dimj, Hijm = eijm + Dijm;

                eh[ij] = eij
                    - dT * (dxi * (Hij * uij - Himj * uimj) + dyi * (Hij * vij - Hijm * vijm));

                uh[ij] = uij
                    - dT
                        * (0.5 * dxi * (uij * uij - uimj * uimj) + dyi * vij * (uij - uijm)
                            + dxi * GA * (eij - eimj) - 0 /* Coriolis force here */);
                vh[ij] = vij
                    - dT
                        * (dxi * uij * (vij - vimj) + 0.5 * dyi * (vij * vij - vijm * vijm)
                            + dyi * GA * (eij - eijm) - 0 /* Coriolis force here */);
            }

            if (Dij + eh[ij] < 0.01) {
                eh[ij] = -Dij;
                uh[ij] = 0.0;
                vh[ij] = 0.0;
            }
        }
    }

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            int ipj = j * nx + i + 1;
            int ijp = (j + 1) * nx + i;

            F dxi = dx[i] * cosv[j];
            F dyi = dy[j];

            F cG = dT * sqrt(GA * D[ij]);

            F ehij = eh[ij], ehipj = eh[ipj], ehijp = eh[ijp];
            F uhij = uh[ij], uhipj = uh[ipj], uhijp = uh[ijp];
            F vhij = vh[ij], vhipj = vh[ipj], vhijp = vh[ijp];
            F Dij = D[ij], Dipj = D[ipj], Dijp = D[ijp];

            if ((j == ny - 1) || (i == nx - 1)) {
                ef[ij] = ehij;
                uf[ij] = uhij;
                vf[ij] = vhij;
                continue;
            }

            if (D[ij] < grnd) {
                ef[ij] = 0.0;
                uf[ij] = 0.0;
                vf[ij] = 0.0;
                continue;
            }

            if ((Dipj < grnd) && (Dijp < grnd)) {
                ehipj = ehij;
                ehijp = ehij;
                uhipj = 0.0;
                uhijp = 0.0;
                vhipj = 0.0;
                vhijp = 0.0;
                Dipj = D[ij];
                Dijp = D[ij];
            } else if (Dipj < grnd) {
                ehipj = ehij;
                uhipj = 0.0;
                vhipj = vhij;
                Dipj = D[ij];
            } else if (Dijp < grnd) {
                ehijp = ehij;
                uhijp = uhij;
                vhijp = 0.0;
                Dijp = D[ij];
            }

            if (j == 0) {
                ef[ij] = ehij + cG * dyi * (ehijp - ehij);
                uf[ij] = uhij + cG * dyi * (uhijp - uhij);
                vf[ij] = vhij + cG * dyi * (vhijp - vhij);
            } else if (i == 0) {
                ef[ij] = ehij + cG * dxi * (ehipj - ehij);
                uf[ij] = uhij + cG * dxi * (uhipj - uhij);
                vf[ij] = vhij + cG * dxi * (vhipj - vhij);
            } else {
                F hhij = ehij + Dij, hhipj = ehipj + Dipj, hhijp = ehijp + Dijp;

                ef[ij] = 0.5 * (e[ij] + ehij)
                    - 0.5 * dT
                        * (dxi * (hhipj * uhipj - hhij * uhij)
                            + dyi * (hhijp * vhijp - hhij * vhij));

                uf[ij] = 0.5 * (u[ij] + uhij)
                    - 0.5 * dT
                        * (0.5 * dxi * (uhipj * uhipj - uhij * uhij) + dyi * vhij * (uhijp - uhij)
                            + dxi * GA * (ehipj - ehij) - 0 /* Coriolis force here */);
                vf[ij] = 0.5 * (v[ij] + vhij)
                    - 0.5 * dT
                        * (dxi * uhij * (vhipj - vhij) + 0.5 * dyi * (vhijp * vhijp - vhij * vhij)
                            + dyi * GA * (ehijp - ehij) - 0 /* Coriolis force here */);
            }

            if (Dij + ef[ij] < 0.001) {
                ef[ij] = -Dij;
                uf[ij] = 0.0f;
                vf[ij] = 0.0f;
            }
        }
    }

    F dxi0 = dx[0] * cosv[0];
    F dxin = dx[nx - 2] * cosv[nx - 2];

    F dyi0 = dy[0];
    F dyin = dy[ny - 2];

    F wx00 = dxi0 / (dxi0 + dyi0);
    F wy00 = dyi0 / (dxi0 + dyi0);
    F wxn0 = dxin / (dxin + dyi0);
    F wyn0 = dyi0 / (dxin + dyi0);
    F wx0n = dxi0 / (dxi0 + dyin);
    F wy0n = dyin / (dxi0 + dyin);
    F wxnn = dxin / (dxin + dyin);
    F wynn = dyin / (dxin + dyin);

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

    std::swap(e, ef);
    std::swap(u, uf);
    std::swap(v, vf);

    for (auto proc : MacCormackProcessor<F>::m_postProcessors) {
        (*proc)(nx, ny, number, dT, e, u, v, D, ef, uf, vf);
        if (!proc->inPlace()) {
            std::swap(e, ef);
            std::swap(u, uf);
            std::swap(v, vf);
        }
    }
}

template <class F> void MacCormackProcessor<F>::finalize() { TsunamiProcessor<F>::finalize(); }

template <class F> static void outputData(F* data, DataArrayPtr<F> out)
{
    memcpy(out->data(), data, out->width() * out->height() * sizeof(F));
}

static void outputData(float* data, DataArrayPtr<double> out)
{
    for (int i = 0; i < out->height(); i++)
        for (int j = 0; j < out->width(); j++)
            out->data(j, i) = data[i * out->width() + j];
}

static void outputData(double* data, DataArrayPtr<float> out)
{
    for (int i = 0; i < out->height(); i++)
        for (int j = 0; j < out->width(); j++)
            out->data(j, i) = data[i * out->width() + j];
}

template <class F> DataArrayPtr<F> MacCormackProcessor<F>::height() const
{
    if (!m_outDataValid) {
        outputData(e, m_outHeight);
        m_outDataValid = true;
    }
    return m_outHeight;
}

template class MacCormackProcessor<float>;
template class MacCormackProcessor<double>;
}
