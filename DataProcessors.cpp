#include "DataProcessors.hpp"

namespace tsunami {

template <class F> std::shared_ptr<DataProcessor<F>> DataProcessor<F>::fromJSON(const json& node)
{
    if (node["type"] == "border-sin-wave") {
        F A = node["amplitide"].get<F>();
        F T = node["period"].get<F>();
        json sides = node["sides"];

        bool top = sides[0].get<bool>();
        bool bottom = sides[1].get<bool>();
        bool left = sides[2].get<bool>();
        bool right = sides[3].get<bool>();

        return std::make_shared<BorderSinWaveProcessor<F>>(A, T, top, bottom, left, right);
    } else if (node["type"] == "offshore-smooth") {
        int period = node["period"].get<F>();
        F gnd = node["ground"].get<F>();
        F bound = node["bound"].get<F>();
        F coeff = node["coeff"].get<F>();

        return std::make_shared<OffshoreSmoothWave<F>>(period, gnd, bound, coeff);
    } else if (node["type"] == "limit-wave-on-deep") {
        int period = node["period"].get<F>();
        F gnd = node["ground"].get<F>();
        F deep = node["deep"].get<F>();

        return std::make_shared<LimitWaveOnDeep<F>>(period, deep, gnd);
    } else if (node["type"] == "remove-single-islands") {
        F gnd = node["ground"].get<F>();

        return std::make_shared<RemoveSingleIslands<F>>(gnd);
    }
    return nullptr;
}

template <class F>
BorderSinWaveProcessor<F>::BorderSinWaveProcessor(
    F amplitude, F period, bool top, bool bottom, bool left, bool right)
    : DataProcessor<F>(true)
    , m_amplitude(amplitude)
    , m_period(period)
    , m_top(top)
    , m_bottom(bottom)
    , m_left(left)
    , m_right(right)
{
}

template <class F>
void BorderSinWaveProcessor<F>::operator()(
    int nx, int ny, int step, F dT, F* e, F* vx, F* vy, F* D, F* oe, F* ovx, F* ovy)
{
    auto ef = [this, step, dT]() -> F {
        F e = 0.5 * m_amplitude * (1.0 - cos(2.0 * step * dT * M_PI / m_period));
        return e;
    };
    auto uf = [](F e, F D) -> F { return e * sqrt(9.81 / D); };

    if (step * dT > m_period)
        return;

    if (m_top) {
        for (int x = 0; x < nx; x++) {
            int idx = x;
            F v = ef();
            e[idx] = v;
            vx[idx] = 0.0;
            vy[idx] = uf(v, D[idx]);
        }
    }

    if (m_bottom) {
        for (int x = 0; x < nx; x++) {
            int idx = (ny - 1) * nx + x;
            F v = ef();
            e[idx] = v;
            vx[idx] = 0.0;
            vy[idx] = -uf(v, D[idx]);
        }
    }

    if (m_left) {
        for (int y = 0; y < ny; y++) {
            int idx = y * nx;
            F v = ef();
            e[idx] = v;
            vx[idx] = uf(v, D[idx]);
            vy[idx] = 0.0;
        }
    }

    if (m_right) {
        for (int y = 0; y < ny; y++) {
            int idx = y * nx + nx - 1;
            F v = ef();
            e[idx] = v;
            vx[idx] = -uf(v, D[idx]);
            vy[idx] = 0.0;
        }
    }
}

template <class F>
OffshoreSmoothWave<F>::OffshoreSmoothWave(int period, F grnd, F bound, F coeff)
    : DataProcessor<F>(false)
    , m_period(period)
    , m_grnd(grnd)
    , m_bound(bound)
    , m_coeff(coeff)
{
}

template <class F>
void OffshoreSmoothWave<F>::operator()(
    int nx, int ny, int step, F dT, F* ef, F* uf, F* vf, F* D, F* ne, F* nu, F* nv)
{
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

            F Hij = ef[ij] + D[ij];

            if ((D[ij] < m_grnd) || (Hij <= 0.0) || (fabs(ef[ij]) / Hij <= m_bound)
                || ((step % m_period) != 0)) {
                ne[ij] = ef[ij];
                nu[ij] = uf[ij];
                nv[ij] = vf[ij];
                continue;
            }

            F Himj = ef[imj] + D[imj];
            F Hijm = ef[ijm] + D[ijm];
            F Hipj = ef[ipj] + D[ipj];
            F Hijp = ef[ijp] + D[ijp];

            F Ae = 0.0f, Au = 0.0f, Av = 0.0f;
            int C = 0;

            if ((i > 0) && (D[imj] >= m_grnd) && (Himj > 0.0)) {
                Ae += ef[imj];
                Au += uf[imj];
                Av += vf[imj];
                C++;
            }
            if ((i < nx - 1) && (D[ipj] >= m_grnd) && (Hipj > 0.0)) {
                Ae += ef[ipj];
                Au += uf[ipj];
                Av += vf[ipj];
                C++;
            }
            if ((j > 0) && (D[ijm] >= m_grnd) && (Hijm > 0.0)) {
                Ae += ef[ijm];
                Au += uf[ijm];
                Av += vf[ijm];
                C++;
            }
            if ((j < ny - 1) && (D[ijp] >= m_grnd) && (Hijp > 0.0)) {
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

            ne[ij] = (1.0 - m_coeff) * ef[ij] + m_coeff * Ae / C;
            nu[ij] = (1.0 - m_coeff) * uf[ij] + m_coeff * Au / C;
            nv[ij] = (1.0 - m_coeff) * vf[ij] + m_coeff * Av / C;

            if (D[ij] + ne[ij] < 0.0) {
                ne[ij] = -D[ij] + 0.01f;
                nu[ij] = 0.0f;
                nv[ij] = 0.0f;
            }
        }
    }
}

template <class F>
LimitWaveOnDeep<F>::LimitWaveOnDeep(int period, F deep, F grnd)
    : DataProcessor<F>(true)
    , m_period(period)
    , m_deep(deep)
    , m_grnd(grnd)
{
}

template <class F>
void LimitWaveOnDeep<F>::operator()(
    int nx, int ny, int step, F dT, F* e, F* u, F* v, F* D, F* ne, F* nu, F* nv)
{
    if ((step % m_period) != 0)
        return;

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            int ipj = j * nx + i + 1;
            int ijp = (j + 1) * nx + i;
            int ipjp = (j + 1) * nx + i + 1;
            int imj = j * nx + i - 1;
            int ijm = (j - 1) * nx + i;
            int imjm = (j - 1) * nx + i - 1;
            int ipjm = (j - 1) * nx + i + 1;
            int imjp = (j + 1) * nx + i - 1;

            ne[ij] = e[ij];
            nu[ij] = u[ij];
            nv[ij] = v[ij];

            if ((i == 0) || (j == 0) || (i == nx - 1) || (j == ny - 1))
                continue;

            if (D[ij] < m_deep)
                continue;

            F cDij = sqrt(9.81 * D[ij]);

            F c = std::max((F)sqrt(9.81 * (ne[ij] + D[ij])), cDij);
            if (std::abs(nu[ij]) > c) {
                F m = nu[ij];
                m = ((D[imjm] >= m_grnd) && (std::abs(u[imjm]) < std::abs(m))) ? u[imjm] : m;
                m = ((D[ijm] >= m_grnd) && (std::abs(u[ijm]) < std::abs(m))) ? u[ijm] : m;
                m = ((D[ipjm] >= m_grnd) && (std::abs(u[ipjm]) < std::abs(m))) ? u[ipjm] : m;
                m = ((D[imj] >= m_grnd) && (std::abs(u[imj]) < std::abs(m))) ? u[imj] : m;
                m = ((D[ipj] >= m_grnd) && (std::abs(u[ipj]) < std::abs(m))) ? u[ipj] : m;
                m = ((D[imjp] >= m_grnd) && (std::abs(u[imjp]) < std::abs(m))) ? u[imjp] : m;
                m = ((D[ijp] >= m_grnd) && (std::abs(u[ijp]) < std::abs(m))) ? u[ijp] : m;
                m = ((D[ipjp] >= m_grnd) && (std::abs(u[ipjp]) < std::abs(m))) ? u[ipjp] : m;
                nu[ij] = m;

                m = ne[ij];
                m = ((D[imjm] >= m_grnd) && (std::abs(e[imjm]) < std::abs(m))) ? e[imjm] : m;
                m = ((D[ijm] >= m_grnd) && (std::abs(e[ijm]) < std::abs(m))) ? e[ijm] : m;
                m = ((D[ipjm] >= m_grnd) && (std::abs(e[ipjm]) < std::abs(m))) ? e[ipjm] : m;
                m = ((D[imj] >= m_grnd) && (std::abs(e[imj]) < std::abs(m))) ? e[imj] : m;
                m = ((D[ipj] >= m_grnd) && (std::abs(e[ipj]) < std::abs(m))) ? e[ipj] : m;
                m = ((D[imjp] >= m_grnd) && (std::abs(e[imjp]) < std::abs(m))) ? e[imjp] : m;
                m = ((D[ijp] >= m_grnd) && (std::abs(e[ijp]) < std::abs(m))) ? e[ijp] : m;
                m = ((D[ipjp] >= m_grnd) && (std::abs(e[ipjp]) < std::abs(m))) ? e[ipjp] : m;
                ne[ij] = m;
            }

            c = std::max((F)sqrt(9.81 * (ne[ij] + D[ij])), cDij);
            if (std::abs(nu[ij]) > c)
                nu[ij] = ((nu[ij] < 0) ? -1.0 : 1.0) * c * std::abs(ne[ij])
                    / std::max(D[ij], D[ij] + ne[ij]);
        }
    }

#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            int ipj = j * nx + i + 1;
            int ijp = (j + 1) * nx + i;
            int ipjp = (j + 1) * nx + i + 1;
            int imj = j * nx + i - 1;
            int ijm = (j - 1) * nx + i;
            int imjm = (j - 1) * nx + i - 1;
            int ipjm = (j - 1) * nx + i + 1;
            int imjp = (j + 1) * nx + i - 1;

            e[ij] = ne[ij];
            u[ij] = nu[ij];
            v[ij] = nv[ij];

            if ((i == 0) || (j == 0) || (i == nx - 1) || (j == ny - 1))
                continue;

            if (D[ij] < m_deep)
                continue;

            F cDij = sqrt(9.81 * D[ij]);

            F c = std::max((F)sqrt(9.81 * (e[ij] + D[ij])), cDij);
            if (std::abs(v[ij]) > c) {
                F m = u[ij];
                m = ((D[imjm] >= m_grnd) && (std::abs(nv[imjm]) < std::abs(m))) ? nv[imjm] : m;
                m = ((D[ijm] >= m_grnd) && (std::abs(nv[ijm]) < std::abs(m))) ? nv[ijm] : m;
                m = ((D[ipjm] >= m_grnd) && (std::abs(nv[ipjm]) < std::abs(m))) ? nv[ipjm] : m;
                m = ((D[imj] >= m_grnd) && (std::abs(nv[imj]) < std::abs(m))) ? nv[imj] : m;
                m = ((D[ipj] >= m_grnd) && (std::abs(nv[ipj]) < std::abs(m))) ? nv[ipj] : m;
                m = ((D[imjp] >= m_grnd) && (std::abs(nv[imjp]) < std::abs(m))) ? nv[imjp] : m;
                m = ((D[ijp] >= m_grnd) && (std::abs(nv[ijp]) < std::abs(m))) ? nv[ijp] : m;
                m = ((D[ipjp] >= m_grnd) && (std::abs(nv[ipjp]) < std::abs(m))) ? nv[ipjp] : m;
                v[ij] = m;

                m = ne[ij];
                m = ((D[imjm] >= m_grnd) && (std::abs(ne[imjm]) < std::abs(m))) ? ne[imjm] : m;
                m = ((D[ijm] >= m_grnd) && (std::abs(ne[ijm]) < std::abs(m))) ? ne[ijm] : m;
                m = ((D[ipjm] >= m_grnd) && (std::abs(ne[ipjm]) < std::abs(m))) ? ne[ipjm] : m;
                m = ((D[imj] >= m_grnd) && (std::abs(ne[imj]) < std::abs(m))) ? ne[imj] : m;
                m = ((D[ipj] >= m_grnd) && (std::abs(ne[ipj]) < std::abs(m))) ? ne[ipj] : m;
                m = ((D[imjp] >= m_grnd) && (std::abs(ne[imjp]) < std::abs(m))) ? ne[imjp] : m;
                m = ((D[ijp] >= m_grnd) && (std::abs(ne[ijp]) < std::abs(m))) ? ne[ijp] : m;
                m = ((D[ipjp] >= m_grnd) && (std::abs(ne[ipjp]) < std::abs(m))) ? ne[ipjp] : m;
                e[ij] = m;
            }

            c = std::max((F)sqrt(9.81 * (e[ij] + D[ij])), cDij);
            if (std::abs(v[ij]) > c)
                v[ij] = ((v[ij] < 0) ? -1.0 : 1.0) * c * std::abs(e[ij])
                    / std::max(D[ij], D[ij] + e[ij]);
        }
    }
}

template <class F>
RemoveSingleIslands<F>::RemoveSingleIslands(F grnd)
    : DataProcessor<F>(true)
    , m_grnd(grnd)
{
}

template <class F>
void RemoveSingleIslands<F>::operator()(
    int nx, int ny, int step, F dT, F* e, F* u, F* v, F* D, F* ne, F* nu, F* nv)
{
    if (step != 0)
        return;

    int cnt = 0;

    for (int j = 0; j < ny; j++) {
        int rcnt = 0;
        for (int i = 0; i < nx; i++) {
            int ij = j * nx + i;
            int ipj = j * nx + i + 1;
            int ijp = (j + 1) * nx + i;
            int ipjp = (j + 1) * nx + i + 1;
            int imj = j * nx + i - 1;
            int ijm = (j - 1) * nx + i;
            int imjm = (j - 1) * nx + i - 1;
            int ipjm = (j - 1) * nx + i + 1;
            int imjp = (j + 1) * nx + i - 1;

            if (D[ij] >= m_grnd)
                continue;

            if ((D[ijm] < m_grnd) && (D[ijp] < m_grnd) && (D[imj] < m_grnd) && (D[ipj] < m_grnd)) {
                D[ij] = m_grnd + 0.01;
                rcnt++;
            }
        }
        cnt += rcnt;
    }
    std::cout << "RemoveSingleIslands: removed " << cnt << " islands" << std::endl;
}

template class DataProcessor<float>;
template class DataProcessor<double>;
}
