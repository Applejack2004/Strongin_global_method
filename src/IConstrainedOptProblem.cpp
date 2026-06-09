#include "IConstrainedOptProblem.hpp"
#include <algorithm>
#include <cmath>

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::MaxFunctionCalculate(vector<double> y)
{
    if (mConstraintIndeces.size() > 0)
    {
        vector<double> x = TransformConstraintPoint(y, 0);

        double res = Compute(mConstraintIndeces[0], x);
        for (uint i = 1; i < mConstraintIndeces.size(); i++)
        {
            x = TransformConstraintPoint(y, i);

            double f = Compute(mConstraintIndeces[i], x);
            if (res < f)
                res = f;
        }
        return res;
    }
    return 0;
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::CalculateRHS(double delta, int m, double Epsilon,
    int maxM)
{
    double rhs = 0;
    double hmin = mOptimumValue;
    double hmax = hmin;
    double d = 0;

    // multidimensional grid
    int* size = new int[mDimension]; // number of points to one dimension
    double* step = new double[mDimension]; // step ofthe grid
    int sumn = 1; // number of trials

    double* a = mLoBound.data();
    double* b = mUpBound.data();

    double multiplyingLength = 1;
    for (int i = 0; i < mDimension; i++)
    {
        d = (b[i] - a[i]);
        size[i] = (int)ceil(d / Epsilon) + 1;
        step[i] = d / (size[i] - 1);
        multiplyingLength = multiplyingLength * d;
        sumn *= (size[i]);
    }

    if ((sumn > maxM) || (sumn <= 0))
    {
        multiplyingLength = multiplyingLength / maxM;
        Epsilon = pow(multiplyingLength, 1.0 / (double)mDimension);
        sumn = 1;
        multiplyingLength = 1;

        for (int i = 0; i < mDimension; i++)
        {
            d = (b[i] - a[i]);
            size[i] = (int)ceil(d / Epsilon) + 1;
            step[i] = d / (size[i] - 1);
            sumn *= (size[i]);
        }
    }

    double* f = new double[sumn]; // function value
    vector<double> yArray(mDimension);

    for (int i = 0; i < sumn; i++)
    {
        double w;
        int z = i;
        // compute coordinates of trial point
        for (int j = 0; j < mDimension; j++)
        {
            w = z % size[j]; // node number
            yArray[j] = a[j] + w * step[j];//left border + node number * step
            z = z / size[j]; // for the next loop iteration
        }
        // carry out the trial
        f[i] = MaxFunctionCalculate(yArray);
        hmax = std::max(f[i], hmax);
        hmin = std::min(f[i], hmin);
    }

    double* h1 = new double[m];
    double* h2 = new double[m];
    int* p = new int[m];
    int* s = new int[m];

    double deltah = (hmax - hmin) / m;

    for (int i = 0; i < m; i++)
    {
        h1[i] = hmin + i * deltah;
        h2[i] = hmin + (i + 1) * deltah;
        p[i] = 0;
        s[i] = 0;
    }

    for (int i = 0; i < sumn; i++)
        for (int j = 0; j < m; j++)
            if ((f[i] >= h1[j]) && (f[i] <= h2[j]))
            {
                p[j]++;
                break;
            }

    s[0] = p[0];
    for (int i = 1; i < m; i++)
    {
        s[i] = s[i - 1] + p[i];
    }

    double smax = s[m - 1];
    double g = delta * smax;
    for (int i = 0; i < m; i++)
    {
        if (s[i] >= g)
        {
            rhs = h2[i];
            break;
        }
    }

    double dm = delta;
    if (dm == 0)
        dm += 0.1;
    dm = dm * (hmax - hmin);

    double criticalValue = MaxFunctionCalculate(mOptimumPoint);

    if (rhs < criticalValue)
    {
        rhs = criticalValue + dm;
    }

    delete[] size;
    delete[] step;
    delete[] f;

    delete[] h1;
    delete[] h2;
    delete[] p;
    delete[] s;

    return rhs;
}

// ------------------------------------------------------------------------------------------------
void IConstrainedOptProblem::SetZoom()
{
    double AccuracyDouble = 0.00000001;
    double* lower = mLoBound.data();
    double* upper = mUpBound.data();
    double* constraintMin = 0;

    double maxDistanceToBoundary = 0;

    for (int k = 0; k < mDimension; k++)
    {
        if (maxDistanceToBoundary < (mOptimumPoint[k] - lower[k]))
            maxDistanceToBoundary = (mOptimumPoint[k] - lower[k]);
        if (maxDistanceToBoundary < (upper[k] - mOptimumPoint[k]))
            maxDistanceToBoundary = (upper[k] - mOptimumPoint[k]);
    }

    if (fabs(maxDistanceToBoundary) < AccuracyDouble)
    {
        mIsZoom = false;
        for (uint j = 0; j < mConstraintIndeces.size(); j++)
        {
            mZoomRatios[j] = 1;
        }

        return;
    }

    for (uint i = 0; i < mConstraintIndeces.size(); i++)
    {
        double minDistanceToBoundary = upper[0] - lower[0];

        constraintMin = mFunctions[mConstraintIndeces[i]].mMinimumPoint.data();

        for (int k = 0; k < mDimension; k++)
        {
            if (minDistanceToBoundary > (constraintMin[k] - lower[k]))
                minDistanceToBoundary = (constraintMin[k] - lower[k]);
            if (minDistanceToBoundary > (upper[k] - constraintMin[k]))
                minDistanceToBoundary = (upper[k] - constraintMin[k]);
        }

        if (fabs(minDistanceToBoundary) < AccuracyDouble)
        {
            mIsZoom = false;
            for (uint j = 0; j < mConstraintIndeces.size(); j++)
            {
                mZoomRatios[j] = 1;
            }

            return;
        }
        else
        {
            mZoomRatios[i] = maxDistanceToBoundary / minDistanceToBoundary;
        }
    }
}

// ------------------------------------------------------------------------------------------------
void IConstrainedOptProblem::SetShift()
{
    double* constraintMin = 0;

    for (uint i = 0; i < mConstraintIndeces.size(); i++)
    {
        constraintMin = mFunctions[mConstraintIndeces[i]].mMinimumPoint.data();

        for (int k = 0; k < mDimension; k++)
        {
            mShift[i][k] = mOptimumPoint[k] - constraintMin[k] * mZoomRatios[i];
        }
    }
}

// ------------------------------------------------------------------------------------------------
void IConstrainedOptProblem::InitProblem()
{
    // --- Стандартная инициализация ---
    mQ.clear();
    mZoomRatios.clear();
    mShift.clear();
    mImprovementCoefficients.clear();
    if (mActiveConstraintNumber == 0) mActiveConstraintNumber = (int)mConstraintIndeces.size();
    for (uint j = 0; j < mConstraintIndeces.size(); j++) mImprovementCoefficients.push_back(10.0);
    mQ.resize(mConstraintIndeces.size());
    mZoomRatios.resize(mConstraintIndeces.size());
    mShift.resize(mConstraintIndeces.size());
    for (uint i = 0; i < mConstraintIndeces.size(); i++) {
        mZoomRatios[i] = 1.0; mQ[i] = 0.0;
        mShift[i].assign(mDimension, 0.0);
    }
    if (mIsZoom) SetZoom();
    if (mIsShift) SetShift();

    // 1. Расчет базового порога для соблюдения площади (30% или 50%)
    if (mIsTotalDelta)
    {
        double q_base = CalculateRHS(mFeasibleDomainFraction);
        for (uint i = 0; i < mConstraintIndeces.size(); i++)
            mQ[i] = q_base;

        // 2. ГЕОМЕТРИЧЕСКИЙ СДВИГ (Для типов OnBorder и OutDomain)
        if ((mProblemType == cptOnFeasibleBorder || mProblemType == cptOutFeasibleDomain) && mConstraintIndeces.size() > 0)
        {
            double x0 = mOptimumPoint[0];
            double y0 = mOptimumPoint[1];
            double best_dx = 0, best_dy = 0;
            double min_dist = 1e100;
            bool found_any = false;

            auto check_feasible = [&](const vector<double>& p) {
                for (uint i = 0; i < mConstraintIndeces.size(); i++) {
                    if (ComputeConstraint(i, p) > 0) return false;
                }
                return true;
                };

            bool start_feasible = check_feasible(mOptimumPoint);

            // Поиск ближайшей границы
            for (int i = 0; i < 360; i++) {
                double angle = i * (3.14159265358979 / 180.0);
                double vx = cos(angle), vy = sin(angle);
                for (double t = 0.001; t < 2.0; t += 0.005) {
                    vector<double> p_test = { x0 + vx * t, y0 + vy * t };
                    if (p_test[0] < -1.0 || p_test[0] > 1.0 || p_test[1] < -1.0 || p_test[1] > 1.0) break;
                    if (check_feasible(p_test) != start_feasible) {
                        double t_low = t - 0.005, t_high = t;
                        for (int b = 0; b < 40; b++) {
                            double t_mid = (t_low + t_high) / 2.0;
                            if (check_feasible({ x0 + vx * t_mid, y0 + vy * t_mid }) == start_feasible) t_low = t_mid;
                            else t_high = t_mid;
                        }
                        if (t_high < min_dist) {
                            min_dist = t_high;
                            best_dx = vx * t_high;
                            best_dy = vy * t_high;
                            found_any = true;
                        }
                        break;
                    }
                }
            }

            if (found_any) {
                // ПРАВКА ЛОГИКИ СДВИГА:
                double multiplier = 1.0;
                if (mProblemType == cptOutFeasibleDomain) {
                    // Если минимум должен быть ВНЕ области, мы двигаем забор "глубже" в синюю зону,
                    // чтобы точка x0 осталась в серой зоне.
                    multiplier = 1.15; // Сдвигаем на 15% дальше границы
                }

                for (uint i = 0; i < mConstraintIndeces.size(); i++) {
                    mShift[i][0] -= (best_dx * multiplier);
                    mShift[i][1] -= (best_dy * multiplier);
                }
            }
        }
    }

    // 3. ОБНОВЛЕНИЕ ИСТИННОГО ЗНАЧЕНИЯ ОПТИМУМА
    if (mProblemType == cptOutFeasibleDomain) {
        // Ландшафт уже сформирован. Безусловный минимум сейчас находится вне допустимой области.
        // Запускаем перебор, чтобы найти настоящий условный минимум внутри области/на границе
        // и сделать его новой мишенью для критерия остановки.
        FindTrueConstrainedOptimum();
    }
    else {
        // Для InDomain и OnBorder старая логика остается нетронутой
        mOptimumValue = ComputeFunction(mOptimumPoint);
    }
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::TransformValue(double val, vector<double> point, int index) const
{
    double res = val;
    double resultCoefficient = 0;

    if (mIsImprovementOfTheObjective)
    {
        for (uint j = 0; j < mConstraintIndeces.size(); j++)
        {
            double val = Compute(mConstraintIndeces[j],
                TransformConstraintPoint(point, j));
            val = TransformConstraintValue(val, j);

            double fVal = std::max(val, 0.0);
            resultCoefficient += mImprovementCoefficients[j] * (fVal * fVal * fVal);
        }
    }

    return res - resultCoefficient;
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::TransformPoint(vector<double> point, int index) const
{
    return point; // Больше не искажаем пространство!
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::TransformConstraintValue(double val, int index) const
{
    double res;

    res = val - mQ[index];

    return res;
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::TransformConstraintPoint(vector<double> point,
    int index) const
{
    vector<double> res(mDimension);

    for (int i = 0; i < mDimension; i++)
    {
        res[i] = (point[i] - mShift[index][i]) /
            mZoomRatios[index];
    }

    return res;
}

// ------------------------------------------------------------------------------------------------
IConstrainedOptProblem::IConstrainedOptProblem(int dim, vector<double> loBound,
    vector<double> upBound, int probIndex,
    EConstrainedProblemType problemType, double fraction, int activeConstrNum)
    : IGeneralOptProblem(dim, loBound, upBound, probIndex)
{
    mFunctions.push_back({ mDimension, {}, 0, {}, 0,
      -1, true, false, false, false });
    mFunctionNumber = mFunctions.size();
    mCriterionIndeces.push_back(0);
    mProblemType = problemType;
    mFeasibleDomainFraction = fraction;
    mActiveConstraintNumber = activeConstrNum;

    if (problemType == cptInFeasibleDomain)
    {
        mIsShift = true;
        mIsZoom = true;
        mIsOnBoundary = false;
        mIsImprovementOfTheObjective = false;
        mIsTotalDelta = true;
    }
    else if (problemType == cptOutFeasibleDomain)
    {
        mIsShift = true;
        mIsZoom = true;
        mIsOnBoundary = true; // Чтобы условный минимум лежал точно на границе!
        mIsImprovementOfTheObjective = false; // А безусловный "провалится" снаружи
        mIsTotalDelta = true;
    }
    else if (problemType == cptOnFeasibleBorder)
    {
        mIsShift = true;
        mIsZoom = true;
        mIsOnBoundary = true; // Минимум лежит точно на границе
        mIsImprovementOfTheObjective = false;
        mIsTotalDelta = true;
    }
    else
    {
        mIsShift = false;
        mIsZoom = false;
        mIsOnBoundary = false;
        mIsImprovementOfTheObjective = false;
        mIsTotalDelta = false;
    }
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::GetOptimumPoint() const
{
    if (mProblemType == cptNormal)
    {
        if (mOptimumPoint.empty())
            throw string("Minimum point is unknown");
        for (uint j = 0; j < mConstraintIndeces.size(); j++)
            if (Compute(mConstraintIndeces[j], mOptimumPoint) > 0)
                throw string("Minimum point is unknown");
    }
    return TransformPoint(IGeneralOptProblem::GetOptimumPoint(), 0);
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::GetOptimumValue() const
{
    if (mProblemType == cptNormal)
    {
        if (mOptimumPoint.empty())
            throw string("Minimum point is unknown");
        for (uint j = 0; j < mConstraintIndeces.size(); j++)
            if (Compute(mConstraintIndeces[j], mOptimumPoint) > 0)
                throw string("Minimum point is unknown");
    }
    return TransformValue(IGeneralOptProblem::GetOptimumValue(),
        IGeneralOptProblem::GetOptimumPoint(), 0);
}

// ------------------------------------------------------------------------------------------------
bool IConstrainedOptProblem::GetStatus(enum EOptFunctionParameter param) const
{
    return IGeneralOptProblem::GetStatus(mCriterionIndeces[0], param);
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::GetMaxPoint() const
{
    vector<double> point = IGeneralOptProblem::GetMaxPoint(mCriterionIndeces[0]);
    return TransformPoint(point, 0);
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::GetMaxValue() const
{
    double val = IGeneralOptProblem::GetMaxValue(mCriterionIndeces[0]);
    return TransformValue(val, GetMaxPoint(), 0);
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::GetLipschitzConstant() const
{
    double lip = IGeneralOptProblem::GetLipschitzConstant(0);
    return lip;
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::ComputeFunction(const vector<double>& y) const
{
    double val = Compute(mCriterionIndeces[0],
        TransformPoint(y, 0));
    return TransformValue(val, y, 0);
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::ComputeFunctionDerivatives(const vector<double>& y) const
{
    vector<double> res = ComputeDerivatives(mCriterionIndeces[0], y);
    return res;
}

// ------------------------------------------------------------------------------------------------
int IConstrainedOptProblem::GetConstraintsNumber() const
{
    return mConstraintIndeces.size();
}

// ------------------------------------------------------------------------------------------------
bool IConstrainedOptProblem::GetConstraintStatus(int index, enum EOptFunctionParameter param) const
{
    return IGeneralOptProblem::GetStatus(mConstraintIndeces[index], param);
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::GetConstraintLipschitzConstant(int index) const
{
    double val = IGeneralOptProblem::GetLipschitzConstant(index);
    return val;
}

// ------------------------------------------------------------------------------------------------
double IConstrainedOptProblem::ComputeConstraint(int index, const vector<double>& y) const
{
    double val = Compute(mConstraintIndeces[index],
        TransformConstraintPoint(y, index));
    return TransformConstraintValue(val, index);
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::ComputeConstraints(const vector<double>& y,
    EConstraintComputationType t, int& index) const
{
    vector<double> res;
    switch (t)
    {
    case cctAllConstraints:
        for (uint i = 0; i < mConstraintIndeces.size(); i++)
        {
            res.push_back(ComputeConstraint(i, y));
        }
        index = mConstraintIndeces.size() - 1;
        break;
    case cctIndexScheme:
        double tmp;
        uint i;
        for (i = 0; i < mConstraintIndeces.size(); i++)
        {
            tmp = ComputeConstraint(i, y);
            if (tmp >= 0)
            {
                res.push_back(tmp);
                break;
            }
        }
        index = i;
        break;
    default:
        throw string("Unknown type of constraints computation");
    }
    return res;
}

// ------------------------------------------------------------------------------------------------
vector<double> IConstrainedOptProblem::ComputeConstraintDerivatives(int index,
    const vector<double>& y) const
{
    vector<double> res = ComputeDerivatives(mConstraintIndeces[index], y);
    return res;
}

// ------------------------------------------------------------------------------------------------
IConstrainedOptProblem::~IConstrainedOptProblem()
{
}

// ------------------------------------------------------------------------------------------------
void IConstrainedOptProblem::FindTrueConstrainedOptimum()
{
    double best_val = 1e100;
    vector<double> best_point = mOptimumPoint;
    bool found = false;

    // Этап 1: Равномерный перебор по крупной сетке (суммарно ~250 000 точек для любой размерности)
    int pts_per_dim = std::max(10, (int)std::pow(250000.0, 1.0 / mDimension));
    long long total_pts = 1;
    for (int d = 0; d < mDimension; d++) total_pts *= pts_per_dim;

    vector<double> pt(mDimension);
    for (long long i = 0; i < total_pts; ++i) {
        long long temp = i;
        for (int d = 0; d < mDimension; ++d) {
            int idx = temp % pts_per_dim;
            temp /= pts_per_dim;
            pt[d] = mLoBound[d] + idx * (mUpBound[d] - mLoBound[d]) / (pts_per_dim - 1);
        }

        // Проверяем, находится ли точка в допустимой области
        bool feasible = true;
        for (uint k = 0; k < mConstraintIndeces.size(); k++) {
            if (ComputeConstraint(k, pt) > 0) {
                feasible = false;
                break; // Нарушено ограничение
            }
        }

        // Если допустима, считаем целевую функцию
        if (feasible) {
            double val = ComputeFunction(pt);
            if (val < best_val) {
                best_val = val;
                best_point = pt;
                found = true;
            }
        }
    }

    // Этап 2: Локальное уточнение (мелкая сетка вокруг найденной точки)
    if (found) {
        double local_best_val = best_val;
        vector<double> local_best_point = best_point;

        // Суммарно ~10 000 точек вокруг окрестности
        int fine_pts = std::max(10, (int)std::pow(10000.0, 1.0 / mDimension));
        long long total_fine = 1;
        for (int d = 0; d < mDimension; d++) total_fine *= fine_pts;

        for (long long i = 0; i < total_fine; ++i) {
            long long temp = i;
            for (int d = 0; d < mDimension; ++d) {
                int idx = temp % fine_pts;
                temp /= fine_pts;

                // Окрестность +/- 2 шага крупной сетки
                double step = (mUpBound[d] - mLoBound[d]) / (pts_per_dim - 1);
                double min_bound = std::max(mLoBound[d], best_point[d] - 2.0 * step);
                double max_bound = std::min(mUpBound[d], best_point[d] + 2.0 * step);

                pt[d] = min_bound + idx * (max_bound - min_bound) / (fine_pts - 1);
            }

            bool feasible = true;
            for (uint k = 0; k < mConstraintIndeces.size(); k++) {
                if (ComputeConstraint(k, pt) > 0) {
                    feasible = false;
                    break;
                }
            }

            if (feasible) {
                double val = ComputeFunction(pt);
                if (val < local_best_val) {
                    local_best_val = val;
                    local_best_point = pt;
                }
            }
        }

        // КРИТИЧЕСКИЙ МОМЕНТ: Подменяем координаты оптимума!
        mOptimumPoint = local_best_point;
        mOptimumValue = local_best_val;
    }
}
