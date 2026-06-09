#pragma once

#include "IGeneralOptProblem.hpp"

class IOptProblem : public IGeneralOptProblem
{
protected:
    uint mFunctionIndex;
    void SetLipschitzConstant(double lipConst);
    void SetFunctionMax(vector<double> maxPoint, double maxValue);
    IOptProblem();
public:
    IOptProblem(int dim, vector<double> loBound, vector<double> upBound,
        vector<double> optPoint, double optVal, int probIndex = -1);

    virtual ~IOptProblem() = default;

    int GetProblemIndex() const { return mProblemIndex; }

    // Делаем метод GetBounds публичным, чтобы к нему был доступ
    void GetBounds(vector<double>& lb, vector<double>& ub) const
    {
        IGeneralOptProblem::GetBounds(lb, ub);
    }

    vector<double> GetOptimumPoint() const;
    double GetOptimumValue() const;
    bool GetStatus(enum EOptFunctionParameter param) const;
    vector<double> GetMaxPoint() const;
    double GetMaxValue() const;
    double GetLipschitzConstant() const;
    double ComputeFunction(const vector<double>& y) const;
    vector<double> ComputeFunctionDerivatives(const vector<double>& y) const;
};