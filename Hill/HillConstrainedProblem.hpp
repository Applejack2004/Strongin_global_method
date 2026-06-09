#pragma once

#include "IConstrainedOptProblem.hpp"
#include "HillProblem.hpp"

class HillConstrainedProblem : public IConstrainedOptProblem
{
protected:
    vector<IOptProblem*> functions;
    int functionCount;

    virtual double Compute(int index, const vector<double>& y) const;
    virtual vector<double> ComputeDerivatives(int index, const vector<double>& y) const;

public:
    HillConstrainedProblem(IOptProblem* objectiveFunc, int numberOfConstraints,
        EConstrainedProblemType problemType, double fraction, int activeConstrNum);

    virtual ~HillConstrainedProblem();
};