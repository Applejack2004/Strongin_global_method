#define _USE_MATH_DEFINES

#include "ShekelConstrainedProblem.hpp"
#include "ShekelProblemFamily.hpp"
#include <math.h>
#include <algorithm>
#include <cmath>
#include <cassert>

// ------------------------------------------------------------------------------------------------
TShekelConstrainedProblem::TShekelConstrainedProblem(IOptProblem* objectiveFunc, int numberOfConstraints,
    EConstrainedProblemType problemType, double fraction, int activeConstrNum)
    : IConstrainedOptProblem(objectiveFunc->GetDimension(), {}, {}, -1, problemType,
        fraction, activeConstrNum)
{
    assert(objectiveFunc != nullptr);

    functions.push_back(objectiveFunc);

    for (int i = 0; i < numberOfConstraints; i++)
    {
        functions.push_back(new TShekelProblem(i));
    }
    functionCount = (int)functions.size();

    mProblemIndex = objectiveFunc->GetProblemIndex();

    mDimension = objectiveFunc->GetDimension();
    objectiveFunc->GetBounds(mLoBound, mUpBound);
    mOptimumPoint = objectiveFunc->GetOptimumPoint();
    mOptimumValue = objectiveFunc->GetOptimumValue();

    mFunctions.clear();
    mFunctions.resize(functionCount);
    mFunctions[0].mDimension = mDimension;
    mFunctions[0].mMinimumPoint = mOptimumPoint;
    mFunctions[0].mMinimumValue = mOptimumValue;
    mFunctions[0].mIsMinimumKnown = true;
    mFunctions[0].mIsDerivativesKnown = true;

    for (int i = 1; i < functionCount; i++)
    {
        mFunctions[i] = { mDimension, functions[i]->GetOptimumPoint(),
                          functions[i]->GetOptimumValue(),
                          {}, 0, -1,
                          true, false, false, true };
    }

    mFunctionNumber = (int)mFunctions.size();
    mCriterionIndeces.clear();
    mConstraintIndeces.clear();
    mCriterionIndeces.push_back(0);

    for (int i = 1; i < functionCount; i++)
    {
        mConstraintIndeces.push_back(i);
    }

    InitProblem();
}

// ------------------------------------------------------------------------------------------------
TShekelConstrainedProblem::~TShekelConstrainedProblem()
{
    for (size_t i = 1; i < functions.size(); ++i)
    {
        delete functions[i];
    }
    functions.clear();
}

// ------------------------------------------------------------------------------------------------
double TShekelConstrainedProblem::Compute(int index, const vector<double>& y) const
{
    assert(index >= 0 && index < (int)functions.size());
    return functions[index]->ComputeFunction(y);
}

// ------------------------------------------------------------------------------------------------
vector<double> TShekelConstrainedProblem::ComputeDerivatives(int index, const vector<double>& y) const
{
    assert(index >= 0 && index < (int)functions.size());
    return functions[index]->ComputeFunctionDerivatives(y);
}