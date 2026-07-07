#pragma once
#include <PlasticConstraint/BoxLagrangianConstraint.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/type/Vec.h>


namespace plasticconstraint::constraint
{

using sofa::core::objectmodel::ComponentState;
using sofa::core::VecCoordId;


template<class DataTypes>
BoxLagrangianConstraint<DataTypes>::BoxLagrangianConstraint(MechanicalState* object)
    : Inherit(object)
    , d_indices(initData(&d_indices, "indices", "indices of the stop constraint"))
    , d_min(initData(&d_min, -100.0_sreal, "min", "minimum value accepted"))
    , d_max(initData(&d_max, 100.0_sreal, "max", "maximum value accepted"))
    , d_force(initData(&d_force, "force", "plastic strain"))
{
    if (d_indices.getValue().empty())
    {
        sofa::type::vector<unsigned int> all;
        for (unsigned int i = 0; i < 389; ++i)
            all.push_back(i);
        d_indices.setValue(all);
    }


    const auto& indices = d_indices.getValue();

    msg_info("BoxLagrangianConstraint") << " am here ";
    msg_info("BoxLagrangianConstraint") << " indices " << indices.size();

    for (size_t i = 0; i < indices.size(); ++i){
        _oldlambda.push_back(0.0);
    }
    d_force.setValue(_oldlambda);
    msg_info("BoxLagrangianConstraint") << "old lambda " << _oldlambda;

}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::buildConstraintMatrix(const sofa::core::ConstraintParams* /*cParams*/, DataMatrixDeriv &c_d, unsigned int &cIndex, const DataVecCoord &/*x*/)
{
    auto c = sofa::helper::getWriteAccessor(c_d);

    if (d_indices.getValue().empty())
    {
        sofa::type::vector<unsigned int> all;
        for (unsigned int i = 0; i < 389; ++i)
            all.push_back(i);
        d_indices.setValue(all);
    }
    const auto& indices = d_indices.getValue();

    for (const auto idx : indices)
    {
        const double h = this->getContext()->getDt();
        MatrixDerivRowIterator c_it = c->writeLine(cIndex++);
        //c_it.setCol(idx, Coord(1,0,1)); //pour tester avec vec3
        c_it.setCol(idx, Coord(1,0,0,0,0,0));
    }
}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* /*cParams*/, sofa::linearalgebra::BaseVector *resV, const DataVecCoord &x, const DataVecDeriv &/*v*/)
{
    const auto constraintIndex = this->d_constraintIndex.getValue();
    const auto& indices = d_indices.getValue();
    const auto& xValue = x.getValue();
    
    for (size_t i = 0; i < indices.size(); ++i)
    {
        resV->set(constraintIndex + i, xValue[indices[i]][0]);
    }
}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::getConstraintResolution(const sofa::core::ConstraintParams *, std::vector<sofa::core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{
    const auto& indices = d_indices.getValue();
    
    for (size_t i = 0; i < indices.size(); ++i)
        resTab[offset++] = new BoxLagrangianConstraintResolution1Dof(d_min.getValue(), d_max.getValue(), d_force.getValue()[i]);
}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::storeLambda(const sofa::core::ConstraintParams* cParams,
                                        sofa::core::MultiVecDerivId res,
                                        const sofa::linearalgebra::BaseVector* lambda)
{

    const auto& indices = d_indices.getValue();
    sofa::type::vector<double> tmp;
    sofa::type::vector<double> deltalambda;
    deltalambda.resize(indices.size());
    
    for (size_t i = 0; i < indices.size(); ++i){

        deltalambda[i] = lambda->element(this->d_constraintIndex.getValue() + i) - _oldlambda[i];
        // msg_info("BoxLagrangianConstraint") << "constraint index" << this->d_constraintIndex.getValue();
        tmp.push_back(d_force.getValue()[i] + deltalambda[i]);
        _oldlambda[i] = lambda->element(this->d_constraintIndex.getValue() + i);
    }
    d_force.setValue(tmp);
}

} //namespace sofa::component::constraint::lagrangian::model
