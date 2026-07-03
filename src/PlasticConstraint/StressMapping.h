#pragma once
#include <PlasticConstraint/config.h.in>
#include <PlasticConstraint/LinearMapping.h>

#include <sofa/linearalgebra/EigenSparseMatrix.h>
#include <sofa/core/Mapping.h>
#include <sofa/core/Mapping.inl>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/vector.h>
#include <sofa/type/Mat.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/objectmodel/Link.h>
//#include <sofa/component/solidmechanics/fem/elastic/TetrahedronFEMForceField.h>

namespace sofa::component::mapping::linear
{

template <class TIn, class TOut>
class StressMapping : public LinearMapping<TIn, TOut>


{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(StressMapping,TIn,TOut), SOFA_TEMPLATE2(LinearMapping,TIn,TOut));
    //sofa::component::forcefield::TetrahedronFEMForceField<TIn>* m_femFF = nullptr;

    typedef LinearMapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;

    typedef typename In::Real			Real;
    typedef typename In::VecCoord		InVecCoord;
    typedef typename In::VecDeriv		InVecDeriv;
    typedef typename In::Coord			InCoord;
    typedef typename In::Deriv			InDeriv;
    typedef typename In::MatrixDeriv	InMatrixDeriv;
    typedef Data<InVecCoord> InDataVecCoord;
    typedef Data<InVecDeriv> InDataVecDeriv;

    typedef typename Out::VecCoord		VecCoord;
    typedef typename Out::VecDeriv		VecDeriv;
    typedef typename Out::Coord			Coord;
    typedef typename Out::Deriv			Deriv;
    typedef typename Out::MatrixDeriv	MatrixDeriv;

    typedef Out OutDataTypes;
    typedef typename OutDataTypes::Real     OutReal;
    typedef typename OutDataTypes::VecCoord OutVecCoord;
    typedef typename OutDataTypes::VecDeriv OutVecDeriv;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;

    enum
    {
        N = OutDataTypes::spatial_dimensions
    };
    enum
    {
        NIn = sofa::defaulttype::DataTypeInfo<InDeriv>::Size
    };
    enum
    {
        NOut = sofa::defaulttype::DataTypeInfo<Deriv>::Size
    };

    typedef type::Mat<N, N, Real> Mat;

protected:
    StressMapping()
        : Inherit()
        , l_inputTopology(initLink("inputTopology", "Link to the input topology"))
        , m_topology(nullptr)
        , d_youngModulus(initData(&d_youngModulus, {Real(5000)}, "youngModulus", "FEM Young's Modulus"))
        , d_poissonRatio(initData(&d_poissonRatio, {Real(0.45)}, "poissonRatio", "FEM Poisson Ratio"))

    {
        Js.resize( 1 );
        Js[0] = &J;
    }

    virtual ~StressMapping()
    {
    }

public:
    /// Return true if the destination model has the same topology as the source model.
    ///
    /// This is the case for mapping keeping a one-to-one correspondence between
    /// input and output DOFs (mostly identity or data-conversion mappings).
    bool sameTopology() const override { return true; }

    void init() override;

    void apply(const core::MechanicalParams* mparams, Data<VecCoord>& out, const Data<InVecCoord>& in) override; 

    void applyJ(const core::MechanicalParams *mparams, Data<VecDeriv>& out, const Data<InVecDeriv>& in) override;

    void applyJT(const core::MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<VecDeriv>& in) override;

    void applyJT(const core::ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<MatrixDeriv>& in) override;

    const sofa::linearalgebra::BaseMatrix* getJ() override;

    void handleTopologyChange() override;

protected:

    typedef linearalgebra::EigenSparseMatrix<TIn, TOut> eigen_type;
    eigen_type J;

    typedef type::vector< linearalgebra::BaseMatrix* > js_type;
    js_type Js;

    SingleLink<StressMapping<TIn, TOut>,
               sofa::core::topology::BaseMeshTopology,
               BaseLink::FLAG_STOREPATH> l_inputTopology;

    sofa::core::topology::BaseMeshTopology* m_topology;
    
    typedef type::VecNoInit<6,Real> VoigtTensor;

    /// structures used to compute vonMises stress
    typedef type::Mat<4, 4, Real> Mat44;
    typedef type::Mat<3, 3, Real> Mat33;

public:

    const js_type* getJs() override;
    Data<InVecCoord> d_initialPoints;
    Data<sofa::type::vector<Real>> d_youngModulus;
    Data<sofa::type::vector<Real>> d_poissonRatio;
    type::vector<Mat44> elemShapeFun;


};

#if !defined(PLASTICCONSTRAINT_STRESSMAPPING_CPP)

extern template class SOFA_COMPONENT_MAPPING_LINEAR_API StressMapping< defaulttype::Vec3Types, defaulttype::Vec6Types >;

#endif

} // namespace sofa::component::mapping::linear
