#pragma once
#include <PlasticConstraint/StressMapping.h>
#include <sofa/core/MappingHelper.h>
#include <sofa/core/visual/VisualParams.h>

#include <algorithm>  // std::sort

namespace sofa::component::mapping::linear
{

using sofa::core::objectmodel::ComponentState;
using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::core::visual::VisualParams;
using sofa::type::vector;

template<class TIn, class TOut>
void StressMapping<TIn, TOut>::init()
{
    Inherit::init();

    m_topology = l_inputTopology.get();
    if (!m_topology)
        this->getContext()->get(m_topology);

    if (!m_topology)
    {
        msg_error() << "No input topology found. StressMapping requires a TetrahedronSetTopology.";
        this->d_componentState.setValue(ComponentState::Invalid);
        return;
    }

    const auto& tetras  = m_topology->getTetrahedra();
    const auto  nTetras = tetras.size();
    const auto  nNodes  = this->fromModel->getSize();

    this->toModel->resize(nTetras);

    //using InVecCoord = typename TIn::VecCoord;
    const InVecCoord& p = this->fromModel->read(core::vec_id::read_access::restPosition)->getValue();
    d_initialPoints.setValue(p);

    //copie & adapte de la fonction reinit dans TetrahedronFEMForceField
    helper::ReadAccessor<Data<InVecCoord>> X0 = d_initialPoints;
    elemShapeFun.resize(nTetras);

    for (size_t i = 0; i < nTetras; ++i)
    {
        Mat44 matVert;
        for (Index k = 0; k < 4; k++) {
            Index ix = tetras[i][k];
            matVert(k, 0) = 1.0;
            for (Index l = 1; l < 4; l++)
                matVert(k, l) = X0[ix][l-1];
        }
        const bool canInvert = type::invertMatrix(elemShapeFun[i], matVert);
        assert(canInvert);
        SOFA_UNUSED(canInvert);
    }
    // elemLambda.resize(nTetras);
    // elemMu.resize(nTetras);

    // Build J: (nTetras * NOut) x (nNodes * NIn)
    // Each output row i has 4 entries of 0.25 at columns t[0..3]
    J.compressedMatrix.resize(nTetras * NOut, nNodes * NIn);
    J.compressedMatrix.reserve(nTetras * 4 * NOut);

    for (size_t i = 0; i < nTetras; ++i)
    {
        const auto& t = tetras[i];

        sofa::type::fixed_array<sofa::Index, 4> nodes = { t[0], t[1], t[2], t[3] };
        std::sort(nodes.begin(), nodes.end());

        for (unsigned r = 0; r < NOut; ++r)
        {
            const auto row = NOut * i + r;
            J.compressedMatrix.startVec(row);
            for (auto nodeIdx : nodes)
            {
                J.compressedMatrix.insertBack(row, NIn * nodeIdx + r) = (OutReal)0.25;
            }
        }
    }
    J.compressedMatrix.finalize();
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::apply(const core::MechanicalParams* /*mparams*/, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    helper::WriteOnlyAccessor<Data<VecDeriv>> out = dOut;
    //helper::ReadAccessor<Data<InVecDeriv>> in  = dIn;

    const InVecCoord& X = this->fromModel->read(core::vec_id::read_access::position)->getValue();
    helper::ReadAccessor<Data<InVecCoord> > X0 =  d_initialPoints;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra(); //tetras correspond au indexed elements du code FEMforcefield
    out.resize(tetras.size());

    InVecCoord U;
    U.resize(X.size());
    for (Index i = 0; i < X0.size(); i++)
        U[i] = X[i] - X0[i];

    for (size_t i = 0; i < tetras.size(); ++i)
    {
        type::Vec<6,Real> vStrain;
        Mat33 gradU;
        const auto& it = tetras[i];
        Mat44& shf = elemShapeFun[i];

        /// compute gradU
        for (Index k = 0; k < 3; k++) {
            for (Index l = 0; l < 3; l++)  {
                gradU(k,l) = 0.0;
                for (Index m = 0; m < 4; m++)
                    gradU(k,l) += shf(l+1,m) * U[it[m]][k];
            }
        }

        Mat33 strain = ((Real)0.5)*(gradU + gradU.transposed());

        for (Index i = 0; i < 3; i++)
            vStrain[i] = strain(i,i);
        vStrain[3] = strain(1,2);
        vStrain[4] = strain(0,2);
        vStrain[5] = strain(0,1);

        Real y = d_youngModulus.getValue()[0];
        Real p = d_poissonRatio.getValue()[0];

        Real lambda = (y * p) / ((1 + p) * (1 - 2*p));
        Real mu = y / (2 * (1 + p));

        /// stress
        VoigtTensor s;

        Real traceStrain = 0.0;
        for (Index k = 0; k < 3; k++) {
            traceStrain += vStrain[k];
            s[k] = vStrain[k]*2*mu;
        }

        for (Index k = 3; k < 6; k++)
            s[k] = vStrain[k]*2*mu;

        for (Index k = 0; k < 3; k++)
            s[k] += lambda*traceStrain;
        
        out[i] = s;
    }
    
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJ(const core::MechanicalParams* /*mparams*/, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    helper::WriteOnlyAccessor<Data<VecDeriv>> out = dOut;
    helper::ReadAccessor<Data<InVecDeriv>> in  = dIn;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra(); //tetras correspond au indexed elements du code FEMforcefield
    out.resize(tetras.size());

    for (size_t i = 0; i < tetras.size(); ++i)
    {
        type::Vec<6,Real> vStrain;
        Mat33 gradU;
        const auto& it = tetras[i];
        Mat44& shf = elemShapeFun[i];

        /// compute gradU
        for (Index k = 0; k < 3; k++) {
            for (Index l = 0; l < 3; l++)  {
                gradU(k,l) = 0.0;
                for (Index m = 0; m < 4; m++)
                    gradU(k,l) += shf(l+1,m) * in[it[m]][k];
            }
        }

        Mat33 strain = ((Real)0.5)*(gradU + gradU.transposed());

        for (Index i = 0; i < 3; i++)
            vStrain[i] = strain(i,i);
        vStrain[3] = strain(1,2);
        vStrain[4] = strain(0,2);
        vStrain[5] = strain(0,1);

        Real y = d_youngModulus.getValue()[0];
        Real p = d_poissonRatio.getValue()[0];

        Real lambda = (y * p) / ((1 + p) * (1 - 2*p));
        Real mu = y / (2 * (1 + p));

        /// stress
        VoigtTensor s;

        Real traceStrain = 0.0;
        for (Index k = 0; k < 3; k++) {
            traceStrain += vStrain[k];
            s[k] = vStrain[k]*2*mu;
        }

        for (Index k = 3; k < 6; k++)
            s[k] = vStrain[k]*2*mu;

        for (Index k = 0; k < 3; k++)
            s[k] += lambda*traceStrain;
        
        out[i] = s;
    }
    // Data<InVecDeriv> testOut;
    // InVecDeriv tmp(this->fromModel->getSize());
    // testOut.setValue(tmp);

    // applyJT(nullptr, testOut, dOut);  // dOut = le stress qu'on vient de calculer

    // const InVecDeriv& res = testOut.getValue();
    // for (size_t n = 0; n < res.size(); ++n)
    //     msg_info() << "node " << n << " displacement from applyJT: " << res[n];
}

template<class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJT(const core::MechanicalParams* /*mparams*/, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{

    helper::WriteAccessor<Data<InVecDeriv>> out = dOut;
    helper::ReadAccessor<Data<VecDeriv>> in  = dIn;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    Real y = d_youngModulus.getValue()[0];
    Real p = d_poissonRatio.getValue()[0];

    for (size_t i = 0; i < tetras.size(); ++i)
    {
        const auto& it   = tetras[i];
        const Mat44& shf = elemShapeFun[i];
        const VoigtTensor& eps = in[i];
        Real trEps = eps[0] + eps[1] + eps[2];
        VoigtTensor s; //epsilon

        for (Index k = 0; k < 3; k++)
                s[k] = (eps[k] + (p/(1-2*p))*trEps) * (y/(1+p));
        for (Index k = 3; k < 6; k++)
            s[k] = eps[k] * (y/(1+p));
    
        Mat33 gradU;
        for (Index k = 0; k < 3; k++)
            gradU(k, k) = s[k]; //remplit la diagonale
        gradU(1,2) = gradU(2,1) = s[3];
        gradU(0,2) = gradU(2,0) = s[4];
        gradU(0,1) = gradU(1,0) = s[5];

        for (Index m = 0; m < 4; m++)
            for (Index k = 0; k < 3; k++)
                for (Index l = 0; l < 3; l++)
                    out[it[m]][k] += shf(l+1, m) * gradU(k, l);
        
    }
 //recuperer sigma (the stress)
 //appliquer la formule de la loi de hooke pour trouver epsilon :
 // tr = s[0] + s[1] + s[2] car sigma est un vecteur colonne de 6 elem qui correspondent a une mat symetrique
 // vec a = (1 + poisson)*s
 // b = poisson * tr
 // vec b = (b b b 0 0 0) (mais sous forme colonne)
 // eps = (1/young) * (vec a - vec b)
 //retrouver U depuis epsilon en utilisant la transposee de shf 
 // U = eps * shfT
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJT(const core::ConstraintParams* /*cparams*/, Data<InMatrixDeriv>& dOut, const Data<MatrixDeriv>& dIn)
{
    InMatrixDeriv& out = *dOut.beginEdit();
    const MatrixDeriv& in = dIn.getValue();

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    Real y = d_youngModulus.getValue()[0];
    Real p = d_poissonRatio.getValue()[0];

    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != in.end(); ++rowIt)
    {
        auto o = out.writeLine(rowIt.index());
        for (typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin(); colIt != rowIt.end(); ++colIt)
        {
            // colIt.index() = index du tétraèdre
            // colIt.val() = VoigtTensor (6 composantes de stress)
            const Index tetraIdx = colIt.index();
            //const VoigtTensor& s = colIt.val();
            const VoigtTensor& eps = colIt.val();
            const auto& it  = tetras[tetraIdx];
            const Mat44& shf = elemShapeFun[tetraIdx];

            Real trEps = eps[0] + eps[1] + eps[2];
            VoigtTensor s;
        
            for (Index k = 0; k < 3; k++)
                s[k] = (eps[k] + (p/(1-2*p))*trEps) * (y/(1+p));
            for (Index k = 3; k < 6; k++)
                s[k] = eps[k] * (y/(1+p));

            Mat33 gradU;
            for (Index k = 0; k < 3; k++)
                gradU(k, k) = s[k]; //remplit la diagonale
            gradU(1,2) = gradU(2,1) = s[3];
            gradU(0,2) = gradU(2,0) = s[4];
            gradU(0,1) = gradU(1,0) = s[5];

            for (Index m = 0; m < 4; m++)
                for (Index k = 0; k < 3; k++)
                    for (Index l = 0; l < 3; l++){
                        InDeriv d;
                        d[k] = shf(l+1, m) * gradU(k, l);
                        o.addCol(it[m], d);
                }
        }
    }
    dOut.endEdit();
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::handleTopologyChange()
{
    if (this->toModel && this->fromModel &&
        this->toModel->getSize() != this->fromModel->getSize())
        this->init();
}

template <class TIn, class TOut>
const sofa::linearalgebra::BaseMatrix* StressMapping<TIn, TOut>::getJ()
{
    return &J;
}

template <class TIn, class TOut>
const typename StressMapping<TIn, TOut>::js_type* StressMapping<TIn, TOut>::getJs()
{
    return &Js;
}

} // namespace sofa::component::mapping::linear