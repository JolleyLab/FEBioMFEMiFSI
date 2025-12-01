#include "linear_elasticity.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

// void LinearElasticity::AssembleElementVector(const FiniteElement &el,
//                                       ElementTransformation &trans,
//                                       const Vector &elfun,
//                                       Vector &elvect)
// {
//    DenseMatrix elmat;
//    AssembleElementMatrix(el, trans, elmat);
//    elvect.SetSize(elmat.Height());
//    elmat.Mult(elfun,elvect);
// }

void LinearElasticity::EvalPK2Stress(DenseMatrix &stress)
{
   EvalCauchyStress(stress);
}

void LinearElasticity::EvalPK1Stress(DenseMatrix &stress)
{
   EvalCauchyStress(stress);
}

void LinearElasticity::EvalCauchyStress(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   stress.Set(e.Trace()*lambda_val,Identity_Mat);
   stress.Add(2.0 * mu_val,e);
}

// void LinearElasticity::AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
// {
//    int dof = el.GetDof();
//    int dim = el.GetDim();
//    double w, L, M;

//    MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

//    Vector shape;
//    DenseMatrix dshape, gshape, pelmat;
//    Vector divshape;


//    dshape.SetSize(dof, dim);
//    gshape.SetSize(dof, dim);
//    pelmat.SetSize(dof);
//    divshape.SetSize(dim*dof);

//    elmat.SetSize(dof * dim);

//    const IntegrationRule *ir = NULL;
//    if (ir == NULL)
//    {
//       int order = 2 * Trans.OrderGrad(&el); // correct order?
//       ir = &IntRules.Get(el.GetGeomType(), order);
//    }

//    elmat = 0.0;

//    for (int i = 0; i < ir -> GetNPoints(); i++)
//    {
//       const IntegrationPoint &ip = ir->IntPoint(i);

//       el.CalcDShape(ip, dshape);

//       Trans.SetIntPoint(&ip);
//       w = ip.weight * Trans.Weight();
//       Mult(dshape, Trans.InverseJacobian(), gshape);
//       MultAAt(gshape, pelmat);
//       gshape.GradToDiv (divshape);

//       M = mu->Eval(Trans, ip);
//       L = lambda->Eval(Trans, ip);

//       if (L != 0.0)
//       {
//          AddMult_a_VVt(L * w, divshape, elmat);
//       }

//       if (M != 0.0)
//       {
//          for (int d = 0; d < dim; d++)
//          {
//             for (int k = 0; k < dof; k++)
//                for (int l = 0; l < dof; l++)
//                {
//                   elmat (dof*d+k, dof*d+l) += (M * w) * pelmat(k, l);
//                }
//          }
//          for (int ii = 0; ii < dim; ii++)
//             for (int jj = 0; jj < dim; jj++)
//             {
//                for (int kk = 0; kk < dof; kk++)
//                   for (int ll = 0; ll < dof; ll++)
//                   {
//                      elmat(dof*ii+kk, dof*jj+ll) +=
//                         (M * w) * gshape(kk, jj) * gshape(ll, ii);
//                   }
//             }
//       }
//    }
// }   

// void LinearElasticity::AssembleElementGrad(const FiniteElement &el,
//                                      ElementTransformation &trans,
//                                      const Vector &elfun,
//                                      DenseMatrix &elmat)
// {
//    AssembleElementMatrix(el,trans,elmat);   
// }

double LinearElasticity::EvalElasticityTensor(int i, int j, int k, int l)
{
   if (formulation_type == form_type::TotalLagrangian)
   {
      // evaluate material elasticity tensor
      return lambda_val*kron_delta(i,j)*kron_delta(k,l) + mu_val*( kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l)*kron_delta(j,k) );

   }
   else if (formulation_type == form_type::UpdatedLagrangian)
   {
      // evaluate spatial elasticity tensor
      return lambda_val*kron_delta(i,j)*kron_delta(k,l) + mu_val*( kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l)*kron_delta(j,k) );
   }
   else
   {
      MFEM_ABORT("Error: Unsupported Formulation Type");
      return 0.0;
   }
}

double LinearElasticity::ElasticityTensorNorm()
{
   MFEM_ASSERT(formulation_type == form_type::UpdatedLagrangian, "Error this method is only supported for UL approach");

   return lambda_val + 2.0*mu_val;
}

void LinearElasticity::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   lambda_val = lambda->Eval(T,ip);
   mu_val = mu->Eval(T,ip);

   mat_consts.SetSize(2);
   mat_consts(0) = mu_val; // c0
   mat_consts(1) = lambda_val; // c1
}

void LinearElasticity::ResetMaterialConstants()
{
   mu_val = mat_consts(0);
   lambda_val = mat_consts(1);
}
                                 