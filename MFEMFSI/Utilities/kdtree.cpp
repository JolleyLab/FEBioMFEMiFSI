// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "kdtree.hpp"

namespace mfem
{

namespace kdtree_v1
{

int KDTreeM::FindClosestPoint(mfem::Vector& pt)
{
   PointS best_candidate;
   best_candidate.sp.xx=pt.GetData();
   //initialize the best candidate
   best_candidate.pos =0;
   best_candidate.dist=Dist(data[0], best_candidate.sp);
   best_candidate.level=0;
   PSearch(data.begin(), data.end(), 0, best_candidate);
   int ind=(best_candidate.sp.xx-gf->GetData())/ndim;
   return ind;
}


void KDTreeM::FindClosestPoint(Vector& pt, int &ind, double& dist)
{
   PointS best_candidate;
   best_candidate.sp.xx=pt.GetData();
   //initialize the best candidate
   best_candidate.pos =0;
   best_candidate.dist=Dist(data[0], best_candidate.sp);
   best_candidate.level=0;
   PSearch(data.begin(), data.end(), 0, best_candidate);

   dist=best_candidate.dist;
   ind=data[best_candidate.pos].ind;
}

void KDTreeM::_FindClosestPoint(Vector& pt,int& ind, double& dist)
{
   PointS best_candidate;
   best_candidate.sp.xx=pt.GetData();


   //initialize the best candidate
   best_candidate.pos =0;
   best_candidate.dist=Dist(data[0], best_candidate.sp);
   double dd;
   for (size_t i=1; i<data.size(); i++)
   {
      dd=Dist(data[i],  best_candidate.sp);
      if (dd<best_candidate.dist)
      {
         best_candidate.pos=i;
         best_candidate.dist=dd;
      }
   }
   dist=best_candidate.dist;
   ind=data[best_candidate.pos].ind;
}

void KDTreeM::PSearch(typename std::vector<NodeND>::iterator itb,
                      typename std::vector<NodeND>::iterator ite,
                      size_t level, PointS& bc)
{
   std::uint8_t dim=(std::uint8_t) (level%ndim);
   size_t siz=ite-itb;
   typename std::vector<NodeND>::iterator mtb=itb+siz/2;

   if (siz>2)
   {
      //median is at itb+siz/2
      level=level+1;
      if ((bc.sp.xx[dim]-bc.dist)>mtb->xx[dim]) //look on the right only
      {
         PSearch(itb+siz/2+1, ite, level, bc);
      }
      else if ((bc.sp.xx[dim]+bc.dist)<mtb->xx[dim]) //look on the left only
      {
         PSearch(itb,itb+siz/2, level, bc);
      }
      else  //check all
      {
         if (bc.sp.xx[dim]<mtb->xx[dim])
         {
            //start with the left portion
            PSearch(itb,itb+siz/2, level, bc);
            //and continue to the right
            if (!((bc.sp.xx[dim]+bc.dist)<mtb->xx[dim]))
            {
               PSearch(itb+siz/2+1, ite, level, bc);
               {
                  //check central one
                  double dd=Dist(*mtb, bc.sp);
                  if (dd<bc.dist) { bc.dist=dd; bc.pos=mtb-data.begin(); bc.level=level; }
               }//end central point check
            }
         }
         else
         {
            //start with the right portion
            PSearch(itb+siz/2+1, ite, level, bc);
            //and continue with left
            if (!((bc.sp.xx[dim]-bc.dist)>mtb->xx[dim]))
            {
               PSearch(itb, itb+siz/2, level, bc);
               {
                  //check central one
                  double dd=Dist(*mtb, bc.sp);
                  if (dd<bc.dist) { bc.dist=dd; bc.pos=mtb-data.begin(); bc.level=level; }
               }//end central point check
            }
         }
      }
   }
   else
   {
      //check the nodes
      double dd;
      for (auto it=itb; it!=ite; it++)
      {
         dd=Dist(*it, bc.sp);
         if (dd<bc.dist) //update bc
         {
            bc.pos=it-data.begin();
            bc.dist=dd;
            bc.level=level;
         }
      }
   }
}

KDTreeClosestPointFinder::KDTreeClosestPointFinder(GridFunction &coords, int dim_, Vector &minbb, Vector &maxbb, double lerr) : dim(dim_)
{
   if (dim==2) { kdt2D=new KDTree2D(); kdt3D=nullptr;}
   if (dim==3) { kdt3D=new KDTree3D(); kdt2D=nullptr;}

   auto coords_tv = coords.GetTrueVector();
   int np = coords_tv.Size()/dim;
   int ind;
   double dist;
   bool flag;

   int ordering = coords.FESpace()->GetOrdering();   

   // build kdtree
   if (dim==2)
   {
      int indx, indx2;
      if (ordering == Ordering::byNODES)
      {
         for (int p=0; p < np; p++)
         {
            if (p == 0)
            {
               indx = p;
               indx2 = np;
            }
            
            kdt2D->AddPoint(coords_tv(indx), coords_tv(indx2), p);

            indx += dim;
            indx2 += dim;
         }
      }
      else if (ordering == Ordering::byVDIM)
      {
         for (int p=0; p < np; p++)
         {
            if (p == 0)
            {
               indx = p;
               indx2 = p+1;
            }

            // add point
            kdt2D->AddPoint(coords_tv(indx), coords_tv(indx2), p);

            indx += dim;
            indx2 += dim;
         }
      }
      else
      {
         MFEM_ABORT("Error: Unsupported ordering type for kdtree search");
      }      
   }
   else if (dim==3)
   {
      int indx, indx2, indx3;
      if (ordering == Ordering::byNODES)
      {
         for (int p=0; p < np; p++)
         {
            if (p == 0)
            {
               indx = p;
               indx2 = np;
               indx3 = 2*np;
            }
            
            kdt3D->AddPoint(coords_tv(indx), coords_tv(indx2), coords_tv(indx3), p);

            indx += dim;
            indx2 += dim;
            indx3 += dim;
         }
      }
      else if (ordering == Ordering::byVDIM)
      {
         for (int p=0; p < np; p++)
         {
            if (p == 0)
            {
               indx = p;
               indx2 = p+1;
               indx3 = p+2;
            }

            // add point
            kdt3D->AddPoint(coords_tv(indx), coords_tv(indx2), coords_tv(indx3), p);

            indx += dim;
            indx2 += dim;
            indx3 += dim;
         }
      }
      else
      {
         MFEM_ABORT("Error: Unsupported ordering type for kdtree search");
      }   
   }

   // build the KDTree
   if (dim==2) { kdt2D->Sort();}
   else {kdt3D->Sort();}
}

KDTreeClosestPointFinder::KDTreeClosestPointFinder(Vector &coords, int ordering, int dim_, Vector &minbb, Vector &maxbb, double lerr) : dim(dim_)
{
   if (dim==2) { kdt2D=new KDTree2D(); kdt3D=nullptr;}
   if (dim==3) { kdt3D=new KDTree3D(); kdt2D=nullptr;}

   int np = coords.Size()/dim;
   int ind;
   double dist;
   bool flag; 

   // build kdtree
   if (dim==2)
   {
      int indx, indx2;
      if (ordering == Ordering::byNODES)
      {
         for (int p=0; p < np; p++)
         {
            if (p == 0)
            {
               indx = p;
               indx2 = np;
            }
            
            kdt2D->AddPoint(coords(indx), coords(indx2), p);

            indx += dim;
            indx2 += dim;
         }
      }
      else if (ordering == Ordering::byVDIM)
      {
         for (int p=0; p < np; p++)
         {
            if (p == 0)
            {
               indx = p;
               indx2 = p+1;
            }

            // add point
            kdt2D->AddPoint(coords(indx), coords(indx2), p);

            indx += dim;
            indx2 += dim;
         }
      }
      else
      {
         MFEM_ABORT("Error: Unsupported ordering type for kdtree search");
      }      
   }
   else if (dim==3)
   {
      int indx, indx2, indx3;
      if (ordering == Ordering::byNODES)
      {
         for (int p=0; p < np; p++)
         {
            if (p == 0)
            {
               indx = p;
               indx2 = np;
               indx3 = 2*np;
            }
            
            kdt3D->AddPoint(coords(indx), coords(indx2), coords(indx3), p);

            indx += dim;
            indx2 += dim;
            indx3 += dim;
         }
      }
      else if (ordering == Ordering::byVDIM)
      {
         for (int p=0; p < np; p++)
         {
            if (p == 0)
            {
               indx = p;
               indx2 = p+1;
               indx3 = p+2;
            }

            // add point
            kdt3D->AddPoint(coords(indx), coords(indx2), coords(indx3), p);

            indx += dim;
            indx2 += dim;
            indx3 += dim;
         }
      }
      else
      {
         MFEM_ABORT("Error: Unsupported ordering type for kdtree search");
      }   
   }

   // build the KDTree
   if (dim==2) { kdt2D->Sort();}
   else {kdt3D->Sort();}
}

void KDTreeClosestPointFinder::FindCPandDist(Vector &pt, Vector &cp, double &dist, int &idx)
{
   int ind;
   double dist_val;
   cp.SetSize(pt.Size());
   cp = 0.0;

   if (dim == 2)
   {
      KDTree2D::PointND pnd;
      KDTree2D::PointND clp;

      pnd.xx[0] = pt(0);
      pnd.xx[1] = pt(1);
      
      // kdt2D->FindClosestPoint(pnd,ind,dist_val);
      kdt2D->FindClosestPoint(pnd, ind, dist_val, clp);

      for (int d = 0; d < dim; d++)
      {
         cp(d) = clp.xx[d];
      }

   }
   else if (dim == 3)
   {
      KDTree3D::PointND pnd;
      KDTree3D::PointND clp;

      pnd.xx[0] = pt(0);
      pnd.xx[1] = pt(1);
      pnd.xx[2] = pt(2);
      
      // kdt3D->FindClosestPoint(pnd,ind,dist_val);
      kdt3D->FindClosestPoint(pnd, ind, dist_val, clp);

      for (int d = 0; d < dim; d++)
      {
         cp(d) = clp.xx[d];
      }
   }
   else
   {
      MFEM_ABORT("Error: KDTree not supported in 1D");
   }
   dist = dist_val;
   idx = ind;
}

}
}