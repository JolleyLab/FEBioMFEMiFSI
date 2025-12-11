#pragma once
#include <string>
#include <vector>
#include <cstdio>
#include <memory>
#include <iomanip>
#include <sstream>

#include "../MFEMFSI/IncompressibleFluid/IncompressibleFluidSolver.hpp"
#include <mfem.hpp>

// Base class
class OutputTimeSeriesBase
{
public:
    virtual ~OutputTimeSeriesBase() { if (file) { fclose(file); file = nullptr; } }

    bool Initialize(const std::string &folder_, const std::string &fname_)
    {
        path = folder_ + "/" + fname_;
        #ifdef MFEM_USE_MPI
        bool root = mfem::Mpi::Root();
        if (root)
        {
            file = fopen(path.c_str(), "w");
            if (!file) return false;
            WriteHeader();
        }
        // return true for all ranks so the caller logic proceeds
        return true;
        #else
        file = fopen(path.c_str(), "w");
        if (!file) return false;
        WriteHeader();
        return true;
        #endif
    }

    void Append(double time)
    {
        // All ranks must participate in the collective computations.
        ComputeValues();
        #ifdef MFEM_USE_MPI
        if (!mfem::Mpi::Root()) return;
        #endif
        if (!file) return;
        WriteRow(time); // Only root writes.
        fflush(file);
    }

protected:
    static constexpr int NUM_PREC  = 16;
    static constexpr int NUM_WIDTH = NUM_PREC + 7;
    static constexpr int NUM_WIDTH_FIRST = NUM_WIDTH - 2;
    static constexpr const char* SEP = "    ";

    virtual void BuildHeader(std::vector<std::string> &cols) { cols.push_back("Time"); }
    virtual const char* GetDescriptor() const { return ""; }
    virtual void ComputeValues() = 0;
    virtual void WriteRow(double time) = 0;

    void WriteHeader()
    {
        const char* desc = GetDescriptor();
        if (desc && desc[0] != '\0')
        {
            fputs("# ", file);
            fputs(desc, file);
            fputc('\n', file);
        }
        std::vector<std::string> cols;
        BuildHeader(cols);
        fputs("# ", file);
        if (!cols.empty())
        {
            std::ostringstream line;
            line.setf(std::ios::right);
            line << std::setw(NUM_WIDTH_FIRST-1) << cols[0]; // minus 1 for '#'
            for (size_t i = 1; i < cols.size(); ++i)
            {
                line << SEP << std::setw(NUM_WIDTH) << cols[i];
            }
            line << '\n';
            fputs(line.str().c_str(), file);
        }
    }

protected:
    std::string path;
    FILE *file = nullptr;
};

// Pressure specs
struct SurfacePressureSpec
{
    std::string name;
    int attribute;
    mfem::Array<int> mask;
    double last_area  = 0.0;
    double last_value = 0.0;
};

// Pressure time series
class SurfacePressureTimeSeries : public OutputTimeSeriesBase
{
public:
    SurfacePressureTimeSeries(mfem::Mesh *mesh_,
                              mfem::GridFunction *pressure_gf_,
                              ifs_solver::IncompressibleFluidSolver *solver_)
        : mesh(mesh_), p(pressure_gf_), solver(solver_) {}

    void AddSurface(const std::string &name_, int bdr_attr_)
    {
        MFEM_ASSERT(mesh != nullptr, "SurfacePressureTimeSeries: mesh is null");
        MFEM_ASSERT(mesh->bdr_attributes.Size() > 0,
                    "SurfacePressureTimeSeries: mesh has no boundary attributes"); 
        int maxb = mesh->bdr_attributes.Max();
        SurfacePressureSpec spec;
        spec.name = name_;
        spec.attribute = bdr_attr_;
        spec.mask.SetSize(maxb);
        spec.mask = 0;
        MFEM_ASSERT(bdr_attr_ >= 1 && bdr_attr_ <= maxb, "Invalid boundary attribute");
        spec.mask[bdr_attr_ - 1] = 1;
        surfaces.push_back(spec);
    }

protected:
    const char* GetDescriptor() const override { return "Pressure"; }

    void BuildHeader(std::vector<std::string> &cols) override
    {
        OutputTimeSeriesBase::BuildHeader(cols);
        for (auto &s : surfaces) cols.push_back(s.name);
    }

    void ComputeValues() override
    {
        for (auto &s : surfaces)
        {
            double area = 0.0;
            double val  = solver->ComputeAvgPressureOnSurf(*p, s.mask, area);
            s.last_area  = area;
            s.last_value = val;
        }
    }

    void WriteRow(double time) override
    {
        std::ostringstream line;
        line.setf(std::ios::right);
        line << std::setw(NUM_WIDTH_FIRST)
             << std::scientific << std::setprecision(NUM_PREC) << time;
        for (auto &s : surfaces)
        {
            line << SEP << std::setw(NUM_WIDTH)
                 << std::scientific << std::setprecision(NUM_PREC) << s.last_value;
        }
        line << '\n';
        fputs(line.str().c_str(), file);
    }

private:
    mfem::Mesh *mesh;
    mfem::GridFunction *p;
    ifs_solver::IncompressibleFluidSolver *solver;
    std::vector<SurfacePressureSpec> surfaces;
};

// Flow rate specs
struct SurfaceFlowRateSpec
{
    std::string name;
    int attribute;
    mfem::Array<int> mask;
    double last_value = 0.0;
};

// Flow rate time series
class SurfaceFlowRateTimeSeries : public OutputTimeSeriesBase
{
public:
   SurfaceFlowRateTimeSeries(mfem::Mesh *mesh_,
                           mfem::GridFunction *vel_gf_,
                           ifs_solver::IncompressibleFluidSolver *solver_,
                           double fluid_density_)
      : mesh(mesh_), v(vel_gf_), solver(solver_), rho(fluid_density_) {}

   void AddSurface(const std::string &name_, int bdr_attr_)
   {
      MFEM_ASSERT(mesh != nullptr, "SurfacePressureTimeSeries: mesh is null");
      MFEM_ASSERT(mesh->bdr_attributes.Size() > 0,
                  "SurfacePressureTimeSeries: mesh has no boundary attributes"); 
      int maxb = mesh->bdr_attributes.Max();
      SurfaceFlowRateSpec spec;
      spec.name = name_;
      spec.attribute = bdr_attr_;
      spec.mask.SetSize(maxb);
      spec.mask = 0;
      MFEM_ASSERT(bdr_attr_ >= 1 && bdr_attr_ <= maxb, "Invalid boundary attribute");
      spec.mask[bdr_attr_ - 1] = 1;
      surfaces.push_back(spec);
   }

protected:
    const char* GetDescriptor() const override { return "FlowRate"; }

    void BuildHeader(std::vector<std::string> &cols) override
    {
        cols.push_back("Time");
        for (auto &s : surfaces) cols.push_back(s.name);
    }

    void ComputeValues() override
    {
        for (auto &s : surfaces)
        {
            double mflow = solver->ComputeMassFlowThruSurf(*v, s.mask, false, rho, true);
            double q = mflow / rho;
            s.last_value = q;
        }
    }

    void WriteRow(double time) override
    {
        std::ostringstream line;
        line.setf(std::ios::right);
        line << std::setw(NUM_WIDTH_FIRST)
             << std::scientific << std::setprecision(NUM_PREC) << time;
        for (auto &s : surfaces)
        {
            line << SEP << std::setw(NUM_WIDTH)
                 << std::scientific << std::setprecision(NUM_PREC) << s.last_value;
        }
        line << '\n';
        fputs(line.str().c_str(), file);
    }

private:
    mfem::Mesh *mesh;
    mfem::GridFunction *v;
    ifs_solver::IncompressibleFluidSolver *solver;
    double rho;
    std::vector<SurfaceFlowRateSpec> surfaces;
};