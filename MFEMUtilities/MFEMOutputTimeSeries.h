#pragma once
#include <string>
#include <vector>
#include <cstdio>
#include <memory>

#include "../MFEMFSI/Navier-Stokes/ins_solver.hpp"
#include <mfem.hpp>

// Base class
class OutputTimeSeriesBase {
public:
    virtual ~OutputTimeSeriesBase() { if (file) { fclose(file); file = nullptr; } }

    bool Initialize(const std::string &folder_, const std::string &fname_) {
        path = folder_ + "/" + fname_;
        file = fopen(path.c_str(), "w");
        if (!file) return false;
        WriteHeader();
        return true;
    }

    void Append(double time) {
        if (!file) return;
        WriteRow(time);
        fflush(file);
    }

protected:
    static constexpr int NUM_PREC = 16;
    static constexpr int NUM_WIDTH = NUM_PREC + 6;
    static constexpr int NUM_WIDTH_FIRST = NUM_WIDTH - 2;

    virtual void BuildHeader(std::vector<std::string> &cols) { cols.push_back("Time"); }
    virtual const char* GetDescriptor() const { return ""; }
    virtual void WriteRow(double time) = 0;

    void WriteHeader() {
        const char* desc = GetDescriptor();
        if (desc && desc[0] != '\0') {
            fputs("# ", file);
            fputs(desc, file);
            fputs("\n", file);
        }
        std::vector<std::string> cols;
        BuildHeader(cols);
        fputs("# ", file);
        if (!cols.empty()) {
            fprintf(file, "%*s", NUM_WIDTH_FIRST, cols[0].c_str());
            for (size_t i = 1; i < cols.size(); ++i) {
                fputs("    ", file);
                fprintf(file, "%*s", NUM_WIDTH, cols[i].c_str());
            }
        }
        fputs("\n", file);
    }

protected:
    std::string path;
    FILE *file = nullptr;
};

// Pressure specs
struct SurfacePressureSpec {
    std::string name;
    int attribute;
    mfem::Array<int> mask;
    double last_area  = 0.0;
    double last_value = 0.0;
};

// Pressure time series
class SurfacePressureTimeSeries : public OutputTimeSeriesBase {
public:
    SurfacePressureTimeSeries(mfem::Mesh *mesh_,
                              mfem::GridFunction *pressure_gf_,
                              ins_solver::INS_Solver *solver_)
        : mesh(mesh_), p(pressure_gf_), solver(solver_) {}

    void AddSurface(const std::string &name_, int bdr_attr_) {
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

    void BuildHeader(std::vector<std::string> &cols) override {
        OutputTimeSeriesBase::BuildHeader(cols);
        for (auto &s : surfaces) cols.push_back(s.name);
    }

    void WriteRow(double time) override {
        char buf[128];
        snprintf(buf, sizeof(buf), "%.*e", NUM_PREC, time);
        fputs(buf, file);
        for (auto &s : surfaces) {
            double area = 0.0;
            double val  = solver->ComputeAvgPressureOnSurf(*p, s.mask, area);
            s.last_area  = area;
            s.last_value = val;
            fputs("    ", file);
            snprintf(buf, sizeof(buf), "%.*e", NUM_PREC, val);
            fputs(buf, file);
        }
        fputs("\n", file);
    }

private:
    mfem::Mesh *mesh;
    mfem::GridFunction *p;
    ins_solver::INS_Solver *solver;
    std::vector<SurfacePressureSpec> surfaces;
};

// Flow rate specs
struct SurfaceFlowRateSpec {
    std::string name;
    int attribute;
    mfem::Array<int> mask;
    double last_value = 0.0;
};

// Flow rate time series
class SurfaceFlowRateTimeSeries : public OutputTimeSeriesBase {
public:
    SurfaceFlowRateTimeSeries(mfem::Mesh *mesh_,
                              mfem::GridFunction *vel_gf_,
                              ins_solver::INS_Solver *solver_,
                              double fluid_density_)
        : mesh(mesh_), v(vel_gf_), solver(solver_), rho(fluid_density_) {}

    void AddSurface(const std::string &name_, int bdr_attr_) {
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

    void BuildHeader(std::vector<std::string> &cols) override {
        cols.push_back("Time");
        for (auto &s : surfaces) cols.push_back(s.name);
    }

    void WriteRow(double time) override {
        char buf[128];
        snprintf(buf, sizeof(buf), "%.*e", NUM_PREC, time);
        fputs(buf, file);
        for (auto &s : surfaces) {
            double mflow = solver->ComputeMassFlowThruSurf(*v, s.mask, false, rho, true);
            double q = mflow / rho;
            s.last_value = q;
            fputs("    ", file);
            snprintf(buf, sizeof(buf), "%.*e", NUM_PREC, q);
            fputs(buf, file);
        }
        fputs("\n", file);
    }

private:
    mfem::Mesh *mesh;
    mfem::GridFunction *v;
    ins_solver::INS_Solver *solver;
    double rho;
    std::vector<SurfaceFlowRateSpec> surfaces;
};