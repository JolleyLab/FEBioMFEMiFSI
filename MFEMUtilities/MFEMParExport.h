#pragma once

#include <XML/XMLWriter.h>
#include <XML/XMLReader.h>
#include <FECore/FECoreBase.h>
#include <FECore/FEParam.h>
#include <FECore/FEModel.h>
#include <FECore/FESolver.h>

#include "../MFEMImSolid/MFEMImmersedElasticDomain.h"

#include <string>
#include <unordered_map>

#include "mfem.hpp"

class MFEMParExport
{
public:
	// constructor
   MFEMParExport();

   // destructor
   virtual ~MFEMParExport();

   virtual void WriteFiles(std::string &feb_filename_, std::string &part_dir_, int np, bool immersed_febio_solid_);

   void SetDomainAttributeMap(std::map<FESolidDomain*, int> &map) { dom_attr_map = &map; };

   void SetFEModel(FEModel* fem){ m_fem = fem; };
   void SetFESolver(FESolver* solver){ m_solver = solver; };
   FEModel* GetFEModel() { return m_fem; }
   FESolver* GetFESolver() { return m_solver; }

private:

   int myid;
   bool immersed_febio_solid;
   std::string feb_filename;
   std::string part_dir;

   FEModel* m_fem;
   FESolver* m_solver;
   XMLWriter xml_writer;

   // reading original feb file
   XMLReader xml_reader;
   bool m_sectionCacheBuilt = false;
   std::unordered_map<std::string, std::string> m_sectionCache;
   void BuildTopLevelSectionCache();
   bool GetCachedSection(const std::string& name, std::string& out) const;

   std::map<FESolidDomain*, int>* dom_attr_map;

   // main file
   void WriteMainFebFile();
   void WriteModuleSection();
   void WriteControlSection();
   void WriteGlobalsSection();
   void WriteMaterialsSection();
   void WriteOutputSection();

   // local files
   void WriteLocalFebFiles();
   void WriteInitialConditionsSection();
   void WriteBoundaryConditionsSection();
   void WriteLoadsSection();
   void WriteLoadDataSection();
   void WriteMeshSection();
   void WriteNodes();
   void WriteElements();
   void WriteMeshDomainsSection();

   // helper functions   
   void WriteSectionComponent(const char* szname, XMLElement &el, FECoreBase &base);
   void ParseTag(XMLTag& tag, XMLElement &el, FECoreBase &base);
   void WriteAttributes(XMLTag& tag, XMLElement &el, FECoreBase &base);
   void WriteParameters(XMLTag& tag, XMLElement &el, FECoreBase &base);   
   void WriteComponent(FECoreBase& base, XMLElement& el);
   void WriteParameterList(FECoreBase& base);
   void WriteParam(FEParam& p);   
   const char* ElementTypeString(int ntype);
   mfem::Mesh* LoadMeshPart(const std::string &prefix);
};