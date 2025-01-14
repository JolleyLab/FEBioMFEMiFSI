#pragma once

namespace mfem {
	class Mesh;
}

class FEMesh;

mfem::Mesh* CreateMFEMMesh2D(FEMesh& mesh);

mfem::Mesh* CreateMFEMMesh3D(FEMesh& mesh);
