#pragma once

namespace mfem {
	class Mesh;
}

class FEMesh;

mfem::Mesh* CreateMFEMMesh(FEMesh& mesh);

