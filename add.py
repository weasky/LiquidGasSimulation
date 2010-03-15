from fipy.meshes.numMesh.mesh2D import Mesh2D

class Cylinderizer(Mesh2D):
   def __init__(self, mesh):
       Mesh2D.__init__(self,
                       mesh.getVertexCoords(),
                       mesh.faceVertexIDs,
                       mesh.cellFaceIDs)

   def _getFaceAreas(self):
       return Mesh2D._getFaceAreas(self) * self.getFaceCenters()[0]

   def getCellVolumes(self):
       return Mesh2D.getCellVolumes(self) * self.getCellCenters()[0]



