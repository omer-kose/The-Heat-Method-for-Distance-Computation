#include <iostream>
#include <vector>



#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"


#include "Eigen/Dense"
#include "Eigen/Sparse"



#include "imgui.h"



using namespace geometrycentral;
using namespace geometrycentral::surface;


//Mesh and Geometry
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;


//Global Properties for Heat Method
struct HeatProperties
{
	Vector<double> DELTA; //Kronecker Delta.
	Vector<double> SOLUTION; //Geodesic Distances
	polyscope::SurfaceVertexColorQuantity* solnColors; //Black to Red
	polyscope::SurfaceGraphQuantity* isolines;
	double maxPhi;
	double vertexRadius;
	double isolinesRadius;
};


HeatProperties ht;
std::vector<Vector3> sourcePositions;
int selectedVertex;

//Polyscope properties
polyscope::SurfaceGraphQuantity* showVerts = nullptr;
polyscope::SurfaceVertexScalarQuantity* distanceColors = nullptr;

//ImGui properties
static const char* meshes[]
{ 
	"bunny.obj", 
	"Armadillo.ply", 
	"catlowpoly.off", 
	"centaur.off",
	"lowresdragon.obj",
	"lucy.obj",
	"Man.obj",
	"horse.off"
};
static int selectedMesh;
static int prevMesh;


double meanEdgeLength()
{
	double mean = 0.0;
	for (const Edge& e : mesh->edges())
	{
		mean += geometry->edgeLength(e);
	}

	return mean / mesh->nEdges();
}


//HEAT METHOD IMPLEMENTATION

double halfedgeCotan(const Halfedge& he)
{
	if (he.isInterior())
	{
		const Vector3& topVertex = geometry->inputVertexPositions[he.next().tipVertex().getIndex()];
		const Vector3& v1 = geometry->inputVertexPositions[he.tailVertex().getIndex()] - topVertex;
		const Vector3& v2 = geometry->inputVertexPositions[he.tipVertex().getIndex()] - topVertex;

		return dot(v1, v2) / norm(cross(v1, v2)); //cos/sin
	}
	else
	{
		return 0.0;
	}
}

double cotan(const Edge& e)
{
	double sum = 0.0;
	//Use adjacentInteriorHalfedges in case of boundary edges
	for (const Halfedge& he : e.adjacentInteriorHalfedges())
	{
		sum += halfedgeCotan(he);
	}

	return sum;
}


//Cotan Laplacian. Note that this is a 2 form Laplacian not 0-form. This is important while solving
//Symmetric Positive Definite Cotan Laplace Matrix

//Important note: If the mesh boundary normals vanish along the boundary, 
// or does not have a boundary, 
//Discrete Laplacian is the negative of Smooth Laplacian. 
//Important while discretizing.
Eigen::SparseMatrix<double> laplaceMatrix()
{
	//Set the Sparse Matrix from triplets
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	Eigen::SparseMatrix<double> spMat(mesh->nVertices(), mesh->nVertices());

	for (size_t i = 0; i < mesh->nVertices(); ++i)
	{
		double totalWeight = 0.0;
		for (const Halfedge& he : mesh->vertex(i).outgoingHalfedges())
		{
			double w = 0.5 * (halfedgeCotan(he) + halfedgeCotan(he.twin()));
			triplets.emplace_back(i, he.tipVertex().getIndex(), -w);
			totalWeight += w;
		}

		//Laplacian is Positive Semi Definite normally. But with this trick
		//It becomes positive definite (and it is indeed symmetric) with this trick. Hence, while solving we can
		//use fast factorizations such as Cholesky Factorization to solve systems fastly.
		triplets.emplace_back(i, i, (totalWeight + 1e-8));
	}

	spMat.setFromTriplets(triplets.begin(), triplets.end());

	return spMat;
}


Eigen::SparseMatrix<double> massMatrix()
{
	//Set the Sparse Matrix from triplets
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	Eigen::SparseMatrix<double> spMat(mesh->nVertices(), mesh->nVertices());

	for (size_t i = 0; i < mesh->nVertices(); ++i)
	{
		triplets.emplace_back(i, i, geometry->vertexDualArea(mesh->vertex(i)));
	}

	spMat.setFromTriplets(triplets.begin(), triplets.end());
	
	return spMat;
}


//Gradient per face. This is computed via geometric gradient
//u is diffused heat
FaceData<Vector3> computeGradientField(const Vector<double>& u)
{
	FaceData<Vector3> vectorField(*mesh);
	for (const Face& f : mesh->faces())
	{
		Vector3 grad = Vector3::zero();
		const Vector3& normal = geometry->faceNormal(f);
		for (const Halfedge& he : f.adjacentHalfedges())
		{
			const Vector3& edgeVector = geometry->inputVertexPositions[he.tipVertex()] - geometry->inputVertexPositions[he.tailVertex()];
			grad += u[he.next().tipVertex().getIndex()] * (cross(normal, edgeVector));
		}

		grad /= 2.0 * geometry->faceArea(f);
		vectorField[f] = -grad.normalize(); //Stores negative normalized gradients

	}

	return vectorField;

}


//Compute Integrated Divergence
//X is the gradient field
Vector<double> computeDivergence(const FaceData<Vector3>& X)
{
	Vector<double> integratedDiv(mesh->nVertices());

	for (const Vertex& v : mesh->vertices())
	{
		double val = 0.0;
		for (const Halfedge& he : v.outgoingHalfedges())
		{
			double cotOne = halfedgeCotan(he);
			double cotTwo = halfedgeCotan(he.next().next());
			const Vector3& edgeOne = geometry->inputVertexPositions[he.tipVertex()] - geometry->inputVertexPositions[he.tailVertex()];
			const Vector3& edgeTwo = geometry->inputVertexPositions[he.next().tipVertex()] - geometry->inputVertexPositions[he.tailVertex()];
			val += cotOne * dot(edgeOne, X[he.face()]) + cotTwo * dot(edgeTwo, X[he.face()]);
		}

		val *= 0.5;
		integratedDiv[v.getIndex()] = val;
	}

	return integratedDiv;

}



void subtractMinimumDistance(Vector<double>& phi)
{
	double min = std::numeric_limits<double>::infinity();
	for (size_t i = 0; i < phi.size(); ++i)
	{
		min = std::min(min, phi[i]);
	}

	for (size_t i = 0; i < phi.size(); ++i)
	{
		phi[i] -= min; //Center the distances
	}

}


//Couple of notes before solving. As stated in cotan laplacian, discrete
//Laplacian is the negative of the smooth counterpart for meshes with no boundaries
//Since we are going to compute Geodesic Distances for Meshes with no boundaries or boundaries
// whose normals vanishes along the boundary, 
//For all the equations that Laplace appears we should use negative of it.

//Second note is about differential forms. Normally Cotan Laplacian has a
//1/dualArea term, meaning that it applies a hodge star operation to a integrated
//value. Thus it is normally a 0-form (L = M^-1*C). But here mine Laplacian does
//not contain M^-1 part so it is a 2 form. Thus while solving equations which
//contains Laplacians be sure to convert all other terms to 2 forms. Which means
//incorporating Mass Matrix M in suitable places.

//Delta is the Kronecker Delta.
Vector<double> computeGeodesics(const Vector<double>& delta)
{
	//Compute short time step which is given in the paper
	double scale = meanEdgeLength();
	double timeStep = scale * scale;
	
	//First compute the Diffused Heat by solving a Backward Euler equation.
	//A : M + t * L
	
	Eigen::SparseMatrix<double> L = laplaceMatrix();
	Eigen::SparseMatrix<double> M = massMatrix();
	Eigen::SparseMatrix<double> A = M + timeStep * L;
	 
	//Our Laplacian is positive definite thus we can solve it very fast
	//Using Cholesky factorization
	Vector<double> u = solvePositiveDefinite(A, delta);

	//Secondly, compute the integrated divergence which is the right side of the
	//Poisson equation. 
	Vector<double> divergence = computeDivergence(computeGradientField(u));


	//Solve for geodesic distances
	L *= -1.0;

	ht.SOLUTION = solvePositiveDefinite(L, divergence);
	

	subtractMinimumDistance(ht.SOLUTION);

	return ht.SOLUTION; //Returning has no purpose here but yeah

}

void loadMesh()
{
	std::string meshPath = "Models/" + std::string(meshes[selectedMesh]);
	std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshPath);


	auto* psMesh =
		polyscope::registerSurfaceMesh("mesh",
			geometry->inputVertexPositions,
			mesh->getFaceVertexList(), polyscopePermutations(*mesh));



	//Data init
	ht.DELTA = Vector<double>::Zero(mesh->nVertices());
	ht.SOLUTION = Vector<double>::Zero(mesh->nVertices());
	psMesh->setSurfaceColor({ 1.0, 0.45, 0.0 }); //Orange

	polyscope::requestRedraw();
}

void callback()
{
	//Selecting the mesh
	ImGui::Combo("Meshes", &selectedMesh, meshes, IM_ARRAYSIZE(meshes));
	//Load the mesh 
	if (prevMesh != selectedMesh)
	{
		prevMesh = selectedMesh;
		loadMesh();
	}
	

	//Could not find a way to process mouse inputs in Polyscope.
	//Select vertices one by one
	ImGui::InputInt("New Source", &selectedVertex);
	if (ImGui::Button("Set Sources"))
	{
		if (selectedVertex >= 0)
		{
			//Selected vertices are heat sources. Their Kronecker Delta is 1.
			ht.DELTA[selectedVertex] = 1.0;
			//Push the vertex position to visualize
			sourcePositions.push_back(geometry->inputVertexPositions[selectedVertex]);
			//Not needed for visualization of vertices but Polyscope wants it anyway
			std::vector<std::array<size_t, 2>> vertInd;
			showVerts = polyscope::getSurfaceMesh("mesh")->addSurfaceGraphQuantity("Heat Sources", sourcePositions, vertInd);
			showVerts->setEnabled(true);
			showVerts->setRadius(0.004);
			showVerts->setColor({1.0, 0.02, 0.02});

		}
		

	}

	if (ImGui::Button("Solve"))
	{
		if (sourcePositions.size() > 0)
		{
			computeGeodesics(ht.DELTA);
			std::vector<double> scalar(mesh->nVertices());
			for (size_t i = 0; i < mesh->nVertices(); ++i)
			{
				scalar[i] = ht.SOLUTION[i];
			}

			distanceColors = polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("Distances", scalar);
			distanceColors->setEnabled(true);
			distanceColors->setColorMap("Hot");
			distanceColors->setIsolinesEnabled(true);
			distanceColors->setIsolineDarkness(0.1);
			distanceColors->setIsolineWidth(0.01, true);
		}

	}

	if (ImGui::Button("Reset"))
	{
		ht.DELTA.setZero();
		ht.SOLUTION.setZero();
		sourcePositions.clear();
		std::vector<std::array<size_t, 2>> vertInd;
		//Not sure this is a good way to reset.
		polyscope::getSurfaceMesh("mesh")->addSurfaceGraphQuantity("Heat Sources", sourcePositions, vertInd);
		//polyscope::getSurfaceMesh("mesh")->setSurfaceColor({ 1.0, 0.45, 0.0 }); //Orange
		if (distanceColors != nullptr)
		{
			distanceColors->setEnabled(false);
		}
	}

}


int main()
{
	loadMesh();

	polyscope::init();
	polyscope::state::userCallback = callback;

	
	polyscope::loadColorMap("Hot", "HotColormapFunction.png");





	polyscope::show(); // pass control to the gui until the user exits


	return 0;
}