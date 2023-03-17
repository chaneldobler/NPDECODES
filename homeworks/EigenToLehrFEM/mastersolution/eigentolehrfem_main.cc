/**
This homework problem consists of reading a simple, gmesh generated, mesh on
the unit square and solving a simple reaction diffusion system using LehrFEM++
*/

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <memory>
#include <utility>


using F_t = std::function<double(const Eigen::Vector2d &)>;

double getArea(const Eigen::MatrixXd &triangle) {
   return std::abs(
       0.5 *
       ((triangle(0, 1) - triangle(0, 0)) * (triangle(1, 2) - triangle(1, 1)) -
        (triangle(0, 2) - triangle(0, 1)) * (triangle(1, 1) - triangle(1, 0))));
}

Eigen::Matrix<double, 2, 3> gradbarycoordinates(const Eigen::MatrixXd &vertices) {
   Eigen::Matrix3d X;
   X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
   X.block<3, 2>(0, 1) = vertices.transpose();
   return X.inverse().block<2, 3>(1, 0);
}

Eigen::Vector3d getElementVector(const Eigen::MatrixXd &V, const F_t &f) {
   Eigen::Vector3d phi_local = Eigen::Vector3d::Zero();
   // Integrate phi against every basis function
   // Intuition: Integrating is like taking the average

   // f(V.col(0))*b_i(0) + f(v.col(1)*b_i(1) + f(V.col(1)*b_i(1)
   for (int i = 0; i < 3; ++i) {
       phi_local(i) = f(V.col(i));
   }
   phi_local *= getArea(V) / 3;
   return phi_local;
}

Eigen::Matrix3d ElementMatrix_Mass_LFE(const Eigen::MatrixXd &V) {
   Eigen::Matrix3d element_matrix;
   element_matrix <<
       2, 1, 1,
       1, 2, 1,
       1, 1, 2;
   element_matrix *= getArea(V)/12.;
   return element_matrix;
}
/* SAM_LISTING_END_1 */

/**
*  @brief Computation of Element Matrix for the Laplacian
*  @param V The vertices of the triangle
*/
Eigen::Matrix3d ElementMatrix_Lapl_LFE(const Eigen::MatrixXd &V) {
   // Argument \texttt{V} same as \texttt{vertices} in \cref{cpp:gradbarycords}.
   // The function returns the $3\times 3$ element matrix as a fixed size
   // \eigen matix.

   Eigen::Matrix<double, 2, 3> X = gradbarycoordinates(V);
   // compute inner products of gradients through matrix multiplication
   return getArea(V) * X.transpose() * X;
}

/**
*  @brief Computation of full Element Matrix
*  @param V The vertices of the triangle
*/
Eigen::Matrix3d getElementMatrix(const Eigen::MatrixXd &V) {
   return ElementMatrix_Lapl_LFE(V) + ElementMatrix_Mass_LFE(V);
}


Eigen::VectorXd solveReactionDiffusion() {
   // Obtain a triangular mesh of the unit square from the collection of
   // test meshes
   lf::mesh::utils::TPTriagMeshBuilder builder(
       std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
   // Set mesh parameters following the Builder pattern
   // Domain is the unit square
   builder.setBottomLeftCorner(Eigen::Vector2d{0., 0.})
       .setTopRightCorner(Eigen::Vector2d{1, 1})
       .setNumXCells(20)
       .setNumYCells(20);
   auto mesh_p = builder.Build();

   // 1. Linear Finite Element Simplicial Lagrangian Space.
   auto fe_space =
       std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
   // TODO: 1.1 Get a const reference to the mesh *from the fe_space variable*.
   const lf::mesh::Mesh &mesh{*(fe_space->Mesh())};

   // TODO: 1.2 Get the dof handler of the finite element space.
   const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

   // TODO: 1.3 Get the dimension of the finite element space.
   const lf::base::size_type N_dofs(dofh.NumDofs());


   // TODO 2.1: Replace the element matrix provider with a lehrfem++ implementation

   // TODO 2.1.1.: Create an ReactionDiffusionElementMatrixProvider
   //  corresponding to the element matrix of the bilinear form
   //  a(u,v) = \int_{\Omega} grad(u)*grad(v) + uv dx
   //  Hint: You'll need to create two MeshFunctionGlobal elements for
   //  constructing this object.

   auto alpha = [](Eigen::Vector2d x){return 1.;};
   auto gamma = [](Eigen::Vector2d x){return 1.;};
   lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};
   lf::mesh::utils::MeshFunctionGlobal mf_gamma{gamma};

   lf::uscalfe::ReactionDiffusionElementMatrixProvider elmat_provider{fe_space, mf_alpha, mf_gamma};

   // TODO 2.2: Replace the Galerkin matrix assembly with a lehrfem++ implementation
   // ---- 2.2 START ----

   // TODO 2.2.1: Replace the triplets with the lehrfem++ assembly COO matrix
   //    std::vector<Eigen::Triplet<double>> triplets;
   lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

   // TODO 2.2.2: Replace the for loop with a lehrfem++ assembly of the Galerkin matrix
   //   And make the COO matrix a sparse eigen matrix.
   //  Hint: Look at the lf::assemble namespace.
   //   You have to use one of its functions.
   //  Test your code once you're done!

   lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider, A);
   Eigen::SparseMatrix<double> A_crs = A.makeSparse();

   //    for (const lf::mesh::Entity* triangle : mesh.Entities(0)) {
   //        // TODO 1.4: Get the corners of the triangle
   //        // Hint: look at the geometry::Corners function.
   ////        Eigen::MatrixXd corners = lf::geometry::Corners(*(triangle->Geometry()));
   ////        Eigen::Matrix3d A_k = getElementMatrix(corners);
   //
   //        // TODO 2.1.2. Comment 1.4 and use your new element matrix provider
   //        //  to get the element matrix. Test your code once you're done !
   //        Eigen::Matrix3d A_k = elmat_provider.Eval(*triangle);
   //
   //        // TODO 1.5: Get the local to global indexing map for this triangle
   //        auto dof_array = dofh.GlobalDofIndices(*triangle);
   //        for (int i = 0; i < 3; i++) {
   //            for (int j = 0; j < 3; j++) {
   //                // TODO 1.6: emplace back to the triplets vector
   //                // The correspondent entry of the local element matrix
   //                // A(global_i, global_j) += A_k(i,j)
   //                triplets.emplace_back(dof_array[i], dof_array[j], A_k(i, j));
   //            }
   //        }
   //    }
   //
   //    Eigen::SparseMatrix<double> A_crs(N_dofs, N_dofs);
   //    A_crs.setFromTriplets(triplets.begin(), triplets.end());
   //    A_crs.makeCompressed();

   // ---- 2.2 END ----

   // Right-hand side vector;
   Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs);

   const double pi = 3.1415926535897;
   auto f = [pi](const Eigen::Vector2d &x) {
       return (1.0 + 8.0 * pi * pi) * std::cos(2.0 * pi * x(0)) *
              std::cos(2.0 * pi * x(1));
   };
   // TODO: 2.3 Replace the load vector assembly with a lehrfem++ implementation
   // TODO: 2.3.1 Create a ScalarLoadElementVectorProvider
   //  corresponding to the element vector of the linear form
   //   l(v) = \int_{\Omega} f(x)v(x) dx
   //  Hint:  You'll need to create a MeshFunctionGlobal element for
   //   constructing this object.

   lf::mesh::utils::MeshFunctionGlobal mf_f{f};
   lf::fe::ScalarLoadElementVectorProvider elvec_provider{fe_space, mf_f};
   // TODO: 2.3.2: Replace the for loop with a lehrfem++ assembly of the Load Vector
   //  Test your code once you're done!
   //  Q: Why is there a difference in the output plot?
   lf::assemble::AssembleVectorLocally(0, dofh, elvec_provider, phi);

   //    for (const lf::mesh::Entity* triangle : mesh.Entities(0)) {
   //        // TODO 1.7: Get the corners of the triangle
   //        Eigen::MatrixXd corners = lf::geometry::Corners(*(triangle->Geometry()));
   //        Eigen::VectorXd phi_k = getElementVector(corners, f);
   //
   //        // TODO 1.8: Get the local to global indexing map for this triangle
   //        auto dof_array = dofh.GlobalDofIndices(*triangle);
   //        for (int i = 0; i < 3; i++) {
   //            // TODO 1.9: phi(global_i) = phi_k(i)
   //            phi[dof_array[i]] = phi_k[i];
   //        }
   //    }

   // Solution vector
   Eigen::VectorXd sol_vec = Eigen::VectorXd::Zero(N_dofs);

   // TODO 1.10: Solve linear system using Eigen's SparseLU solver
   // A_crs * sol_vec = phi
   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
   solver.compute(A_crs);
   sol_vec = solver.solve(phi);
   return sol_vec;
}

Eigen::MatrixXd solveEvolutionReactionDiffusion(Eigen::VectorXd mu_0, double dt, int iterations) {
   // TODO: 2.4 (Bonus):
   //  Solve  M*mu' = (-0.1*A+M)mu,
   //  where -A+M is the galerkin matrix spawned from the bilinear form
   //   a(u,v) = \int_{\Omega} -0.1*grad(u)*grad(v) + uv dx
   //  and M is the galerkin matrix spawned from the bilinear form
   //   m(u,v) = \int_{\Omega} uv dx
   //  Hint: Use the ReactionDiffusionElementMatrixProvider for both
   //   Galerkin matrices.

   // TODO: 2.4.0 Copy code from solveReactionDiffusion until TODO 1.3
   lf::mesh::utils::TPTriagMeshBuilder builder(
       std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
   builder.setBottomLeftCorner(Eigen::Vector2d{0., 0.})
       .setTopRightCorner(Eigen::Vector2d{1, 1})
       .setNumXCells(20)
       .setNumYCells(20);
   auto mesh_p = builder.Build();
   auto fe_space =
       std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
   const lf::mesh::Mesh &mesh{*(fe_space->Mesh())};
   const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
   const lf::base::size_type N_dofs(dofh.NumDofs());

   lf::assemble::COOMatrix<double> M_A(N_dofs, N_dofs);

   // TODO: 2.4.1 Build -0.1*A+M
   auto alpha = [](Eigen::Vector2d x){return -0.1;};
   auto gamma = [](Eigen::Vector2d x){return 1.;};
   lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};
   lf::mesh::utils::MeshFunctionGlobal mf_gamma{gamma};
   lf::uscalfe::ReactionDiffusionElementMatrixProvider M_A_provider{fe_space, mf_alpha, mf_gamma};
   lf::assemble::AssembleMatrixLocally(0, dofh, dofh, M_A_provider, M_A);

   Eigen::SparseMatrix<double> M_A_crs = M_A.makeSparse();

   lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);
   // TODO: 2.4.2 Build M
   auto alpha_0 = [](Eigen::Vector2d x){return 0.;};
   lf::mesh::utils::MeshFunctionGlobal mf_alpha_0{alpha_0};
   lf::uscalfe::ReactionDiffusionElementMatrixProvider M_provider{fe_space, mf_alpha_0, mf_gamma};
   lf::assemble::AssembleMatrixLocally(0, dofh, dofh, M_provider, M);

   Eigen::SparseMatrix<double> M_crs = M.makeSparse();

   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
   solver.compute(M_crs);

   Eigen::VectorXd mu = mu_0;
   Eigen::MatrixXd solution(iterations, mu.size());
   solution.row(0) = mu;

   for(int i = 1; i < iterations; ++i) {
       // TODO: 2.4.3. Compute phi = (-A+M)_crs * mu
       Eigen::VectorXd phi = M_A_crs * mu;

       // TODO: 2.4.4 Solve  M_crs*dmu_dt = (-A+M)_crs * mu,
       Eigen::VectorXd dmu_dt = solver.solve(phi);

       mu += dmu_dt * dt;
       solution.row(i) = mu;
       std::cout << "Iterations " << i << std::endl;
   }

   return solution;
}

int main() {
   bool BONUS = false;

   Eigen::MatrixXd solution;
   if(BONUS) {
       const int iterations = 200;
       const double dt = 0.001;
       Eigen::VectorXd mu_0 = solveReactionDiffusion();
       solution = solveEvolutionReactionDiffusion(mu_0, dt, iterations);
   } else {
       solution = solveReactionDiffusion();
   }
   std::ofstream solution_file(CURRENT_SOURCE_DIR "/solution.txt");
   if (solution_file.is_open()) {
       solution_file << solution << std::endl;
       std::cout << "Solution saved!" << std::endl;
   }
}
