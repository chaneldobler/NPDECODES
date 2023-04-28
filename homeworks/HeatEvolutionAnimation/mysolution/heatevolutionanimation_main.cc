/**
 * @file heatevolutionanimation_main.cc
 * @brief NPDE homework HeatEvolutionAnimation code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
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

#include "heatevolutionanimation.h"

#define N 20
#define T 2.
#define steps 100

#define PI 3.1415926535897

std::shared_ptr<lf::mesh::Mesh> build_mesh() {
    // Obtain a triangular mesh of the unit square from the collection of
    // test meshes
    lf::mesh::utils::TPTriagMeshBuilder builder(
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
    // Set mesh parameters following the Builder pattern
    // Domain is the unit square
    builder.setBottomLeftCorner(Eigen::Vector2d{-1., -1.})
        .setTopRightCorner(Eigen::Vector2d{1, 1})
        .setNumXCells(N)
        .setNumYCells(N);
    return builder.Build();
}

int main() {
    const double tau = T / steps;

    // Source function
    auto source = [](const double t) {
        return [t](const Eigen::Vector2d &x) {
            Eigen::Vector2d c{sin(PI*t), cos(PI*t)};
            return (x - 0.5*c).norm() < 0.5 ? 1. : 0.;
        };
    };

    std::shared_ptr<lf::mesh::Mesh> mesh_p = build_mesh();
    const lf::mesh::Mesh &mesh{*mesh_p};

    // TODO 1. Setup of Linear Finite Element Simplicial Lagrangian Space.
    // TODO 1.1: Create finite element simplicial space from mesh_p
    auto fe_space = /*TODO Fill me*/
    // Obtain local->global index mapping for current finite element space
    // TODO 1.2: Get the dofhandler from the fe_space
    const lf::assemble::DofHandler &dofh{/*TODO Fill me*/};

    // Dimension of finite element space`
    const lf::base::size_type N_dofs(dofh.NumDofs());

    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
    lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);


    // TODO 2. Build lhs matrix
    // TODO 2.1 Build Stiffness matrix A.
    auto one = [](const Eigen::Vector2d&){return 1.;};
    auto zero = [](const Eigen::Vector2d&){return 0.;};

    lf::mesh::utils::MeshFunctionGlobal mf_one{one};
    lf::mesh::utils::MeshFunctionGlobal mf_zero{zero};

    //  Use the ReactionDiffusionElementMatrixProvider and the
    //  AssembleMatrixLocally functions to build the matrix.


    Eigen::SparseMatrix<double> A_crs = A.makeSparse();

    // TODO 2.2 Build Mass Matrix M.

    //  Use the ReactionDiffusionElementMatrixProvider and the
    //  AssembleMatrixLocally functions to build the matrix.


    Eigen::SparseMatrix<double> M_crs = M.makeSparse();

    // TODO 2.3 Create the final lhs matrix
    Eigen::SparseMatrix<double> lhs = /*TODO Fill me*/

    // TODO 4. Boundary Conditions
    // TODO 4.1 Create a CodimMeshDataSet<bool> from the mesh pointer for the
    //  vertices that belong to the boundary.
    //  Hint: flagEntitiesOnBoundary
//        auto bd_flags =
    // TODO 4.2 Create a predicate (lambda that returns a boolean) that takes
    //  a vertex (of type lf::mesh::Entity) as argument and returns if it belongs
    //  to the boundary or not
//        auto bd_predicate =

    // TODO 4.3: Uncomment the following section. You also have to store in
    //  the variable fill_me the correct value.

    // Ensure BCs: note: this is not efficient, could be better done if we
    // could have the matrix in RowMajor Format, but Lehrfem++ only allow for
    // ColumnMajor as this is the default of Eigen
    // Boundary Conditions
    /*
    const double fill_me = 1.;
    for (const lf::mesh::Entity* vertex : mesh.Entities(2)) {
        if(bd_predicate(*vertex)) {
            int global_idx = dofh.GlobalDofIndices(*vertex)[0];
            lhs.row(global_idx) *= 0.;
            lhs.coeffRef(global_idx, global_idx) = fill_me;
        }
    }
    */

    // TODO 2.4 Evaluate the initial conditions.
    // Initial Conditions
    Eigen::VectorXd mu = Eigen::VectorXd::Zero(N_dofs);

    auto ICs = [](const Eigen::Vector2d &x) {
        if(x.norm() < 1./2.) {
            return 1.;
        }
        return 0.;
    };
    lf::mesh::utils::MeshFunctionGlobal mf_ICs{ICs};

    // Evaluate ICs on vertices of the mesh
    mu = /*TODO Fill me*/


    Eigen::MatrixXd solution(steps, N_dofs);
    solution.row(0) = mu;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(lhs);

    // TODO 3. Timestepping
    Eigen::VectorXd rhs, phi;
    double t = 0.;
    for(int i = 1; i < steps; ++i) {
        t += tau;

        lf::mesh::utils::MeshFunctionGlobal mf_source{source(t)};

        // TODO 3.3: Expand mf_source into the vector phi
        // phi = TODO Fill me

        // TODO 3.1 Create the rhs vector
        rhs = /*TODO Fill me*/

        // TODO 4.4: Uncomment and fill the global_idx with the correct
        //  expression.

        //  Ensure BCs by setting every entry of the rhs corresponding to a
        //  boundary vertex to 0.
        /*for (const lf::mesh::Entity* vertex : mesh.Entities(2)) {
            if(bd_predicate(*vertex)) {
                int global_idx = TODO Fill me
                rhs[global_idx] = 0.;
            }
        }*/

        // TODO 3.2 Iterate. After solving this, you should be able to build and
        //  run the code.
        mu = /*TODO Fill me*/

        solution.row(i) = mu;
    }

    std::ofstream solution_file(CURRENT_SOURCE_DIR "/solution.txt");
    solution_file << solution << std::endl;

    return 0;
}

