/**
 * @file heatevolutionanimation_main.cc
 * @brief NPDE homework HeatEvolutionAnimation code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include <fstream>
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
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    // TODO 1.2: Get the dofhandler from the fe_space
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
    // Dimension of finite element space`
    const lf::base::size_type N_dofs(dofh.NumDofs());

    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
    lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);


    // Build Stiffness matrix A
    auto one = [](const Eigen::Vector2d&){return 1.;};
    auto zero = [](const Eigen::Vector2d&){return 0.;};
    lf::mesh::utils::MeshFunctionGlobal mf_one{one};
    lf::mesh::utils::MeshFunctionGlobal mf_zero{zero};

    lf::uscalfe::ReactionDiffusionElementMatrixProvider elmat_provider_A{fe_space, mf_one, mf_zero};
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider_A, A);
    Eigen::SparseMatrix<double> A_crs = A.makeSparse();

    // Build Mass Matrix M
    lf::uscalfe::ReactionDiffusionElementMatrixProvider elmat_provider_M{fe_space, mf_zero, mf_one};
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_provider_M, M);
    Eigen::SparseMatrix<double> M_crs = M.makeSparse();

    // Timestepping

    // Initial Conditions
    Eigen::VectorXd mu = Eigen::VectorXd::Zero(N_dofs);
    /*auto ICs = [](const Eigen::Vector2d &x) {
        if(x.norm() < 1./2.) {
            return 1.;
        }
        return 0.;

    };
    lf::mesh::utils::MeshFunctionGlobal mf_ICs{ICs};

    mu = lf::fe::NodalProjection(*fe_space, mf_ICs);*/




    Eigen::MatrixXd solution(steps, N_dofs);
    solution.row(0) = mu;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double> lhs = M_crs + tau*A_crs;

    // Ensure BCs: note: this is not efficient, could be better done if we
    // could have the matrix in RowMajor Format, but Lehrfem++ only allow for
    // ColumnMajor as this is the default of Eigen
    // Boundary Conditions
    auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2);
    auto bd_predicate = [bd_flags](const lf::mesh::Entity& e){return bd_flags(e);};


    for (const lf::mesh::Entity* vertex : mesh.Entities(2)) {
        if(bd_predicate(*vertex)) {
            int global_idx = dofh.GlobalDofIndices(*vertex)[0];
            lhs.row(global_idx) *= 0.;
            lhs.coeffRef(global_idx, global_idx) = 1.;
        }
    }
    solver.compute(lhs);

    Eigen::VectorXd phi, rhs;
    double t = 0.;
    for(int i = 1; i < steps; ++i) {
        t += tau;

        lf::mesh::utils::MeshFunctionGlobal mf_source{source(t)};

        // Evaluate phi on every vertex
        phi = lf::fe::NodalProjection(*fe_space, mf_source);

        rhs = M_crs*mu + tau*phi;

        // Ensure BCs.
        for (const lf::mesh::Entity* vertex : mesh.Entities(2)) {
            if(bd_predicate(*vertex))
                rhs[dofh.GlobalDofIndices(*vertex)[0]] = 0.;
        }

        mu = solver.solve(rhs);
        solution.row(i) = mu;
    }

    std::ofstream solution_file(CURRENT_SOURCE_DIR "/solution.txt");
    solution_file << solution << std::endl;

    return 0;
}
