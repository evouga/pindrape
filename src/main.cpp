#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "../include/MeshConnectivity.h"
#include "../include/ElasticShell.h"
#include "../include/MidedgeAngleTanFormulation.h"
#include "../include/MidedgeAngleSinFormulation.h"
#include "../include/MidedgeAverageFormulation.h"
#include "../include/StVKMaterial.h"
#include "../include/TensionFieldStVKMaterial.h"
#include "../include/NeoHookeanMaterial.h"
#include "../include/RestState.h"
#include "../include/StVKMaterial.h"
#include "igl/readOBJ.h"
#include <set>
#include <vector>
#include "polyscope/surface_vector_quantity.h"
#include "../optimization/include/NewtonDescent.h"
#include "../include/types.h"
#include "igl/writePLY.h"

double thickness;
double poisson;
double density;

std::vector<int> pinnedVerts = {128, 194};

void lameParameters(double& alpha, double& beta)
{
    double young = 100000;
    alpha = young * poisson / (1.0 - poisson * poisson);
    beta = young / 2.0 / (1.0 + poisson);
}

template <class SFF>
void run_simulation(const LibShell::MeshConnectivity &mesh,
                    const Eigen::MatrixXd &rest_pos, Eigen::MatrixXd &cur_pos,
                    const std::unordered_set<int> *fixed_verts,
                    double thickness, double lame_alpha, double lame_beta,
                    LibShell::HessianProjectType proj_type) {
  // initialize default edge DOFs (edge director angles)
  Eigen::VectorXd init_edge_DOFs;
  SFF::initializeExtraDOFs(init_edge_DOFs, mesh, cur_pos);

  // initialize the rest geometry of the shell
  LibShell::MonolayerRestState rest_state;

  // set uniform thicknesses
  rest_state.thicknesses.resize(mesh.nFaces(), thickness);

  // initialize first fundamental forms to those of input mesh
  LibShell::ElasticShell<SFF>::firstFundamentalForms(mesh, rest_pos,
                                                     rest_state.abars);
                                               
  int nverts = rest_pos.rows();
  int nfaces = mesh.nFaces();
  Eigen::VectorXd vertexWeights(nverts);
  vertexWeights.setZero();
  for(int i=0; i<nfaces; i++)
  {
    double A = 0.5 * std::sqrt(rest_state.abars[i].determinant());
    for(int j=0; j<3; j++)
    {
      vertexWeights[mesh.faceVertex(i, j)] += density * thickness * A / 3.0;
    }
  }

  // initialize second fundamental forms to those of input mesh
  rest_state.bbars.resize(mesh.nFaces());
  for (int i = 0; i < mesh.nFaces(); i++) {
    rest_state.bbars[i].setZero();
  }

  rest_state.lameAlpha.resize(mesh.nFaces(), lame_alpha);
  rest_state.lameBeta.resize(mesh.nFaces(), lame_beta);

  std::shared_ptr<LibShell::MaterialModel<SFF>> mat = std::make_shared<LibShell::NeoHookeanMaterial<SFF>>();

  // projection matrix
  Eigen::SparseMatrix<double> P;
  std::vector<Eigen::Triplet<double>> Pcoeffs;
  int nedges = mesh.nEdges();
  int nedgedofs = SFF::numExtraDOFs;
  // we only allow fixed vertices in the current implementation
  Eigen::VectorXd fixed_dofs(3 * cur_pos.rows());
  fixed_dofs.setZero();
  int nfree = 0;
  for (int i = 0; i < cur_pos.rows(); i++) {
    if (!fixed_verts || !fixed_verts->count(i)) {
      Pcoeffs.push_back({nfree, 3 * i, 1.0});
      Pcoeffs.push_back({nfree + 1, 3 * i + 1, 1.0});
      Pcoeffs.push_back({nfree + 2, 3 * i + 2, 1.0});
      nfree += 3;
    } else {
      fixed_dofs.segment<3>(3 * i) = cur_pos.row(i).transpose();
    }
  }
  for (int i = 0; i < nedges * nedgedofs; i++) {
    Pcoeffs.push_back(
        Eigen::Triplet<double>(nfree, 3 * cur_pos.rows() + i, 1.0));
    nfree++;
  }

  P.resize(nfree, 3 * cur_pos.rows() + nedges * nedgedofs);
  P.setFromTriplets(Pcoeffs.begin(), Pcoeffs.end());

  int totalDOFs = 3 * cur_pos.rows() + nedges * nedgedofs;

  // project the current position
  auto pos_edgedofs_to_variable = [&](const Eigen::MatrixXd &pos,
                                      const Eigen::VectorXd &edge_DOFs) {
    Eigen::VectorXd var(nfree);
    int n = 0;
    for (int i = 0; i < pos.rows(); i++) {
      if (!fixed_verts || !fixed_verts->count(i)) {
        var.segment<3>(n) = pos.row(i).transpose();
        n += 3;
      }
    }
    var.tail(nedges * nedgedofs) = edge_DOFs;
    return var;
  };

  auto variable_to_pos_edgedofs = [&](const Eigen::VectorXd &var) {
    Eigen::MatrixXd pos(cur_pos.rows(), 3);
    int n = 0;
    for (int i = 0; i < cur_pos.rows(); i++) {
      if (!fixed_verts || !fixed_verts->count(i)) {
        pos.row(i) = var.segment<3>(n).transpose();
        n += 3;
      } else {
        pos.row(i) = fixed_dofs.segment<3>(3 * i).transpose();
      }
    }
    Eigen::VectorXd edge_DOFs = var.tail(nedges * nedgedofs);
    return std::pair<Eigen::MatrixXd, Eigen::VectorXd>{pos, edge_DOFs};
  };

  // energy, gradient, and hessian
  auto obj_func = [&](const Eigen::VectorXd &var, Eigen::VectorXd *grad,
                      Eigen::SparseMatrix<double> *hessian, bool psd_proj) {
    Eigen::MatrixXd pos;
    Eigen::VectorXd edge_DOFs;
    std::vector<Eigen::Triplet<double>> hessian_triplets;
    std::tie(pos, edge_DOFs) = variable_to_pos_edgedofs(var);

    double energy = LibShell::ElasticShell<SFF>::elasticEnergy(
        mesh, pos, edge_DOFs, *mat, rest_state, grad,
        hessian ? &hessian_triplets : nullptr, psd_proj ? proj_type : LibShell::HessianProjectType::kNone);
    
    double g = -9.8;
    for(int i=0; i<nverts; i++)
    {
        energy += g * vertexWeights[i] * pos(i, 1);
        if(grad)
            (*grad)[3*i+1] += g * vertexWeights[i];
    }

    if (grad) {
      if (fixed_verts) {
        *grad = P * (*grad);
      }
    }

    if (hessian) {
      hessian->resize(totalDOFs, totalDOFs);
      hessian->setFromTriplets(hessian_triplets.begin(),
                               hessian_triplets.end());
      if (fixed_verts) {
        *hessian = P * (*hessian) * P.transpose();
      }
    }

    return energy;
  };

  auto find_max_step = [&](const Eigen::VectorXd &x,
                           const Eigen::VectorXd &dir) { return 1.0; };

  Eigen::VectorXd x0 = pos_edgedofs_to_variable(cur_pos, init_edge_DOFs);
  OptSolver::TestFuncGradHessian(obj_func, x0);
  
  int num_steps = 10000;
  double grad_tol = 1e-6;
  double f_tol = 0;
  double x_tol = 0;
  bool is_swap = true;

  OptSolver::NewtonSolver(obj_func, find_max_step, x0, num_steps, grad_tol,
                          x_tol, f_tol, proj_type != LibShell::HessianProjectType::kNone, true, is_swap);

  std::tie(cur_pos, init_edge_DOFs) = variable_to_pos_edgedofs(x0);
}


void robustReadMesh(const std::string &filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    std::string fullname = filename;
    if(!igl::readOBJ(fullname, V, F))
    {
        fullname = std::string("../") + filename;
        if(!igl::readOBJ(fullname, V, F))
        {
            fullname = std::string("../") + filename;
            if(!igl::readOBJ(fullname, V, F))
            {
                std::cerr << "Couldn't load: " << filename << std::endl;
                exit(-1);
            }
        }
    }
}

int main(int argc, char* argv[])
{
    // set up material parameters
    thickness = 0.001;
    poisson = 0.3;    
    density = 500.0;

    // load mesh

    Eigen::MatrixXd origV;
    Eigen::MatrixXd curV;
    Eigen::MatrixXi F;

    double lameAlpha, lameBeta;
    lameParameters(lameAlpha, lameBeta);
    
    robustReadMesh("meshes/rest_mesh.obj", origV, F);
    robustReadMesh("meshes/initial_mesh.obj", curV, F);
    
    polyscope::init();
    auto restMesh = polyscope::registerSurfaceMesh("Rest mesh", origV, F);
    restMesh->setEnabled(false);
    auto curMesh = polyscope::registerSurfaceMesh("Start mesh", curV, F);
    curMesh->setEnabled(false);
    
    Eigen::MatrixXd optV = curV;

    LibShell::HessianProjectType proj_type = LibShell::HessianProjectType::kMaxZero;
    LibShell::MeshConnectivity mesh(F);
    
    std::unordered_set<int> pinnedVertsSet;    
    for(auto it : pinnedVerts)
        pinnedVertsSet.insert(it);
    
    run_simulation<LibShell::MidedgeAngleSinFormulation>(
              mesh, origV, optV, &pinnedVertsSet, thickness, lameAlpha, lameBeta,
              proj_type);
    
    polyscope::registerSurfaceMesh("Final mesh", optV, F);
    igl::writePLY("final.ply", optV, F);
   
    polyscope::show();
}
