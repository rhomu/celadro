/*
 * This file is part of CELADRO, Copyright (C) 2016-17, Romain Mueller
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "header.hpp"
#include "model.hpp"
#include "derivatives.hpp"
#include "tools.hpp"
#include <algorithm>
using namespace std;

void Model::Pre()
{
  // we make the system relax (without activity)
  if(relax_time>0)
  {
    vector<double> save_alpha(nphases, 0); swap(alpha,  save_alpha);
    vector<double> save_beta(nphases,0);  swap(beta,save_beta);
    vector<double> save_zetaS(nphases, 0); swap(zetaS,  save_zetaS);
    vector<double> save_zetaQ(nphases, 0); swap(zetaQ,  save_zetaQ);
    vector<double> save_Dnem(nphases, 0); swap(Dnem,  save_Dnem);
    vector<double> save_Dpol(nphases, 0); swap(Dpol,  save_Dpol);
    vector<double> save_Jnem(nphases, 0); swap(Jnem,  save_Jnem);
    vector<double> save_Jpol(nphases, 0); swap(Jpol,  save_Jpol);
    vector<double> save_Knem(nphases, 0); swap(Knem,  save_Knem);
    vector<double> save_Wnem(nphases, 0); swap(Wnem,  save_Wnem);
    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    for(unsigned i=0; i<relax_time*nsubsteps; ++i)
      for(unsigned i=0; i<=npc; ++i) Update(i==0);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    swap(alpha, save_alpha);
    swap(beta,save_beta);
    swap(zetaS, save_zetaS);
    swap(zetaQ, save_zetaQ);
    swap(Jnem, save_Jnem);
    swap(Jpol, save_Jpol);
    swap(Dnem, save_Dnem);
    swap(Dpol, save_Dpol);
    swap(Knem, save_Knem);
    swap(Wnem, save_Wnem);
  }

  if(BC==5 || BC==7 || BC==8) ConfigureWalls(1);
  if(BC==6) ConfigureWalls(0);
  if(rm_confine_after_relax){
    std::cout<<"\nremove confinement after relaxation"<<std::endl;
    wall_type = "nonadhesive";
    ConfigureWalls(1);
  }
  if (wall_type =="nonadhesive"){
    wall_kappa = wall_omega = 0.;
    std::cout<<"\nwall_omega and wall_kappa are forced to be zero for nonadhesive walls"<<std::endl;
  }
  if (obstacle_type =="nonadhesive"){
    obstacle_kappa = obstacle_omega = 0.;
    std::cout<<"\nobstacle_omega and obstacle_kappa are forced to be zero for nonadhesive obstacles"<<std::endl;
  }
}

void Model::Post()
{}

void Model::PreRunStats()
{
  // packing fraction
  {
    // total cell area
    double packing = 0.;
    for(unsigned n=0; n<nphases; ++n)
      packing += R[n]*R[n];
    packing *= Pi;

    // divide by available area
    double area = 0;
    for(unsigned k=0; k<N; ++k)
      area += 1.-walls[k] - obstacles[k];
    packing /= area;
    cout << "Packing fraction = " << packing << endl;
  }
}

void Model::RuntimeStats()
{
  // TBD
}

void Model::RuntimeChecks()
{
  // check that the area is more or less conserved (20%)
  for(unsigned n=0; n<nphases; ++n)
    if(abs(1.-area[n]/(Pi*R[n]*R[n]))>.2)
      throw warning_msg("area is not conserved.");

  for(unsigned n=0; n<nphases; ++n)
  {
    // check that the cells are not leaking, i.e. that at least 90% of the
    // contributions to the area comes from inside the cell (>1/2).
    double a = 0.;
    // compute area of points outside the cell (<1/2)
    for(const auto v : phi[n]) if(v<.5) a += v*v;
    // check that it is less than 5%
    if(a/area[n]>.90)
      throw warning_msg("your cells are leaking!");

    // check that the phase fields stay between 0 and 1
    for(const auto& p : phi[n])
      if(p<-0.5 or p>1.5)
        throw warning_msg("phase-field is not in [0,1]!");
  }
}

void Model::UpdateSumsAtNode(unsigned n, unsigned q)
{
  const auto k = GetIndexFromPatch(n, q);
  const auto p = phi[n][q];

  sum[k]     += p;
  square[k]  += p*p;
  thirdp[k]  += p*p*p;
  fourthp[k] += p*p*p*p;
  sumA[k]    += p*p*area[n];
  sumS00zeta[k]  += p*S00[n]*zetaS[n];
  sumS01zeta[k]  += p*S01[n]*zetaS[n];
  sumQ00zeta[k]  += p*Q00[n]*zetaQ[n];
  sumQ01zeta[k]  += p*Q01[n]*zetaQ[n];
  sumS00[k]  += p*S00[n];
  sumS01[k]  += p*S01[n];
  sumQ00[k]  += p*Q00[n];
  sumQ01[k]  += p*Q01[n];
  P0[k]      += p*polarization[n][0];
  P1[k]      += p*polarization[n][1];
  U0[k]      += p*com_velocity[n][0];
  U1[k]      += p*com_velocity[n][1];
}

void Model::UpdatePotAtNode(unsigned n, unsigned q)
{
  // supress  = x/sqrt(1 + eps*x*x)
  auto supress = [](double x,double eps){
      return x/sqrt(1 + eps*x*x);
  };

  const auto  k  = GetIndexFromPatch(n, q);
  const auto& s  = neighbors[k];
  const auto& sq = neighbors_patch[q];

  const auto p  = phi[n][q];
  const auto a  = area[n];
  const auto ll = laplacian(phi[n], sq);
  const auto ls = laplacian(sum, s);
  
  const double internal = (
      // CH term
      + gam[n]*(8*p*(1-p)*(1-2*p)/lambda - 2*lambda*ll)
      // area conservation term
      - 4*mu[n]/a0[n]*(1-a/a0[n])*p
    );
  double interactions = (
      // repulsion term
      + 2*kappa/lambda*p*(square[k]-p*p)
      // adhesion term
      - 2*omega*lambda*supress(ls-ll,1.0)
    );
  if (wall_type == "hard"){
      interactions += (    
      // repulsion with walls+ 
      2.0*wall_kappa/lambda*p*walls[k]*walls[k]
      // adhesion with walls
      - 2.0*wall_omega*lambda*supress(walls_laplace[k],1.0)
      );
  }
  if (obstacle_type == "hard"){
      interactions += (    
      // repulsion with obstacles 
      2.0*obstacle_kappa/lambda*p*obstacles[k]*obstacles[k]
      // adhesion with obstacles
      - 2.0*obstacle_omega*lambda*supress(obstacles_laplace[k],1.0)
      );
  }
  // delta F / delta phi_i
  V[n][q] = internal + interactions;

}

void Model::UpdateVelocityFieldAtNode(unsigned n,unsigned q){
  const auto k = GetIndexFromPatch(n,q);
  velocity_field_x[k] += phi[n][q]*vx[n][q];
  velocity_field_y[k] += phi[n][q]*vy[n][q];
}
void Model::UpdateForcesAtNode(unsigned n, unsigned q)
{
  // step function
  auto step = [](double x){
      if (x>0.0){
        return 1.0;
      }
      return 0.0;
  };
  const auto  k  = GetIndexFromPatch(n, q);
  const auto& sq = neighbors_patch[q];
  const auto p  = phi[n][q];
  const auto dx  = derivX(phi[n], sq);
  const auto dy  = derivY(phi[n], sq);

  fp_x[n][q] = V[n][q]*dx;
  fp_y[n][q] = V[n][q]*dy; 
  double fpol_x = 0.0, fpol_y = 0.0;
  double fnem_x = 0.0, fnem_y = 0.0;
  double fshape_x = 0.0, fshape_y = 0.0;
  double dphi_dot_p = dx * polarization[n][0] + dy * polarization[n][1];

  overlap[n] += step(-dphi_dot_p)*p*p*(square[k] + walls[k]*walls[k] + obstacles[k]*obstacles[k] -p*p);
  //polarization density distributed at the front of cells or at whole cells
  if (pol_distribution ==1){  
          //polarization distributed on the front of a cell
          fpol_x = -step(-dphi_dot_p)*dphi_dot_p*polarization[n][0];
          fpol_y = -step(-dphi_dot_p)*dphi_dot_p*polarization[n][1];
  }
  else if(pol_distribution == 0){
        //polarization distributed on the whole cell
      fpol_x = p*polarization[n][0];
      fpol_y = p*polarization[n][1];
  }else{
      throw error_msg("unknown polarization distribution methodology");
  }
  
  // dipolar force density
  //fnem_x = -zetaS[n]*(S00[n]*dx + S01[n]*dy) - zetaQ[n]*(Q00[n]*dx + Q01[n]*dy);
  //fnem_y = -zetaS[n]*(S01[n]*dx - S00[n]*dy) - zetaQ[n]*(Q01[n]*dx - Q00[n]*dy); 
  // dipolar force density - \nabla \cdot (Q) then expand this using product rule
  const auto& sq_k = neighbors[k]; //stencil for global index k

  double fnem_self_x = - zetaQ[n]*(Q00[n]*dx + Q01[n]*dy);
  double fnem_self_y = - zetaQ[n]*(Q01[n]*dx - Q00[n]*dy); 
  double fshape_self_x = -zetaS[n]*(S00[n]*dx + S01[n]*dy);
  double fshape_self_y = -zetaS[n]*(S01[n]*dx - S00[n]*dy);

  
  double fnem_other_x = -(derivX(sumQ00zeta ,sq_k) + derivY(sumQ01zeta,sq_k)) - fnem_self_x; 
  double fshape_other_x = -(derivX(sumS00zeta,sq_k) + derivY(sumS01zeta,sq_k))- fshape_self_x;
 
  double fnem_other_y = -(derivX(sumQ01zeta,sq_k) - derivY(sumQ00zeta,sq_k)) - fnem_self_x; 
  double fshape_other_y = -(derivX(sumS01zeta,sq_k) - derivY(sumS00zeta,sq_k)) - fshape_self_y;

  if (wall_type == "nonadhesive"){
     fpol_x *= (1.0-walls[k]);
     fpol_y *= (1.0-walls[k]);
     fnem_self_x *= (1.0-walls[k]);
     fnem_self_y *= (1.0-walls[k]);
     fshape_self_x *=(1.0-walls[k]);
     fshape_self_y *=(1.0-walls[k]);
  }

  if (obstacle_type == "nonadhesive"){
     fpol_x *= (1.0-obstacles[k]);
     fpol_y *= (1.0-obstacles[k]);
     fnem_self_x *= (1.0-obstacles[k]);
     fnem_self_y *= (1.0-obstacles[k]);
     fshape_self_x *=(1.0-obstacles[k]);
     fshape_self_y *=(1.0-obstacles[k]);
  } 


  if (self_deformation == 1){
    //exclude the self-elongation term
    fnem_x = fnem_self_x + fnem_other_x;
    fnem_y = fnem_self_y + fnem_other_y;
    fshape_x = fshape_self_x + fshape_other_x;
    fshape_y = fshape_self_y + fshape_other_y;
  }
  else if (self_deformation == 0){
    fnem_x = fnem_other_x;
    fnem_y = fnem_other_y;
    fshape_x = fshape_other_x;
    fshape_y = fshape_other_y;
  }

  // weighted averaging field
  fp_field_x[k] += p*V[n][q]*dx;
  fp_field_y[k] += p*V[n][q]*dy; 
  fpol_field_x[k] += p*fpol_x;
  fpol_field_y[k] += p*fpol_y;
  fdipole_field_x[k] += p*(fnem_x + fshape_x);
  fdipole_field_y[k] += p*(fnem_y + fshape_y);
  
  vx[n][q] = (fp_x[n][q] + fpol_x + fnem_x + fshape_x)/xi[n];
  vy[n][q] = (fp_y[n][q] + fpol_y + fnem_y + fshape_y)/xi[n];


  Fint[n]  += { fp_x[n][q] + fnem_x + fshape_x, fp_y[n][q] + fnem_y + fshape_y}; //interal force in fp vanishes after integration
  Fpassive[n] += {fp_x[n][q], fp_y[n][q]};
  Fshape[n]   += {fshape_x,fshape_y };
  Fnem[n]     += {fnem_x,fnem_y};
  Fpol[n]     += {fpol_x,fpol_y};
  double rep = 2.0*(kappa + wall_kappa + obstacle_kappa)/lambda*p*(square[k]-p*p);
  Frep[n]     += {rep*dx,rep*dy};
 
  // store derivatives
  phi_dx[n][q] = dx;
  phi_dy[n][q] = dy;

  // nematic torques
  tau[n]       += sumQ00[k]*Q01[n] - sumQ01[k]*Q00[n];
  vorticity[n] += U0[k]*dy - U1[k]*dx;
}

void Model::UpdatePhaseFieldAtNode(unsigned n, unsigned q, bool store)
{
  const auto k = GetIndexFromPatch(n, q);

  // compute dphi
  dphi[n][q] =
    // free energy
    - V[n][q]
    // advection term
    - vx[n][q]*phi_dx[n][q] - vy[n][q]*phi_dy[n][q];
    ;
  

  // store values
  if(store)
  {
    dphi_old[n][q] = dphi[n][q];
    phi_old[n][q]  = phi[n][q];
  }

  // predictor-corrector
  {
    double p = phi_old[n][q]
               + time_step*.5*(dphi[n][q] + dphi_old[n][q]);

    // update for next call
    phi[n][q]    = p;
    com_x[n] += com_x_table[GetXPosition(k)]*p;
    com_y[n] += com_y_table[GetYPosition(k)]*p;
    area[n]  += p*p;

  }

  // reinit values: we do reinit values here for the simple reason that it is
  // faster than having a supplementary loop afterwards. There is a race
  // condition in principle here but since we are setting evth back to 0 it
  // should be fine. Note that this should be done before the patches are updated
  ReinitSumsAtNode(k);
}

void Model::UpdateNematic(unsigned n, bool store)
{
  // euler-marijuana update
  if(store)
    theta_nem_old[n] = theta_nem[n] + sqrt_time_step*Dnem[n]*random_normal();

  double F00 = 0, F01 = 0;
  switch(align_nematic_to)
  {
    case 0: // center of mass velocity
    {
      const auto ff = com_velocity[n];
      F00 =   ff[0]*ff[0]
            - ff[1]*ff[1];
      F01 = 2*ff[0]*ff[1];
      break;
    }
    case 1: // interaction force
    {
      const auto ff = Fint[n];
      F00 =   ff[0]*ff[0]
            - ff[1]*ff[1];
      F01 = 2*ff[0]*ff[1];
      break;
    }
    case 2: // deformation tensor
    {
      F00 = S00[n];
      F01 = S01[n];
      break;
    }
    case 3: // polarisation
      const auto ff = polarization[n];
      F00 =   ff[0]*ff[0]
            - ff[1]*ff[1];
      F01 = 2*ff[0]*ff[1];
      break;
     
  }
  const auto strength = pow(F01*F01 + F00*F00, 0.25);

  theta_nem[n] = theta_nem_old[n] - time_step*(
      + Knem[n]*tau[n]
      + Jnem[n]*strength*atan2(F00*Q01[n]-F01*Q00[n], F00*Q00[n]+F01*Q01[n]))
      + Wnem[n]*vorticity[n];
  Q00[n] = Snem[n]*cos(2*theta_nem[n]);
  Q01[n] = Snem[n]*sin(2*theta_nem[n]);
  
}

void Model::UpdatePolarization(unsigned n, bool store)
{
  // euler-marijuana update
  if(store){
    theta_pol_old[n] = theta_pol[n] + sqrt_time_step*Dpol[n]*random_normal();
    polarization_old[n] = {polarization[n][0] + sqrt_time_step*Dpol[n]*random_normal(),
                           polarization[n][1] + sqrt_time_step*Dpol[n]*random_normal()};
  }

  vec<double, 2> ff = {0, 0};
  double t = 0.5*atan2(S01[n],S00[n]);
  vec<double, 2> diff = {0.,0.};
  switch(align_polarization_to)
  {
    case 0: //velocity
      ff = com_velocity[n];
      theta_pol[n] = theta_pol_old[n] - time_step*(Jpol[n]*ff.abs()*atan2(ff[0]*polarization[n][1]-ff[1]*polarization[n][0], ff*polarization[n]));
      polarization[n] = {alpha[n]*Spol[n]*cos(theta_pol[n]), alpha[n]*Spol[n]*sin(theta_pol[n]) };
      break;
    case 1: //interaction force
      ff = Fint[n];
      theta_pol[n] = theta_pol_old[n] - time_step*(Jpol[n]*ff.abs()*atan2(ff[0]*polarization[n][1]-ff[1]*polarization[n][0], ff*polarization[n]));
      polarization[n] = { alpha[n]*Spol[n]*cos(theta_pol[n]), alpha[n]*Spol[n]*sin(theta_pol[n]) };
      break;
    case 2: // shape(abs(polarization) fixed)
      ff = {cos(t),sin(t)};
      //make sure the angle between n and p are smaller thant 90 degree
      if (ff[0]*polarization[n][0] + ff[1]*polarization[n][1] < 0.0 ){
        ff = -sqrt(S01[n]*S01[n]+S00[n]*S00[n])*ff;
      }
      else{
        ff = sqrt(S01[n]*S01[n]+S00[n]*S00[n])*ff;
      }
      theta_pol[n] = theta_pol_old[n] - time_step*(Jpol[n]*ff.abs()*atan2(ff[0]*polarization[n][1]-ff[1]*polarization[n][0], ff*polarization[n]));
      polarization[n] = { alpha[n]*Spol[n]*cos(theta_pol[n]), alpha[n]*Spol[n]*sin(theta_pol[n]) };
      break;
    case 3: //shape
      ff = {cos(t),sin(t)};
      //make sure the angle between n and p are smaller thant 90 degree
      if (ff[0]*polarization[n][0] + ff[1]*polarization[n][1] < 0.0 ){
        ff = -sqrt(S01[n]*S01[n]+S00[n]*S00[n])*ff;
      }
      else{
        ff = sqrt(S01[n]*S01[n]+S00[n]*S00[n])*ff;
      }
      diff = - time_step*Jpol[n]*(polarization[n] - beta[n]*ff);
      polarization[n] = polarization_old[n] - time_step*Jpol[n]*(polarization[n] - beta[n]*ff);
      break;
    case 4: //interaction force
      ff = Frep[n];
      theta_pol[n] = theta_pol_old[n] - time_step*(Jpol[n]*ff.abs()*atan2(ff[0]*polarization[n][1]-ff[1]*polarization[n][0], ff*polarization[n]));
      polarization[n] = { alpha[n]*Spol[n]*cos(theta_pol[n]), alpha[n]*Spol[n]*sin(theta_pol[n]) };
      break;
  }
  polarization[n] *= exp(-eta*overlap[n]);
}

void Model::ComputeCoM(unsigned n)
{
  // the strategy to deal with the periodic boundary conditions is to compute
  // all the integrals in Fourier space and come back at the end. This way the
  // periodicity of the domain is automatically taken into account.
  const auto mx = arg(com_x[n]/static_cast<double>(N)) + Pi;
  const auto my = arg(com_y[n]/static_cast<double>(N)) + Pi;
  com[n] = { mx/2./Pi*Size[0], my/2./Pi*Size[1] };
}

void Model::UpdatePatch(unsigned n)
{
  // obtain the new location of the patch min and max
  const coord com_grd { unsigned(round(com[n][0])), unsigned(round(com[n][1])) };
  const coord new_min = ( com_grd + Size - patch_margin ) % Size;
  const coord new_max = ( com_grd + patch_margin - coord {1u, 1u} ) % Size;
  coord displacement  = ( Size + new_min - patch_min[n] ) % Size;

  // I guess there is somehthing better than this...
  if(displacement[0]==Size[0]-1u) displacement[0] = patch_size[0]-1u;
  if(displacement[1]==Size[1]-1u) displacement[1] = patch_size[1]-1u;

  // update offset and patch location
  offset[n]    = ( offset[n] + patch_size - displacement ) % patch_size;
  patch_min[n] = new_min;
  patch_max[n] = new_max;
}

void Model::UpdateStructureTensorAtNode(unsigned n, unsigned q)
{
  const auto  dx = phi_dx[n][q];
  const auto  dy = phi_dy[n][q];
  S00[n] += 0.5*dy*dy - 0.5*dx*dx;
  S01[n] += -dx*dy;
}

void Model::UpdateNeighbourCells(unsigned n,unsigned q){
  unsigned k = GetIndexFromPatch(n, q);
  double threshold = 0.1;
  if ((phi[n][q] > threshold) && (occupation[k] > -1)){
     nbr_cells[n].insert(occupation[k]);
     nbr_cells[occupation[k]].insert(n);
  }
  if ((phi[n][q] >= threshold) &&(occupation[k] < 0)){
    occupation[k] = n;
  }
}

void Model::ReinitSumsAtNode(unsigned k)
{
  sum[k] = 0;
  square[k] = 0;
  thirdp[k] = 0;
  fourthp[k] = 0;
  sumA[k] = 0;
  sumS00[k] = 0;
  sumS01[k] = 0;
  sumQ00[k] = 0;
  sumQ01[k] = 0;
  sumS00zeta[k] = 0;
  sumS01zeta[k] = 0;
  sumQ00zeta[k] = 0;
  sumQ01zeta[k] = 0;
  U0[k] = 0;
  U1[k] = 0;
}

void Model::UpdateStress(unsigned k){
}

#ifndef _CUDA_ENABLED
void Model::Update(bool store, unsigned nstart)
{
  // Compute all global sums
  for(unsigned n=nstart; n<nphases; ++n)
  { 
    // update only patch (in parallel, each node to a different core)
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned q=0; q<patch_N; ++q){
       UpdateSumsAtNode(n, q);    
    }
  }
  
  // Compute stresses
  //
  // We need another loop because the passive force involves a double sum over all
  // the cells.
  for(unsigned n=nstart; n<nphases; ++n)
  {
    // update only patch (in parallel, each node to a different core)
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned q=0; q<patch_N; ++q)
      UpdatePotAtNode(n, q);
  }
  
  //reset fields for calculation
  for(unsigned k = 0;k <N;k++){
      fp_field_x[k] = fp_field_y[k] = 0.0;
      fpol_field_x[k] = fpol_field_y[k] = 0.0;
      fdipole_field_x[k] = fdipole_field_y[k] = 0.0;
      velocity_field_x[k] = velocity_field_y[k] = 0.0;
      occupation[k] = -1;
      stress_xx[k] = stress_yy[k] = stress_xy[k] = 0.0;
  }
  // compute the neighbour of a cell and the rearrangement
  // std::set is not used because it cannot be serialized properly 
  if (neighbour_tracking){
    std::vector<std::set<unsigned> > nbr_old = nbr_cells;
    for (unsigned n = nstart; n<nphases; n++){
        nbr_cells[n].clear();
    }
    for (unsigned n = nstart; n<nphases;n++){
      for (unsigned q = 0 ; q<patch_N;q++){
        UpdateNeighbourCells(n,q);
      }
    }
    /* 
    for (unsigned n = nstart; n < nphases;n++){
      std::sort(nbr_cells.begin(),nbr_cells.end());
      nbr_cells.erase(std::unique(nbr_cells.begin(),nbr_cells.end()),nbr_cells.end());
    }
    */
    for (unsigned n = nstart; n<nphases; n++){
        if (nbr_old[n] == nbr_cells[n]){
          isNbrChanged[n] = false;
        }
        else{
          isNbrChanged[n] = true;
        }
    }
  }
  // Compute induced force (density) and velocity
  for(unsigned n=nstart; n<nphases; ++n)
  {
    Fpol[n] = Fshape[n] = Fpassive[n] = Fnem[n] = Fint[n] = Frep[n] = {0, 0};
    tau[n] = vorticity[n] = overlap[n] =0;
    // update in restricted patch only
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned q=0; q<patch_N; ++q)
      UpdateForcesAtNode(n, q);
    // normalise and compute total forces and vel
    tau[n]     /= lambda;
    com_velocity[n] = (Fpassive[n] + Fnem[n] + Fshape[n] + Fpol[n])/(xi[n]*area[n]);
  }
  // update stress, this step must be after the force density calculation
  for (unsigned k = 0;k<N;k++){
    UpdateStress(k);
  }
  
  for (unsigned n = nstart; n<nphases;n++){
    //update velocity field
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for (unsigned q = 0; q<patch_N;q++){
      UpdateVelocityFieldAtNode(n,q);
    }
  }
  // Predictor-corrector function for updating the phase fields
  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=nstart; n<nphases; ++n)
  {
    com_x[n] = com_y[n] = area[n] = S00[n] = S01[n]  = 0;
    // only update fields in the restricted patch of field n
    for(unsigned q=0; q<patch_N; ++q)
    {
      UpdatePhaseFieldAtNode(n, q, store);
      UpdateStructureTensorAtNode(n, q);
    }
    // update polarisation
    UpdatePolarization(n, store);
    // update Q-tensor
    UpdateNematic(n, store);
    // update center of mass
    ComputeCoM(n);
    // update patch boundaries
    UpdatePatch(n);
  }
}

#endif
