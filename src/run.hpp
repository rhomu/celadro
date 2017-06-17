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

#ifndef RUN_HPP_
#define RUN_HPP_

/** Phase field simulation of cellular monolayers
 *
 * We follow reference doi:10.1038/srep11745. Each cell is represented by a
 * unique phase field and the dynamics is given by
 *
 *    \partial_t \phi + v \partial_x \phi = -\frac 1 2 V
 *
 * where V = \delta F / \delta \phi is computed from a free energy. The
 * velocity v contains both a passive and an acitve contribution. As it stands
 * this is a dry model, but hydrodynamics can be easily added later. Please
 * refer to the original paper for the definition of the coefficients and more
 * detailed discussion.
 *
 * Solution of the evolution equation is performed using a predictor corrector
 * method (with an arbitrary high number of correction steps).
 *
 * -----------------------------------------------------------------------------
 *
 * In order to reduce the amount of computations associated with a high number
 * of different phases we use the fact that the phases are usually well
 * located (blobs that represent cells) and that it is unnecessary to take a
 * phase field into account far outside this domain. Hence we track each cell
 * and update only the region around it.
 *
 * The domains are obtained by computing the center of mass (com) of each cell
 * and by adding a 'margin' of nodes on each side.
 *
 * The walls are implemented as if they were another cell... more to come here.
 */

// =============================================================================
// Declarations

/** Helper function
 *
 * Update the fields in a square domain that is entirely contained inside the
 * domain, i.e. that is not wrapping around the borders.
 * */
template<typename Ret, typename ...Args>
void UpdateSubDomain(Ret (*)(unsigned, unsigned, Args...),
                     unsigned, unsigned, unsigned, unsigned,
                     unsigned, Args&&... args);
/** Parallel version */
template<typename Ret, typename ...Args>
void UpdateSubDomainP(Ret (*)(unsigned, unsigned, Args...),
                      unsigned, unsigned, unsigned, unsigned,
                      unsigned, Args&&... args);
/** Helper function
 *
 * This function is used to updated the fields only in a restricted domain
 * around the cell center. One needs to be careful because of the periodic
 * boundary conditions. The template argument is the function used to update
 * the fields at each node (called ***AtNode()).
 * */
template<typename R, typename ...Args>
void UpdateDomain(R (*)(unsigned, unsigned, Args...),
                  unsigned, Args&&... args);
/** Parallel version */
template<typename R, typename ...Args>
void UpdateDomainP(R (*)(unsigned, unsigned, Args...),
                   unsigned, Args&&... args);

/** Small helper function to add a cell at a certain point */
void AddCell(unsigned, const std::array<unsigned, 2>&);

/** Compute center of mass of a given phase field */
inline void ComputeCoM(unsigned);

/** Initialize all qties that depend on the number of phases
 *
 * We put this in a separate function because it is reused when a cell
 * divides.
 * */
void InitializeFields();

// functions from base class Model
void Initialize();
void Step();
void Configure();
void Pre();
void Post() {};
void PreRunStats();
void RuntimeChecks();
void RuntimeStats();

/** Configure boundary conditions, i.e. the walls */
void ConfigureWalls();

/** Predictor-corrector function for updating the system */
void Update();

/** Subfunction for update */
inline void UpdateAtNode(unsigned, unsigned);
/** Subfunction for update */
inline void UpdateFieldsAtNode(unsigned, unsigned);
/** Subfunction for update */
inline void UpdateStructureTensorAtNode(unsigned, unsigned);
/** Subfunction for update */
inline void UpdateFrictionForceAtNode(unsigned, unsigned);
/** Subfunction for update */
inline void SquareAndSumAtNode(unsigned, unsigned);
/** Subfunction for update */
inline void ReinitSquareAndSumAtNode(unsigned);

/** Make a cell divide
 *
 * The strategy for implementing division is to chose a division axis randomly
 * then divide the given cell in two while fixing all the other cells. We then
 * let the two new cells relax while fixing the other cells such that they are
 * able to create a common interface.
 * */
void Divide(unsigned i);

/** Update polarisation of a given field
 *
 * This function updates the polarisation of the cell which give the direction
 * of the active velocity of the cell. We follow reference 10.1101/095133 and
 * define the dynamics as
 *
 *    d theta / dt = J_r torque + 2 D_r eta
 *
 * where eta is gaussian white noise with zero mean and unit variance, see
 * paper for more details. Note that we use the euler-maruyama method, instead
 * of a predictor-corrector method.
 * */
void UpdatePolarization(unsigned);

/** Update friction force
 *
 * This function updates the friction force and needs to be in a separate
 * loop because it must be computed after the passive part of the velocity has
 * been computed fully. See paper for more details.
 * */
void UpdateFrictionForce();

/** Compute shape parameters
 *
 * This function effectively computes the second moment of area, which ca n be used to
 * fit the shape of a cell to an ellipse.
 * */
void ComputeShape(unsigned);

/** Update the window for tracking */
void UpdateWindow(unsigned);

#if 0
/** Serialization of parameters */
template<class Archive>
void serialize_params(Archive& ar)
{
  ar & auto_name(gamma)
     & auto_name(mu)
     & auto_name(nphases)
     & auto_name(lambda)
     & auto_name(kappa)
     & auto_name(alpha)
     & auto_name(R)
     & auto_name(xi)
     & auto_name(omega)
     & auto_name(init_config)
     & auto_name(zeta)
     & auto_name(D)
     & auto_name(J)
     & auto_name(f)
     & auto_name(f_walls)
     & auto_name(wall_thickness)
     & auto_name(wall_kappa)
     & auto_name(wall_omega)
     & auto_name(walls);

  ar & auto_name(tracking)
     & auto_name(margin);
}

/** Serialization of the current frame */
template<class Archive>
void serialize_frame(Archive& ar)
{
  ar & auto_name(phi)
     & auto_name(area)
     & auto_name(com)
     & auto_name(S_order)
     & auto_name(S_angle)
     & auto_name(pol)
     & auto_name(velp)
     & auto_name(velf)
     & auto_name(velc)
     & auto_name(vel);
  if(tracking) ar
     & auto_name(domain_min)
     & auto_name(domain_max);
}
#endif

#endif//MODELS_PHASES_HPP_
