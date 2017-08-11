void UpdateAtNode(
    double **phi,
    double **phi_old,
    double **V,
    double **potential,
    double **potential_old,
    double *gam,
    double *mu,
    double *area,
    double *area_cnt,
    double *theta,
    vec<double, 2> *pol,
    vec<double, 2> *velp,
    vec<double, 2> *velc,
    vec<double, 2> *velf,
    coord *patch_min,
    coord *patch_max,
    coord *offset,
    vec<double, 2> *com,
    complex<double> *com_x,
    complex<double> *com_y,
    double *sum,
    double *square,
    double *P,
    double *Theta,
    double *Q00,
    double *Q01,
    double *Px,
    double *Py,
    complex<double> *com_x_table,
    complex<double> *com_y_table,
    stencil *neighbors_patch,
    coord patch_size,
    coord patch_margin,
    coord Size,
    double alpha,
    double time_step,
    double xi,
    double C1,
    double C2,
    double J,
    double D,
    bool store
    )
{
  const unsigned n = 0;
  const unsigned q = 0;

  const auto& sq   = neighbors_patch[q];
  const coord qpos = { q/patch_size[1], q%patch_size[1] };
  const coord dpos = ( (qpos + offset[n])%patch_size + patch_min[n] )%Size;
  const auto   k   = dpos[1] + Size[1]*dpos[0];

  // ---------------------------------------------------------------------------
  // Update

  {
    const auto  p  = phi[n][q];
    const auto  a  = area[n];
    const auto  dx = derivX(phi[n], sq);
    const auto  dy = derivY(phi[n], sq);
    const auto  ll = laplacian(phi[n], sq);

    potential[n][q] = (
      // free energy term
      -.5*V[n][q]
      -.5*(
        + C1*gam[n]*p*(1.-p)*(1.-2.*p)
        - 2.*mu[n]*(1.-a/C2)*2.*p
        - 2.*gam[n]*ll
      )
      // advection term
      - (alpha*pol[n][0]+velp[n][0]+velc[n][0]+velf[n][0])*dx/xi
      - (alpha*pol[n][1]+velp[n][1]+velc[n][1]+velf[n][1])*dy/xi
      );
  }

  // store values
  if(store)
  {
    potential_old[n][q] = potential[n][q];
    phi_old[n][q]       = phi[n][q];
  }

  // predictor-corrector
  {
    double p = phi_old[n][q]
               + time_step*.5*(potential[n][q] + potential_old[n][q]);

    // update for next call
    phi[n][q]    = p;
    area_cnt[n] += p*p;
    com_x[n]    += com_x_table[dpos[0]]*p;
    com_y[n]    += com_y_table[dpos[1]]*p;
  }

  // reinit values
  sum[k]    = 0;
  square[k] = 0;
  P[k]      = 0;
  Theta[k]  = 0;
  Q00[k]    = 0;
  Q01[k]    = 0;
  Px[k]     = 0;
  Py[k]     = 0;

  if(q==0)
  {
    // -------------------------------------------------------------------------
    // UpdatePolarization

    if(store)
    {
      // Polarization...
      array<double, 2> v = {
        velp[n][0] + velf[n][0] + velc[n][0],
        velp[n][1] + velf[n][1] + velc[n][1]
      };

      // ...the norm of the passive velocity
      const double ni = sqrt(abs(v[0]*v[0]+v[1]*v[1]));

      // alignement torque
      double torque = -ni*atan2(v[0]*pol[n][1]-v[1]*pol[n][0],
          v[0]*pol[n][0]+v[1]*pol[n][1]);

      // ...euler-marijuana update
      theta[n] += time_step*J*torque + sqrt(time_step)*D*random_normal();

      // update polarisation and contractility
      pol[n] = { cos(theta[n]), sin(theta[n]) };
    }

    // -------------------------------------------------------------------------
    // ComputeCoM

    const auto mx = arg(com_x[n]/static_cast<double>(Size[0]*Size[1])) + Pi;
    const auto my = arg(com_y[n]/static_cast<double>(Size[0]*Size[1])) + Pi;
    com[n] = { mx/2./Pi*Size[0], my/2./Pi*Size[1] };

    // -------------------------------------------------------------------------
    // ComputeShape

    // obtain the new location of the patch min and max
    const coord new_min = {
      (static_cast<unsigned>(round(com[n][0])) + Size[0] - patch_margin[0])%Size[0],
      (static_cast<unsigned>(round(com[n][1])) + Size[1] - patch_margin[1])%Size[1]
    };
    const coord new_max = {
      (static_cast<unsigned>(round(com[n][0])) + patch_margin[0] - 1u)%Size[0],
      (static_cast<unsigned>(round(com[n][1])) + patch_margin[1] - 1u)%Size[1]
    };
    coord displacement = (Size + new_min - patch_min[n])%Size;
    // I guess there is somehthing better than this...
    if(displacement[0]==Size[0]-1u) displacement[0] = patch_size[0]-1u;
    if(displacement[1]==Size[1]-1u) displacement[1] = patch_size[1]-1u;

    // update offset and patch location
    offset[n]    = ( offset[n] + patch_size - displacement )%patch_size;
    patch_min[n] = new_min;
    patch_max[n] = new_max;

    // -------------------------------------------------------------------------
    // UpdatePatch

    com_x[n] = com_y[n] = area[n] = 0.;
  }
}
