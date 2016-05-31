double rhoFloorInFluidElement         = 1e-20;
double uFloorInFluidElement           = 1e-20;
double temperatureFloorInFluidElement = 1e-20;
double bSqrFloorInFluidElement        = 1e-20;
double adiabaticIndex                 = 5./3;
double ConductionAlpha                = .25;
double ViscosityAlpha                 = .25;


void setFluidElement(const int oldEval,
                     const array primVars[10],
                     const array gCov[4][4],
                     const array gCon[4][4],
                     const array &g,
                     const array &alpha,
                     array &gammaLorentzFactor,
                     array uCon[4],
                     array uCov[4],
                     array bCon[4],
                     array bCov[4],
                     array &bSqr,
                     array &bNorm,
                     array &q,
                     array &deltaP,
                     array NUp[4],
                     array TUpDown[4],
                    )
{
  array rho         = af::max(primVars[0],rhoFloorInFluidElement);
  array u           = af::max(primVars[1],uFloorInFluidElement);
  array u1          = primVars[2];
  array u2          = primVars[3];
  array u3          = primVars[4];
  array B1          = primVars[5];
  array B2          = primVars[6];
  array B3          = primVars[7];
  array qTilde      = primVars[8];
  array deltaPTilde = primVars[9];

  array pressure    = (adiabaticIndex - 1.)*u;
  array temperature = af::max(pressure/rho,temperatureFloorInFluidElement);

  array soundSpeed  = af::sqrt( adiabaticIndex*pressure
                               /(rho+adiabaticIndex*u)
                              );
  
  gammaLorentzFactor =
    af::sqrt(1 + gCov[1][1] * u1 * u1
               + gCov[2][2] * u2 * u2
               + gCov[3][3] * u3 * u3

             + 2*(  gCov[1][2] * u1 * u2
                  + gCov[1][3] * u1 * u3
                  + gCov[2][3] * u2 * u3
                 )
            );
  if (oldEval)
  {
    gammaLorentzFactor.eval();
  }
  /* Reads:
   * -----
   * gCov[1][1], gCov[2][2], gCov[3][3], gCov[1][2], gCov[1][3], gCov[2][3]: 6
   * u1, u2, u3 : 3
   *
   * Writes:
   * ------
   * gammaLorentzFactor : 1 */

  uCon[0] = gammaLorentzFactor/alpha;
  if (oldEval)
  {
    uCon[0].eval();
  }
  /* Reads:
   * -----
   * gammaLorentzFactor, alpha : 2
   *
   * Writes:
   * ------
   * uCon[0] : 1 */


  uCon[1] = u1 - gammaLorentzFactor*gCon[0][1]*alpha;
  if (oldEval)
  {
    uCon[1].eval();
  }
  /* Reads:
   * -----
   * u1, gammaLorentzFactor, gCon[0][1], alpha : 4
   *
   * Writes:
   * ------
   * uCon[1] : 1 */

  uCon[2] = u2 - gammaLorentzFactor*gCon[0][2]*alpha;
  if (oldEval)
  {
    uCon[2].eval();
  }
  /* Reads:
   * -----
   * u2, gammaLorentzFactor, gCon[0][2], alpha : 4
   *
   * Writes:
   * ------
   * uCon[2] : 1 */

  uCon[3] = u3 - gammaLorentzFactor*gCon[0][3]*alpha;
  if (oldEval)
  {
    uCon[3].eval();
  }
  /* Reads:
   * -----
   * u3, gammaLorentzFactor, gCon[0][3], alpha : 4
   *
   * Writes:
   * ------
   * uCon[3] : 1 */

  for (int mu=0; mu < NDIM; mu++)
  {
    uCov[mu] =  gCov[mu][0] * uCon[0]
              + gCov[mu][1] * uCon[1]
              + gCov[mu][2] * uCon[2]
              + gCov[mu][3] * uCon[3];
    if (oldEval)
    {
      uCov[mu].eval();
    }
  } 
  /* Reads:
   * -----
   * gCov[mu][0], gCov[mu][1], gCov[mu][2], gCov[mu][3]: 16
   * uCon[0], uCon[1], uCon[2], uCon[3]: 4 x 4 = 16
   *
   * Writes:
   * ------
   * uCov[mu] : 4 */

  bCon[0] =  B1*uCov[1] + B2*uCov[2] + B3*uCov[3];
  if (oldEval)
  {
    bCon[0].eval();
  }
  /* Reads:
   * -----
   * B1, B2, B3, uCov[1], uCov[2], uCov[3] : 6
   *
   * Writes:
   * ------
   * bCon[0] : 1 */

  bCon[1] = (B1 + bCon[0] * uCon[1])/uCon[0];
  if (oldEval)
  {
    bCon[1].eval();
  }
  /* Reads:
   * -----
   * B1, bCon[0], uCon[1], uCon[0]: 4
   *
   * Writes:
   * ------
   * bCon[1] : 1 */

  bCon[2] = (B2 + bCon[0] * uCon[2])/uCon[0];
  if (oldEval)
  {
    bCon[2].eval();
  }
  /* Reads:
   * -----
   * B2, bCon[0], uCon[2], uCon[0]: 4
   *
   * Writes:
   * ------
   * bCon[2] : 1 */

  bCon[3] = (B3 + bCon[0] * uCon[3])/uCon[0];
  if (oldEval)
  {
    bCon[3].eval();
  }
  /* Reads:
   * -----
   * B3, bCon[0], uCon[3], uCon[0]: 4
   *
   * Writes:
   * ------
   * bCon[3] : 1 */

  for (int mu=0; mu < NDIM; mu++)
  {
    bCov[mu] =  gCov[mu][0] * bCon[0]
              + gCov[mu][1] * bCon[1]
              + gCov[mu][2] * bCon[2]
              + gCov[mu][3] * bCon[3];
    if (oldEval)
    {
      bCov[mu].eval();
    }
  }
  /* Reads:
   * -----
   * gCov[mu][0], gCov[mu][1], gCov[mu][2], gCov[mu][3]: 16
   * bCon[0], bCon[1], bCon[2], bCon[3]: 4 x 4 = 16
   *
   * Writes:
   * ------
   * bCov[mu] : 4 */

  bSqr =  bCon[0]*bCov[0] + bCon[1]*bCov[1]
        + bCon[2]*bCov[2] + bCon[3]*bCov[3] + bSqrFloorInFluidElement;
  if (oldEval)
  {  
    bSqr.eval();
  }
  /* Reads:
   * -----
   * bCon[0], bCon[1], bCon[2], bCon[3]: 4
   * bCov[0], bCov[1], bCov[2], bCov[3]: 4
   *
   * Writes:
   * ------
   * bSqr : 1 */

  bNorm = af::sqrt(bSqr);
  if (oldEval)
  {
    bNorm.eval();
  }
  /* Reads:
   * -----
   * bSqr: 1
   *
   * Writes:
   * ------
   * bNorm : 1 */

  // Note: this needs to be before setFluidElementParameters
  // because the closure relation uses q, deltaP!
  q = qTilde * temperature * af::sqrt(rho*ConductionAlpha*soundSpeed*soundSpeed);
  if (oldEval)
  {
    q.eval();
  }
  /* Reads:
   * -----
   * qTilde, rho, u : 3
   *
   * Writes:
   * ------
   * q : 1 */

  deltaP = deltaPTilde * af::sqrt(temperature * rho * ViscosityAlpha*soundSpeed*soundSpeed);
  if (oldEval)
  {
    deltaP.eval();
  }
  /* Reads:
   * -----
   * deltaPTilde, rho, u : 3
   *
   * Writes:
   * ------
   * deltaP : 1 */

  for (int mu=0; mu < NDIM; mu++)
  {
    NUp[mu] = rho * uCon[mu];

    for (int nu=0; nu < NDIM; nu++)
    {
      TUpDown[mu][nu] =   (rho + u + pressure + bSqr)*uCon[mu]*uCov[nu]
                        + (pressure + 0.5*bSqr)*DELTA(mu, nu)
                        - bCon[mu] * bCov[nu];

      
      TUpDown[mu][nu] += q/bNorm * (uCon[mu]*bCov[nu] + bCon[mu]*uCov[nu]);

      TUpDown[mu][nu] += (- deltaP)       
                         * (  bCon[mu] * bCov[nu]/bSqr
                              - (1./3.)*(DELTA(mu, nu) + uCon[mu]*uCov[nu])
                           );
      
      if (oldEval)
      {
        TUpDown[mu][nu].eval();
      }
      /* Reads:
       * -----
       * rho, u, bSqr, q, deltaP: 5 x 16 = 80
       * uCon[mu], uCov[nu], bCon[mu], bCov[nu]: 4 x 16 = 64
       *
       * Writes:
       * ------
       * TUpDown[mu][nu] : 16 */
    }
    if (oldEval)
    {
      NUp[mu].eval();
    }
    /* Reads:
     * -----
     * rho : 1 x 4 = 4
     * uCon[mu] : 4
     *
     * Writes:
     * ------
     * NUp[mu] : 4 */
  }
  /* Total reads : 271
   * Total writes: 40 */ 

  numReads  = 271;
  numWrites = 40;
  
  if (oldEval == 0)
  {
    std::vector<array> arraysThatNeedEval{gammaLorentzFactor, 
                                          uCon[0], uCon[1], uCon[2], uCon[3],
                                          uCov[0], uCov[1], uCov[2], uCov[3],
                                          bCon[0], bCon[1], bCon[2], bCon[3],
                                          bCov[0], bCov[1], bCov[2], bCov[3],
                                          bSqr, bNorm, q, deltaP
                                          TUpDown[0][0], TUpDown[0][1], 
                                          TUpDown[0][2], TUpDown[0][3],
                                          TUpDown[1][0], TUpDown[1][1], 
                                          TUpDown[1][2], TUpDown[1][3],
                                          TUpDown[2][0], TUpDown[2][1], 
                                          TUpDown[2][2], TUpDown[2][3],
                                          TUpDown[3][0], TUpDown[3][1], 
                                          TUpDown[3][2], TUpDown[3][3],
                                          NUp[0], NUp[1], NUp[2], NUp[3]
                                         };
    af::eval(arraysThatNeedEval.size(), &arraysThatNeedEval[0]);
  }
}

