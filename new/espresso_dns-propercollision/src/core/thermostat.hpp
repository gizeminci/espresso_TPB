/*
 Copyright (C) 2010,2011,2012,2013 The ESPResSo project
 Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 Max-Planck-Institute for Polymer Research, Theory Group
 
 This file is part of ESPResSo.
 
 ESPResSo is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 ESPResSo is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _THERMOSTAT_H
#define _THERMOSTAT_H
/** \file thermostat.hpp
 
 */

#include <cmath>
#include "utils.hpp"
#include "particle_data.hpp"
#include "random.hpp"
#include "global.hpp"
#include "integrate.hpp"
#include "cells.hpp"
#include "lb.hpp"
#include "dpd.hpp"
#include "virtual_sites.hpp"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


/** \name Thermostat switches*/
/************************************************************/
/*@{*/

#define THERMO_OFF        0
#define THERMO_LANGEVIN   1
#define THERMO_DPD        2
#define THERMO_NPT_ISO    4
#define THERMO_LB         8
#define THERMO_INTER_DPD  16
#define THERMO_GHMC       32

/*@}*/

#if (!defined(FLATNOISE) && !defined(GAUSSRANDOMCUT) && !defined(GAUSSRANDOM))
#define FLATNOISE
#endif

/************************************************
 * exported variables
 ************************************************/

/** Switch determining which thermostat to use. This is a or'd value
 of the different possible thermostats (defines: \ref THERMO_OFF,
 \ref THERMO_LANGEVIN, \ref THERMO_DPD \ref THERMO_NPT_ISO). If it
 is zero all thermostats are switched off and the temperature is
 set to zero.  */
extern int thermo_switch;

/** temperature. */
extern double temperature;

/** Langevin friction coefficient gamma. */
extern double langevin_gamma;

/** Friction coefficient for nptiso-thermostat's inline-function friction_therm0_nptiso */
extern double nptiso_gamma0;
/** Friction coefficient for nptiso-thermostat's inline-function friction_thermV_nptiso */
extern double nptiso_gammav;

/** Number of NVE-MD steps in GHMC Cycle*/
extern int ghmc_nmd;
/** Phi parameter for GHMC partial momenum update step */
extern double ghmc_phi;



/************************************************
 * functions
 ************************************************/


/** initialize constants of the thermostat on
 start of integration */
void thermo_init();

/** very nasty: if we recalculate force when leaving/reentering the integrator,
 a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
 numbers are drawn twice, resulting in a different variance of the random force.
 This is corrected by additional heat when restarting the integrator here.
 Currently only works for the Langevin thermostat, although probably also others
 are affected.
 */
void thermo_heat_up();

/** pendant to \ref thermo_heat_up */
void thermo_cool_down();

#ifdef NPT
/** add velocity-dependend noise and friction for NpT-sims to the particle's velocity
 @param dt_vj  j-component of the velocity scaled by time_step dt
 @return       j-component of the noise added to the velocity, also scaled by dt (contained in prefactors) */
inline double friction_therm0_nptiso(double dt_vj) {
    extern double nptiso_pref1, nptiso_pref2;
    if(thermo_switch & THERMO_NPT_ISO)
#if defined (FLATNOISE)
        return ( nptiso_pref1*dt_vj + nptiso_pref2*(d_random()-0.5) );
#elif defined (GAUSSRANDOMCUT)
    return ( nptiso_pref1*dt_vj + nptiso_pref2*gaussian_random_cut() );
#elif defined (GAUSSRANDOM)
    return ( nptiso_pref1*dt_vj + nptiso_pref2*gaussian_random() );
#else
#error No Noise defined
#endif
    return 0.0;
}

/** add p_diff-dependend noise and friction for NpT-sims to \ref nptiso_struct::p_diff */
inline double friction_thermV_nptiso(double p_diff) {
    extern double nptiso_pref3, nptiso_pref4;
    if(thermo_switch & THERMO_NPT_ISO)
#if defined (FLATNOISE)
        return ( nptiso_pref3*p_diff + nptiso_pref4*(d_random()-0.5) );
#elif defined (GAUSSRANDOMCUT)
    return ( nptiso_pref3*p_diff + nptiso_pref4*gaussian_random_cut() );
#elif defined (GAUSSRANDOM)
    return ( nptiso_pref3*p_diff + nptiso_pref4*gaussian_random() );
#else
#error No Noise defined
#endif
    return 0.0;
}
#endif

/* Linear interpolation to calculate the fluid velocity on each particle */



inline void fluid_velocity(Particle *p, double vfxyz[3]){
    
    extern double * velu;
    extern double * velv;
    extern double * velw;
    extern int size;
    extern int TOTALSIZE;
    // Velocity values at the corners of the grid surrounding the point
    double u000,u100,u010,u001,u101,u011,u110,u111;
    double v000,v100,v010,v001,v101,v011,v110,v111;
    double w000,w100,w010,w001,w101,w011,w110,w111;
    
    // distance between 2 grid points.
    // 0.0981747704 comes from the DNS_Dataset and it should change with changing box dimension
    double incr = 0.0981747704*box_l[0]/6.283185306;
    //int size = SIZE;
    int i,j,k;

//printf("memory allocation for velocity components\n");
//double * velu = (double *)malloc(size*size*size*sizeof(double));
//printf("VEL_U memory allocation for velocity components\n");

//double * velv = (double *)malloc(size*size*size*sizeof(double));
//printf("VEL_V memory allocation for velocity components\n");

//double * velw = (double *)malloc(size*size*size*sizeof(double));
//printf("VEL_W memory allocation for velocity components\n");

    // Initialize all parameters
    vfxyz[0]=0; vfxyz[1]=0; vfxyz[2]=0;
    u000 = 0; u100 = 0; u010 = 0; u001 = 0; u101 = 0; u011 = 0; u110 = 0; u111 = 0;
    v000 = 0; v100 = 0; v010 = 0; v001 = 0; v101 = 0; v011 = 0; v110 = 0; v111 = 0;
    w000 = 0; w100 = 0; w010 = 0; w001 = 0; w101 = 0; w011 = 0; w110 = 0; w111 = 0;

/*    // Open and Read binary data
    FILE *fp, *fq, *fr;

    fp = fopen("u.dat", "rb");
    fq = fopen("v.dat", "rb");
    fr = fopen("w.dat", "rb");

  //printf("Files Should be Opened\n"); 
 
    if (!fp) {printf("Unable to open u_velocity file!");}
    else if (!fq) {printf("Unable to open v_velocity file!");}
    else if (!fr) {printf("Unable to open w_velocity file!");}

      //printf("we have data now assign it\n"); 

    fread(velu, sizeof(double), TOTALSIZE, fp);
//printf("vel_u variables should be initialized\n");
    fread(velv, sizeof(double), TOTALSIZE, fq);
    fread(velw, sizeof(double), TOTALSIZE, fr);

fclose(fp);
fclose(fq);
fclose(fr);

free(velu);
free(velv);
free(velw);

*/
#define veluu(i,j,k) (velu[size*size*i + size*j + k])
#define velvv(i,j,k) (velv[size*size*i + size*j + k])
#define velww(i,j,k) (velu[size*size*i + size*j + k])





//printf("vel_* variables should be initialized\n");
    /*if (p->r.p[0]>box_l[0] || p->r.p[1]>box_l[0] || p->r.p[2]>box_l[0] || p->r.p[0]<0.1 || p->r.p[1]<0.1 || p->r.p[2]<0.1) {
printf("particle is out of boundary\n");
}*/
/*
 int a;
 double new_pos[3];
 double tmp;
 for (a=0;a<3;a++)
 {
  new_pos[a] = p->r.p[a]; 

  if (p->r.p[a]>box_l[a])
      new_pos[a] = p->r.p[a]-box_l[a]; 
  else if (p->r.p[a]<0)
      new_pos[a] = p->r.p[a]+box_l[a];

  p->r.p[a] = new_pos[a]; 
 }
*/
    //printf("particle:%d position: %e, %e, %e, NOW CALCULATE GRID POINTS \n",p->p.identity, p->r.p[0], p->r.p[1], p->r.p[2]);

    // 8 vertices of the grid that surrounds the point
    int x0 = p->r.p[0]/incr;
    int x1 = x0+1;
    
    int y0 = p->r.p[1]/incr;
    int y1 = y0+1;
    
    int z0 = p->r.p[2]/incr;
    int z1 = z0+1;

    //printf("position places, %d, %d, %d, %d, %d, %d\n", x0, x1, y0, y1, z0, z1);
    //printf("%d:particle:%d position: %e, %e, %e position places: %d, %d, %d, %d, %d, %d\n",this_node, p->p.identity, p->r.p[0], p->r.p[1], p->r.p[2], x0, x1, y0, y1, z0, z1);

    // Find the weights for each dimension
    double wx = (p->r.p[0]-(x0*incr))/incr;
    double wy = (p->r.p[1]-(y0*incr))/incr;
    double wz = (p->r.p[2]-(z0*incr))/incr;

    if (p->r.p[0]>box_l[0]) {
     x0=0;
     x1=1;
    }

    if (p->r.p[1]>box_l[1]) {
     y0=0;
     y1=1;
    }

    if (p->r.p[2]>box_l[2]) {
     z0=0;
     z1=1;
    }

    // Look up the values of the 8 points of the grid
    
    // HERE I,J,K VALUES ARE IN REVERSE ORDER INORDE TO MAINTAIN THE CORRECT CONFIGURATION !!!
    // DOUBLE CHECK WHEN THE DNS DATASET CHANGED
    
    u000 = veluu(z0,x0,y0);
    u100 = veluu(z1,x0,y0);
    u010 = veluu(z0,x1,y0);
    u001 = veluu(z0,x0,y1);
    u101 = veluu(z1,x0,y1);
    u011 = veluu(z0,x1,y1);
    u110 = veluu(z1,x1,y0);
    u111 = veluu(z1,x1,y1);
    
    v000 = velvv(z0,x0,y0);
    v100 = velvv(z1,x0,y0);
    v010 = velvv(z0,x1,y0);
    v001 = velvv(z0,x0,y1);
    v101 = velvv(z1,x0,y1);
    v011 = velvv(z0,x1,y1);
    v110 = velvv(z1,x1,y0);
    v111 = velvv(z1,x1,y1);
    
    w000 = velww(z0,x0,y0);
    w100 = velww(z1,x0,y0);
    w010 = velww(z0,x1,y0);
    w001 = velww(z0,x0,y1);
    w101 = velww(z1,x0,y1);
    w011 = velww(z0,x1,y1);
    w110 = velww(z1,x1,y0);
    w111 = velww(z1,x1,y1);

    // Compute the velocity components u, v, w at point
    
    vfxyz[0] = u000*(1-wx)*(1-wy)*(1-wz)+u100*wx*(1-wy)*(1-wz)+u010*(1-wx)*wy*(1-wz)+u001*(1-wx)*(1-wy)*wz+u101*wx*(1-wy)*wz+u011*(1-wx)*wy*wz+u110*wx*wy*(1-wz)+u111*wx*wy*wz;
    vfxyz[1] = v000*(1-wx)*(1-wy)*(1-wz)+v100*wx*(1-wy)*(1-wz)+v010*(1-wx)*wy*(1-wz)+v001*(1-wx)*(1-wy)*wz+v101*wx*(1-wy)*wz+v011*(1-wx)*wy*wz+v110*wx*wy*(1-wz)+v111*wx*wy*wz;
    vfxyz[2] = w000*(1-wx)*(1-wy)*(1-wz)+w100*wx*(1-wy)*(1-wz)+w010*(1-wx)*wy*(1-wz)+w001*(1-wx)*(1-wy)*wz+w101*wx*(1-wy)*wz+w011*(1-wx)*wy*wz+w110*wx*wy*(1-wz)+w111*wx*wy*wz;


    printf("%e %e %e %e %e %e \n",p->r.p[0], p->r.p[1], p->r.p[2], vfxyz[0], vfxyz[1], vfxyz[2]);

}



/** overwrite the forces of a particle with
 the friction term, i.e. \f$ F_i= -\gamma v_i + \xi_i\f$.
 */
inline void friction_thermo_langevin(Particle *p)
{
    extern double langevin_pref1, langevin_pref2, langevin_pref3;
#ifdef LANGEVIN_PER_PARTICLE
    double langevin_pref1_temp, langevin_pref2_temp;
#endif
    
    double vfxyz[3];
    double folded_pos[3];
    int img[3]; 
    
    int j;
#ifdef MASS
    double massf = sqrt(PMASS(*p));
#else
    double massf = 1;
#endif
    
    
#ifdef VIRTUAL_SITES
#ifndef VIRTUAL_SITES_THERMOSTAT
    if (ifParticleIsVirtual(p))
    {
        for (j=0;j<3;j++)
            p->f.f[j]=0;
        return;
    }
#endif
    
#ifdef THERMOSTAT_IGNORE_NON_VIRTUAL
    if (!ifParticleIsVirtual(p))
    {
        for (j=0;j<3;j++)
            p->f.f[j]=0;
        return;
    }
#endif
#endif
    
    //for ( j = 0 ; j < 3 ; j++) {
#ifdef EXTERNAL_FORCES
    if (!(p->l.ext_flag & COORD_FIXED(0)) && !(p->l.ext_flag & COORD_FIXED(1)) && !(p->l.ext_flag & COORD_FIXED(1)))
#endif
    {
#ifdef LANGEVIN_PER_PARTICLE
        
#if defined (FLATNOISE)
        if(p->p.gamma >= 0.) {
            langevin_pref1_temp = -p->p.gamma/time_step;
            
            if(p->p.T >= 0.)
                langevin_pref2_temp = sqrt(24.0*p->p.T*p->p.gamma/time_step);
            else
                langevin_pref2_temp = sqrt(24.0*temperature*p->p.gamma/time_step);
            
            p->f.f[0] = langevin_pref1_temp*(p->m.v[0]-(langevin_pref3*vfxyz[0]))*PMASS(*p)+ langevin_pref2_temp*(d_random()-0.5)*massf;
            p->f.f[1] = langevin_pref1_temp*(p->m.v[1]-(langevin_pref3*vfxyz[1]))*PMASS(*p)+ langevin_pref2_temp*(d_random()-0.5)*massf;
            p->f.f[2] = langevin_pref1_temp*(p->m.v[2]-(langevin_pref3*vfxyz[2]))*PMASS(*p)+ langevin_pref2_temp*(d_random()-0.5)*massf;
        }
        else {
            if(p->p.T >= 0.)
                langevin_pref2_temp = sqrt(24.0*p->p.T*langevin_gamma/time_step);
            else
                langevin_pref2_temp = langevin_pref2;
            
            p->f.f[0] = langevin_pref1*(p->m.v[0]-(langevin_pref3*vfxyz[0]))*PMASS(*p)+ langevin_pref2_temp*(d_random()-0.5)*massf;
            p->f.f[1] = langevin_pref1*(p->m.v[1]-(langevin_pref3*vfxyz[1]))*PMASS(*p)+ langevin_pref2_temp*(d_random()-0.5)*massf;
            p->f.f[2] = langevin_pref1*(p->m.v[2]-(langevin_pref3*vfxyz[2]))*PMASS(*p)+ langevin_pref2_temp*(d_random()-0.5)*massf;
        }
#elif defined (GAUSSRANDOMCUT)
        if(p->p.gamma >= 0.) {
            langevin_pref1_temp = -p->p.gamma/time_step;
            
            if(p->p.T >= 0.)
                langevin_pref2_temp = sqrt(2.0*p->p.T*p->p.gamma/time_step);
            else
                langevin_pref2_temp = sqrt(2.0*temperature*p->p.gamma/time_step);
            
            p->f.f[0] = langevin_pref1_temp*p->m.v[0]*PMASS(*p) + langevin_pref2_temp*gaussian_random_cut()*massf;
            p->f.f[1] = langevin_pref1_temp*p->m.v[1]*PMASS(*p) + langevin_pref2_temp*gaussian_random_cut()*massf;
            p->f.f[2] = langevin_pref1_temp*p->m.v[2]*PMASS(*p) + langevin_pref2_temp*gaussian_random_cut()*massf;
        }
        else {
            if(p->p.T >= 0.)
                langevin_pref2_temp = sqrt(2.0*p->p.T*langevin_gamma/time_step);
            else
                langevin_pref2_temp = langevin_pref2;
            
            p->f.f[0] = langevin_pref1*p->m.v[0]*PMASS(*p) + langevin_pref2_temp*gaussian_random_cut()*massf;
            p->f.f[1] = langevin_pref1*p->m.v[1]*PMASS(*p) + langevin_pref2_temp*gaussian_random_cut()*massf;
            p->f.f[2] = langevin_pref1*p->m.v[2]*PMASS(*p) + langevin_pref2_temp*gaussian_random_cut()*massf;
        }
#elif defined (GAUSSRANDOM)
        if(p->p.gamma >= 0.) {
            langevin_pref1_temp = -p->p.gamma/time_step;
            
            if(p->p.T >= 0.)
                langevin_pref2_temp = sqrt(2.0*p->p.T*p->p.gamma/time_step);
            else
                langevin_pref2_temp = sqrt(2.0*temperature*p->p.gamma/time_step);
            
            p->f.f[0] = langevin_pref1_temp*p->m.v[0]*PMASS(*p) + langevin_pref2_temp*gaussian_random()*massf;
            p->f.f[1] = langevin_pref1_temp*p->m.v[1]*PMASS(*p) + langevin_pref2_temp*gaussian_random()*massf;
            p->f.f[2] = langevin_pref1_temp*p->m.v[2]*PMASS(*p) + langevin_pref2_temp*gaussian_random()*massf;
        }
        else {
            if(p->p.T >= 0.)
                langevin_pref2_temp = sqrt(2.0*p->p.T*langevin_gamma/time_step);
            else
                langevin_pref2_temp = langevin_pref2;
            
            p->f.f[0] = langevin_pref1*p->m.v[0]*PMASS(*p) + langevin_pref2_temp*gaussian_random()*massf;
            p->f.f[1] = langevin_pref1*p->m.v[1]*PMASS(*p) + langevin_pref2_temp*gaussian_random()*massf;
            p->f.f[2] = langevin_pref1*p->m.v[2]*PMASS(*p) + langevin_pref2_temp*gaussian_random()*massf;
        }
#else
#error No Noise defined
#endif
        //}
#else
        
#if defined (FLATNOISE)
        //      printf("Check3\n");
        // only this is executed, not the previous checks,
        // for circular flow and agglomeration
        fluid_velocity(p, vfxyz);
        p->f.f[0] = langevin_pref1*(p->m.v[0]-(langevin_pref3*vfxyz[0]))*PMASS(*p)+langevin_pref2*(d_random()-0.5)*massf;
        p->f.f[1] = langevin_pref1*(p->m.v[1]-(langevin_pref3*vfxyz[1]))*PMASS(*p)+langevin_pref2*(d_random()-0.5)*massf;
        p->f.f[2] = langevin_pref1*(p->m.v[2]-(langevin_pref3*vfxyz[2]))*PMASS(*p)+langevin_pref2*(d_random()-0.5)*massf;
        
#elif defined (GAUSSRANDOMCUT)
        p->f.f[0] = langevin_pref1*p->m.v[0]*PMASS(*p) + langevin_pref2*gaussian_random_cut()*massf;
        p->f.f[1] = langevin_pref1*p->m.v[1]*PMASS(*p) + langevin_pref2*gaussian_random_cut()*massf;
        p->f.f[2] = langevin_pref1*p->m.v[2]*PMASS(*p) + langevin_pref2*gaussian_random_cut()*massf;
#elif defined (GAUSSRANDOM)
        p->f.f[0] = langevin_pref1*p->m.v[0]*PMASS(*p) + langevin_pref2*gaussian_random()*massf;
        p->f.f[1] = langevin_pref1*p->m.v[1]*PMASS(*p) + langevin_pref2*gaussian_random()*massf;
        p->f.f[2] = langevin_pref1*p->m.v[2]*PMASS(*p) + langevin_pref2*gaussian_random()*massf;
#else
#error No Noise defined
#endif
        
#endif
    }
#ifdef EXTERNAL_FORCES
    else {
        p->f.f[0] = 0;
        p->f.f[1] = 0;
        p->f.f[2] = 0;
    }
#endif
    
    
   //printf("%d: Particle: %d pos=(%.3e,%.3e,%.3e), vel=(%.3e,%.3e,%.3e), force=(%.3e,%.3e,%.3e)\n",this_node,p->p.identity,p->r.p[0],p->r.p[1],p->r.p[2],vfxyz[0],vfxyz[1],vfxyz[2],p->f.f[0],p->f.f[1],p->f.f[2]);
    
      memcpy(folded_pos, p->r.p, 3*sizeof(double));
      memcpy(img, p->l.i, 3*sizeof(int));
      fold_position(folded_pos, img);

    
    ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
    THERMO_TRACE(fprintf(stderr,"%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",this_node,p->p.identity,p->f.f[0],p->f.f[1],p->f.f[2]));
}

#ifdef ROTATION
/** set the particle torques to the friction term, i.e. \f$\tau_i=-\gamma w_i + \xi_i\f$.
 The same friction coefficient \f$\gamma\f$ is used as that for translation.
 */
inline void friction_thermo_langevin_rotation(Particle *p)
{
    extern double langevin_pref2;
    
    int j;
#ifdef VIRTUAL_SITES
#ifndef VIRTUAL_SITES_THERMOSTAT
    if (ifParticleIsVirtual(p))
    {
        for (j=0;j<3;j++)
            p->f.torque[j]=0;
        return;
    }
#endif
    
#ifdef THERMOSTAT_IGNORE_NON_VIRTUAL
    if (!ifParticleIsVirtual(p))
    {
        for (j=0;j<3;j++)
            p->f.torque[j]=0;
        return;
    }
#endif
#endif
    for ( j = 0 ; j < 3 ; j++)
    {
#if defined (FLATNOISE)
#ifdef ROTATIONAL_INERTIA
        p->f.torque[j] = -langevin_gamma*p->m.omega[j] *p->p.rinertia[j] + langevin_pref2*sqrt(p->p.rinertia[j]) * (d_random()-0.5);
#else
        p->f.torque[j] = -langevin_gamma*p->m.omega[j] + langevin_pref2*(d_random()-0.5);
#endif
#elif defined (GAUSSRANDOMCUT)
#ifdef ROTATIONAL_INERTIA
        p->f.torque[j] = -langevin_gamma*p->m.omega[j] *p->p.rinertia[j] + langevin_pref2*sqrt(p->p.rinertia[j]) * gaussian_random_cut();
#else
        p->f.torque[j] = -langevin_gamma*p->m.omega[j] + langevin_pref2*gaussian_random_cut();
#endif
#elif defined (GAUSSRANDOM)
#ifdef ROTATIONAL_INERTIA
        p->f.torque[j] = -langevin_gamma*p->m.omega[j] *p->p.rinertia[j] + langevin_pref2*sqrt(p->p.rinertia[j]) * gaussian_random();
#else
        p->f.torque[j] = -langevin_gamma*p->m.omega[j] + langevin_pref2*gaussian_random();
#endif
#else
#error No Noise defined
#endif
    }
    ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: LANG f = (%.3e,%.3e,%.3e)\n",this_node,p->f.f[0],p->f.f[1],p->f.f[2]));
    THERMO_TRACE(fprintf(stderr,"%d: Thermo: P %d: force=(%.3e,%.3e,%.3e)\n",this_node,p->p.identity,p->f.f[0],p->f.f[1],p->f.f[2]));
}
#endif


#endif


