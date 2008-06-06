//----------------------------------------------------------------------------
// photonmap.c
// An example implementation of the photon map data structure
//
// Henrik Wann Jensen - February 2001
// converted to C by Tim Davis
//----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pm.h"


static void balance_segment(
    Photon **pbal,
    Photon **porg,
    int index,
    int start,
    int end );

static void median_split(
    Photon **p,
    int start,
    int end,
    int median,
    int axis );
  

static PhotonMap pmap;
static PhotonMap *pm = &pmap;

static float MAX_R;

/* This is the constructor for the photon map.
 * To create the photon map it is necessary to specify the
 * maximum number of photons that will be stored */

void init_Photon_map( long max_phot )

{
  int i;
  double angle;

  pm -> stored_photons = 0;
  pm -> prev_scale = 1;
  pm -> max_photons = max_phot;

  pm -> photons = (Photon*)malloc( sizeof( Photon ) * ( pm -> max_photons+1 ) );

  if (pm -> photons == NULL) {
    fprintf(stderr,"Out of memory initializing photon map\n");
    exit(-1);
  }

  pm -> bbox_min[0] = pm -> bbox_min[1] = pm -> bbox_min[2] = 1e8f;
  pm -> bbox_max[0] = pm -> bbox_max[1] = pm -> bbox_max[2] = -1e8f;
  
  //----------------------------------------
  // initialize direction conversion tables
  //----------------------------------------

  for (i=0; i<256; i++) {
    angle = (double)(i)*(1.0/256.0)*M_PI;
    pm -> costheta[i] = cos( angle );
    pm -> sintheta[i] = sin( angle );
    pm -> cosphi[i]   = cos( 2.0*angle );
    pm -> sinphi[i]   = sin( 2.0*angle );
  }
}


void delete_Photon_map(void)

{
  free( pm -> photons );
}


/* photon_dir returns the direction of a photon */

void pm_photon_dir( float *dir, Photon *p )

{
  dir[0] = pm -> sintheta[p->theta]*pm -> cosphi[p->phi];
  dir[1] = pm -> sintheta[p->theta]*pm -> sinphi[p->phi];
  dir[2] = pm -> costheta[p->theta];
}


/* irradiance_estimate computes an irradiance estimate
* at a given surface position */

void pm_irradiance_estimate(
  float irrad[3],                // returned irradiance
  float pos[3],                  // surface position
  float normal[3],               // surface normal at pos
  float max_dist,                // max distance to look for photons
  int nphotons )                 // number of photons to use

{
  int i;
  NearestPhotons np;
  float pdir[3];
  float tmp;

  irrad[0] = irrad[1] = irrad[2] = 0.0;

  np.dist2 = (float*)malloc( sizeof(float)*(nphotons+1) );
  np.index = (Photon**)malloc( sizeof(Photon*)*(nphotons+1) );

  np.pos[0] = pos[0]; np.pos[1] = pos[1]; np.pos[2] = pos[2];
  np.max = nphotons;
  np.found = 0;
  np.got_heap = 0;
  np.dist2[0] = max_dist*max_dist;

  // locate the nearest photons
  pm_locate_photons( &np, 1 );

  // if less than 8 photons return
  if (np.found<8)
    return;


  // sum irradiance from all photons
  for (i=1; i<=np.found; i++) {
    Photon *p = np.index[i];
    // the photon_dir call and following if can be omitted (for speed)
    // if the scene does not have any thin surfaces
    pm_photon_dir( pdir, p );
    if ( (pdir[0]*normal[0]+pdir[1]*normal[1]+pdir[2]*normal[2]) < 0.0f ) {
      irrad[0] += p->power[0];
      irrad[1] += p->power[1];
      irrad[2] += p->power[2];
    }
  }

  tmp=(1.0f/M_PI)/(np.dist2[0]);    // estimate of density

  irrad[0] *= tmp;
  irrad[1] *= tmp;
  irrad[2] *= tmp;
}


/* locate_photons finds the nearest photons in the
 * photon map given the parameters in np
*/
//******************************************
void pm_locate_photons(
  NearestPhotons *np,
  long index ) 
//******************************************
{
  Photon *p = &(pm -> photons[index]);
  float dist1, dist2;

  if (index<pm -> half_stored_photons) {
    dist1 = np->pos[ p->plane ] - p->pos[ p->plane ];

    if (dist1>0.0) { // if dist1 is positive search right plane
      pm_locate_photons( np, 2*index+1 );
      if ( dist1*dist1 < np->dist2[0] )
        pm_locate_photons( np, 2*index );
    } else {         // dist1 is negative search left first
      pm_locate_photons( np, 2*index );
      if ( dist1*dist1 < np->dist2[0] )
        pm_locate_photons( np, 2*index+1 );
    }
  }

  // compute squared distance between current photon and np->pos

  dist1 = p->pos[0] - np->pos[0];
  dist2 = dist1*dist1;
  dist1 = p->pos[1] - np->pos[1];
  dist2 += dist1*dist1;
  dist1 = p->pos[2] - np->pos[2];
  dist2 += dist1*dist1;
  
  if ( dist2 < np->dist2[0] ) {
    // we found a photon  [:)] Insert it in the candidate list

    if ( np->found < np->max ) {
      // heap is not full; use array
      np->found++;
      np->dist2[np->found] = dist2;
      np->index[np->found] = p;
    } else {
      int j,parent, half_found, k;

      if (np->got_heap==0) { // Do we need to build the heap?
        // Build heap
        float dst2;
        Photon *phot;
        half_found = np->found>>1;
        for ( k=half_found; k>=1; k--) {
          parent=k;
          phot = np->index[k];
          dst2 = np->dist2[k];
          while ( parent <= half_found ) {
            j = parent+parent;
            if (j<np->found && np->dist2[j]<np->dist2[j+1])
              j++;
            if (dst2>=np->dist2[j])
              break;
            np->dist2[parent] = np->dist2[j];
            np->index[parent] = np->index[j];
            parent=j;
          }
          np->dist2[parent] = dst2;
          np->index[parent] = phot;
        }
        np->got_heap = 1;
      }

      // insert new photon into max heap
      // delete largest element, insert new and reorder the heap

      parent=1;
      j = 2;
      while ( j <= np->found ) {
        if ( j < np->found && np->dist2[j] < np->dist2[j+1] )
          j++;
        if ( dist2 > np->dist2[j] )
          break;
        np->dist2[parent] = np->dist2[j];
        np->index[parent] = np->index[j];
        parent = j;
        j += j;
      }
      np->index[parent] = p;
      np->dist2[parent] = dist2;

      np->dist2[0] = np->dist2[1];
    }
  }
}


/* store puts a photon into the flat array that will form
 * the final kd-tree.
 *
 * Call this function to store a photon.
*/

void pm_store(
   float power[3],
   float pos[3],
   float dir[3] )

{
  int i, theta, phi;
  Photon *node;

  if (pm -> stored_photons > pm -> max_photons)
    return;

  pm -> stored_photons++;
  node = &(pm -> photons[pm -> stored_photons]);

  for (i=0; i<3; i++) {
    node->pos[i] = pos[i];

    if (node->pos[i] < pm -> bbox_min[i])
      pm -> bbox_min[i] = node->pos[i];
    if (node->pos[i] > pm -> bbox_max[i])
      pm -> bbox_max[i] = node->pos[i];

    node->power[i] = power[i];
  }

  theta = (int)( acos(dir[2])*(256.0/M_PI) );
  if (theta>255)
    node->theta = 255;
  else
   node->theta = (unsigned char)theta;

  phi = (int)( atan2(dir[1],dir[0])*(256.0/(2.0*M_PI)) );
  if (phi>255)
    node->phi = 255;
  else if (phi<0)
    node->phi = (unsigned char)(phi+256);
  else
    node->phi = (unsigned char)phi;
}


/* scale_photon_power is used to scale the power of all
 * photons once they have been emitted from the light
 * source. scale = 1/(#emitted photons).
 * Call this function after each light source is processed.
*/
//********************************************************
void pm_scale_photon_power( float scale )
//********************************************************
{
  int i;

  for (i=pm -> prev_scale; i<=pm -> stored_photons; i++) {
    pm -> photons[i].power[0] *= scale;
    pm -> photons[i].power[1] *= scale;
    pm -> photons[i].power[2] *= scale;
  }
  pm -> prev_scale = pm -> stored_photons;
}


/* balance creates a left balanced kd-tree from the flat photon array.
 * This function should be called before the photon map
 * is used for rendering.
 */

void pm_balance(void)

{ 
  int i;
    int d, j=1, foo=1;
    Photon foo_photon;

  if (pm -> stored_photons>1) {
    // allocate two temporary arrays for the balancing procedure
    Photon **pa1 = (Photon**)malloc(sizeof(Photon*)*(pm -> stored_photons+1));
    Photon **pa2 = (Photon**)malloc(sizeof(Photon*)*(pm -> stored_photons+1));

    for (i=0; i<=pm -> stored_photons; i++)
      pa2[i] = &(pm -> photons[i]);

    balance_segment( pa1, pa2, 1, 1, pm -> stored_photons );
    free(pa2);

    // reorganize balanced kd-tree (make a heap)
    foo_photon = pm -> photons[j];

    for (i=1; i<=pm -> stored_photons; i++) {
      d=pa1[j]-pm -> photons;
      pa1[j] = NULL;
      if (d != foo)
        pm -> photons[j] = pm -> photons[d];
      else {
        pm -> photons[j] = foo_photon;

        if (i<pm -> stored_photons) {
          for (;foo<=pm -> stored_photons; foo++)
            if (pa1[foo] != NULL)
              break;
          foo_photon = pm -> photons[foo];
          j = foo;
        }
        continue;
      }
      j = d;
    }
    free(pa1);
  }

  pm -> half_stored_photons = pm -> stored_photons/2-1;
}


#define swap(ph,a,b) { Photon *ph2=ph[a]; ph[a]=ph[b]; ph[b]=ph2; }

// median_split splits the photon array into two separate
// pieces around the median with all photons below the
// the median in the lower half and all photons above
// than the median in the upper half. The comparison
// criteria is the axis (indicated by the axis parameter)
// (inspired by routine in "Algorithms in C++" by Sedgewick)

void median_split(
  Photon **p,
  int start,               // start of photon block in array
  int end,                 // end of photon block in array
  int median,              // desired median number
  int axis )               // axis to split along

{
  int left = start;
  int right = end;
  int i, j;

  while ( right > left ) {
    float v = p[right]->pos[axis];
    i=left-1;
    j=right;
    for (;;) {
      while ( p[++i]->pos[axis] < v )
        ;
      while ( p[--j]->pos[axis] > v && j>left )
        ;
      if ( i >= j )
        break;
      swap(p,i,j);
    }

    swap(p,i,right);
    if ( i >= median )
      right=i-1;
    if ( i <= median )
      left=i+1;
  }
}

  
// See "Realistic image synthesis using Photon Mapping" chapter 6
// for an explanation of this function

void balance_segment(
  Photon **pbal,
  Photon **porg,
  int index,
  int start,
  int end )

{
  int axis=2;

  //--------------------
  // compute new median
  //--------------------

  int median=1;
  while ((4*median) <= (end-start+1))
    median += median;

  if ((3*median) <= (end-start+1)) {
    median += median;
    median += start-1;
  } else        
    median = end-median+1;

  //--------------------------
  // find axis to split along
  //--------------------------

  if ((pm -> bbox_max[0]-pm -> bbox_min[0])>(pm -> bbox_max[1]-pm -> bbox_min[1]) &&
(pm -> bbox_max[0]-pm -> bbox_min[0])>(pm -> bbox_max[2]-pm -> bbox_min[2]))
    axis=0;
  else if ((pm -> bbox_max[1]-pm -> bbox_min[1])>(pm -> bbox_max[2]-pm -> bbox_min[2]))
    axis=1;

  //------------------------------------------
  // partition photon block around the median
  //------------------------------------------

  median_split( porg, start, end, median, axis );

  pbal[ index ] = porg[ median ];
  pbal[ index ]->plane = axis;

  //----------------------------------------------
  // recursively balance the left and right block
  //----------------------------------------------

  if ( median > start ) {
    // balance left segment
    if ( start < median-1 ) {
      float tmp=pm -> bbox_max[axis];
      pm -> bbox_max[axis] = pbal[index]->pos[axis];
      balance_segment( pbal, porg, 2*index, start, median-1 );
      pm -> bbox_max[axis] = tmp;
    } else {
      pbal[ 2*index ] = porg[start];
    }
  }

  if ( median < end ) {
    // balance right segment
    if ( median+1 < end ) {
      float tmp = pm -> bbox_min[axis];         
      pm -> bbox_min[axis] = pbal[index]->pos[axis];
      balance_segment( pbal, porg, 2*index+1, median+1, end );
      pm -> bbox_min[axis] = tmp;
    } else {
      pbal[ 2*index+1 ] = porg[end];
    }
  }     
}

