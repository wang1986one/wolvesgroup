#ifndef PM_H
#define PM_H

/* This is the photon.  The power is not compressed so the size is 28 bytes */
typedef struct Photon {
  float pos[3];                 // photon position
  short plane;                  // splitting plane for kd-tree
  unsigned char theta, phi;     // incoming direction
  float power[3];               // photon power (uncompressed)
} Photon;


/* This structure is used only to locate the nearest photons */
typedef struct NearestPhotons {
  long max;
  long found;
  int got_heap;
  float pos[3];
  float *dist2;
  Photon **index;
} NearestPhotons;

typedef struct {

  Photon *photons;

  long stored_photons;
  long half_stored_photons;
  long max_photons;
  long prev_scale;

  float costheta[256];
  float sintheta[256];
  float cosphi[256];
  float sinphi[256];

  float bbox_min[3];            // use bbox_min;
  float bbox_max[3];            // use bbox_max;
} PhotonMap;

void init_Photon_map( long max_phot );

void delete_Photon_map(void);

  void pm_store(
          float power[3],          // photon power
          float pos[3],            // photon position
          float dir[3] );          // photon direction

  void pm_scale_photon_power(
          float scale );           // 1/(number of emitted photons)

  void pm_balance(void);           // balance the kd-tree (before use!)

  void pm_irradiance_estimate(
    float irrad[3],                // returned irradiance
    float pos[3],                  // surface position
    float normal[3],               // surface normal at pos
    float max_dist,                // max distance to look for photons
    int nphotons );                // number of photons to use

  void pm_locate_photons(
    NearestPhotons *np,            // np is used to locate the photons
    long index );                   // call with index = 1

  void pm_photon_dir(
    float *dir,                    // direction of photon (returned)
    Photon *p );                   // the photon

#endif
