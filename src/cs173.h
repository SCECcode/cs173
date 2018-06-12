/**
 * @file cs173.h
 * @brief Main header file for CS173 library.
 * @author SCEC <>
 * @version 1.0
 *
 * Delivers the CS173 model which consists ..
 *
 */

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "proj_api.h"

// Constants
#ifndef M_PI
	/** Defines pi */
	#define M_PI 3.14159265358979323846
#endif

/** Defines a return value of success */
#define SUCCESS 0
/** Defines a return value of failure */
#define FAIL 1

// Structures
/** Defines a point (latitude, longitude, and depth) in WGS84 format */
typedef struct cs173_point_t {
	/** Longitude member of the point */
	double longitude;
	/** Latitude member of the point */
	double latitude;
	/** Depth member of the point */
	double depth;
} cs173_point_t;

/** Defines the material properties this model will retrieve. */
typedef struct cs173_properties_t {
	/** P-wave velocity in meters per second */
	double vp;
	/** S-wave velocity in meters per second */
	double vs;
	/** Density in g/m^3 */
	double rho;
	/** Qp */
	double qp;
	/** Qs */
	double qs;
} cs173_properties_t;

/** The CS173 configuration structure. */
typedef struct cs173_configuration_t {
	/** The zone of UTM projection */
	int utm_zone;
	/** The model directory */
	char model_dir[128];
        /** GTL on or off (1 or 0) */
        int gtl;
	/** Number of x points */
	int nx;
	/** Number of y points */
	int ny;
	/** Number of z points */
	int nz;
	/** Depth in meters */
	double depth;
	/** Top left corner easting in UTM projection */
	double top_left_corner_e;
	/** Top left corner northing in UTM projection */
	double top_left_corner_n;
	/** Top right corner easting in UTM projection */
	double top_right_corner_e;
	/** Top right corner northing in UTM projection */
	double top_right_corner_n;
	/** Bottom left corner easting in UTM projection */
	double bottom_left_corner_e;
	/** Bottom left corner northing in UTM projection */
	double bottom_left_corner_n;
	/** Bottom right corner easting in UTM projection */
	double bottom_right_corner_e;
	/** Bottom right corner northing in UTM projection */
	double bottom_right_corner_n;
	/** Z interval for the data */
	double depth_interval;
	/** The data access seek method, fast-X, or fast-Y */
	char seek_axis[128];
	/** The data seek direction, bottom-up, or top-down */
	char seek_direction[128];
	/**  using vs or vp*/
	char density[128];
        /** Brocher 2005 scaling polynomial coefficient 10^0 */
        double p0;
        /** Brocher 2005 scaling polynomial coefficient 10^1 */
        double p1;
        /** Brocher 2005 scaling polynomial coefficient 10^2 */
        double p2;
        /** Brocher 2005 scaling polynomial coefficient 10^3 */
        double p3;
        /** Brocher 2005 scaling polynomial coefficient 10^4 */
        double p4;
        /** Brocher 2005 scaling polynomial coefficient 10^5 */
        double p5;
} cs173_configuration_t;

/** The model structure which points to available portions of the model. */
typedef struct cs173_model_t {
	/** A pointer to the Vs data either in memory or disk. Null if does not exist. */
	void *vs;
	/** Vs status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int vs_status;
	/** A pointer to the Vp data either in memory or disk. Null if does not exist. */
	void *vp;
	/** Vp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int vp_status;
	/** A pointer to the rho data either in memory or disk. Null if does not exist. */
	void *rho;
	/** Rho status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int rho_status;
	/** A pointer to the Qp data either in memory or disk. Null if does not exist. */
	void *qp;
	/** Qp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int qp_status;
	/** A pointer to the Qs data either in memory or disk. Null if does not exist. */
	void *qs;
	/** Qs status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int qs_status;
} cs173_model_t;


// UCVM API Required Functions

#ifdef DYNAMIC_LIBRARY

/** Initializes the model */
int model_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int model_finalize();
/** Returns version information */
int model_version(char *ver, int len);
/** Queries the model */
int model_query(cs173_point_t *points, cs173_properties_t *data, int numpts);

int (*get_model_init())(const char *, const char *);
int (*get_model_query())(cs173_point_t *, cs173_properties_t *, int);
int (*get_model_finalize())();
int (*get_model_version())(char *, int);


#endif

// CS173 Related Functions

/** Initializes the model */
int cs173_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int cs173_finalize();
/** Returns version information */
int cs173_version(char *ver, int len);
/** Queries the model */
int cs173_query(cs173_point_t *points, cs173_properties_t *data, int numpts);

// Non-UCVM Helper Functions
/** Reads the configuration file. */
int read_configuration(char *file, cs173_configuration_t *config);
/** Prints out the error string. */
void cs173_print_error(char *err);
/** Retrieves the value at a specified grid point in the model. */
void cs173_read_properties(int x, int y, int z, cs173_properties_t *data);
/** Attempts to malloc the model size in memory and read it in. */
int cs173_try_reading_model(cs173_model_t *model);
/** Calculates density from Vs. */
double cs173_calculate_density(double vs);
/** Calculates density from Vp. */
double cs173_nafe_drake_rho(double vp);

// Interpolation Functions
/** Linearly interpolates two cs173_properties_t structures */
void cs173_linear_interpolation(double percent, cs173_properties_t *x0, cs173_properties_t *x1, cs173_properties_t *ret_properties);
/** Bilinearly interpolates the properties. */
void cs173_bilinear_interpolation(double x_percent, double y_percent, cs173_properties_t *four_points, cs173_properties_t *ret_properties);
/** Trilinearly interpolates the properties. */
void cs173_trilinear_interpolation(double x_percent, double y_percent, double z_percent, cs173_properties_t *eight_points,
							 cs173_properties_t *ret_properties);


/** Configuration parameters. */
extern cs173_configuration_t *cs173_configuration;
/** Proj.4 latitude longitude, WGS84 projection holder. */
extern projPJ cs173_latlon;
/** Proj.4 UTM projection holder. */
extern projPJ cs173_utm;

/** The cosine of the rotation angle used to rotate the box and point around the bottom-left corner. */
extern double cs173_cos_rotation_angle;
/** The sine of the rotation angle used to rotate the box and point around the bottom-left corner. */
extern double cs173_sin_rotation_angle;

/** The height of this model's region, in meters. */
extern double cs173_total_height_m;
/** The width of this model's region, in meters. */
extern double cs173_total_width_m;

