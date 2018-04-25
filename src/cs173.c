/**
 * @file cs173.c
 * @brief Main file for CS173 library.
 * @author - SCEC <>
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Delivers the prototype CS173 model which consists of ..
 *
 */

#include "limits.h"
#include "cs173.h"
#include "proj_api.h"

/**
 * Initializes the CS173 plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param dir The directory in which UCVM has been installed.
 * @param label A unique identifier for the velocity model.
 * @return Success or failure, if initialization was successful.
 */
int cs173_init(const char *dir, const char *label) {
	int tempVal = 0;
	char configbuf[512];
	double north_height_m = 0, east_width_m = 0, rotation_angle = 0;

	// Initialize variables.
	cs173_configuration = calloc(1, sizeof(cs173_configuration_t));
	cs173_velocity_model = calloc(1, sizeof(cs173_model_t));

	// Configuration file location.
	sprintf(configbuf, "%s/model/%s/data/config", dir, label);

	// Read the cs173_configuration file.
	if (read_cs173_configuration(configbuf, cs173_configuration) != SUCCESS)
		return FAIL;

	// Set up the iteration directory.
	sprintf(cs173_iteration_directory, "%s/model/%s/data/%s/", dir, label, cs173_configuration->model_dir);

	// Can we allocate the model, or parts of it, to memory. If so, we do.
	tempVal = try_reading_model(cs173_velocity_model);

	if (tempVal == SUCCESS) {
		fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
		fprintf(stderr, "hard disk may result in slow performance.\n");
	} else if (tempVal == FAIL) {
		print_error("No model file was found to read from.");
		return FAIL;
	}

	// We need to convert the point from lat, lon to UTM, let's set it up.
	if (!(cs173_latlon = pj_init_plus("+proj=latlong +datum=WGS84"))) {
		print_error("Could not set up latitude and longitude projection.");
		return FAIL;
	}
	if (!(cs173_utm = pj_init_plus("+proj=utm +zone=11 +ellps=clrk66 +datum=NAD27 +units=m +no_defs"))) {
		print_error("Could not set up UTM projection.");
		return FAIL;
	}

	// In order to simplify our calculations in the query, we want to rotate the box so that the bottom-left
	// corner is at (0m,0m). Our box's height is total_height_m and total_width_m. We then rotate the
	// point so that is is somewhere between (0,0) and (total_width_m, total_height_m). How far along
	// the X and Y axis determines which grid points we use for the interpolation routine.

	// Calculate the rotation angle of the box.
	north_height_m = cs173_configuration->top_left_corner_n - cs173_configuration->bottom_left_corner_n;
	east_width_m = cs173_configuration->top_left_corner_e - cs173_configuration->bottom_left_corner_e;

	// Rotation angle. Cos, sin, and tan are expensive computationally, so calculate once.
	rotation_angle = atan(east_width_m / north_height_m);

	cos_rotation_angle = cos(rotation_angle);
	sin_rotation_angle = sin(rotation_angle);

	total_height_m = sqrt(pow(cs173_configuration->top_left_corner_n - cs173_configuration->bottom_left_corner_n, 2.0f) +
						  pow(cs173_configuration->top_left_corner_e - cs173_configuration->bottom_left_corner_e, 2.0f));
	total_width_m  = sqrt(pow(cs173_configuration->top_right_corner_n - cs173_configuration->top_left_corner_n, 2.0f) +
						  pow(cs173_configuration->top_right_corner_e - cs173_configuration->top_left_corner_e, 2.0f));

	// Let everyone know that we are initialized and ready for business.
	cs173_is_initialized = 1;

	return SUCCESS;
}

#define PROJ_GEO_IPJ "+proj=latlong +datum=WGS84"
#define PROJ_GEO_OPJ "+proj=utm +zone=11 +ellps=WGS84"
int to_utm(double *lon, double *lat) {
	projPJ ipj= pj_init_plus(PROJ_GEO_IPJ);
	projPJ opj = pj_init_plus(PROJ_GEO_OPJ);
	int p = pj_transform(ipj, opj, 1, 1, lon, lat, NULL );
	return p;
}

/**
 * Queries CS173 at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, density, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @return SUCCESS or FAIL.
 */
int cs173_query(cs173_point_t *points, cs173_properties_t *data, int numpoints) {
	int i = 0;
	double point_utm_e = 0, point_utm_n = 0;
	double temp_e = 0, temp_n = 0; // holding either in deg or utm
	int load_x_coord = 0, load_y_coord = 0, load_z_coord = 0;
	double x_percent = 0, y_percent = 0, z_percent = 0;
	cs173_properties_t surrounding_points[8];

//MEI, original->	int zone = 10;
	int zone = cs173_configuration->utm_zone;
	int longlat2utm = 0;

	for (i = 0; i < numpoints; i++) {

		// We need to be below the surface to service this query.
		if (points[i].depth < 0) {
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;
			continue;
		}

		temp_e = points[i].longitude; // * DEG_TO_RAD;
		temp_n = points[i].latitude; // * DEG_TO_RAD;

/****
		// Still need to use utm_geo
		utm_geo_(&temp_e, &temp_n, &point_utm_e, &point_utm_n, &zone, &longlat2utm);
***/
                temp_e *= DEG_TO_RAD;
                temp_n *= DEG_TO_RAD;

                to_utm(&temp_e, &temp_n); // rewritten with utm

                point_utm_e = temp_e;
                point_utm_n = temp_n;

		// Point within rectangle.
		point_utm_n -= cs173_configuration->bottom_left_corner_n;
		point_utm_e -= cs173_configuration->bottom_left_corner_e;

		temp_e = point_utm_e;
		temp_n = point_utm_n;

		// We need to rotate that point, the number of degrees we calculated above.
		point_utm_e = cos_rotation_angle * temp_e - sin_rotation_angle * temp_n;
		point_utm_n = sin_rotation_angle * temp_e + cos_rotation_angle * temp_n;

		// Which point base point does that correspond to?
		load_x_coord = floor(point_utm_e / total_width_m * (cs173_configuration->nx - 1));
		load_y_coord = floor(point_utm_n / total_height_m * (cs173_configuration->ny - 1));

		// And on the Z-axis?
		load_z_coord = (cs173_configuration->depth / cs173_configuration->depth_interval - 1) -
					   floor(points[i].depth / cs173_configuration->depth_interval);

		// Are we outside the model's X and Y boundaries?
		if (load_x_coord > cs173_configuration->nx - 2 || load_y_coord > cs173_configuration->ny - 2 || load_x_coord < 0 || load_y_coord < 0) {
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;
			continue;
		}

		// Get the X, Y, and Z percentages for the bilinear or trilinear interpolation below.
		double x_interval=(cs173_configuration->nx > 1) ?
                     total_width_m / (cs173_configuration->nx-1):total_width_m;
                double y_interval=(cs173_configuration->ny > 1) ?
                     total_height_m / (cs173_configuration->ny-1):total_height_m;

                x_percent = fmod(point_utm_e, x_interval) / x_interval;
                y_percent = fmod(point_utm_n, y_interval) / y_interval;
                z_percent = fmod(points[i].depth, cs173_configuration->depth_interval) / cs173_configuration->depth_interval;

		if (load_z_coord < 1) {
			// We're below the model boundaries. Bilinearly interpolate the bottom plane and use that value.
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;
			continue;
		} else {
			// Read all the surrounding point properties.
			read_properties(load_x_coord,     load_y_coord,     load_z_coord,     &(surrounding_points[0]));	// Orgin.
			read_properties(load_x_coord + 1, load_y_coord,     load_z_coord,     &(surrounding_points[1]));	// Orgin + 1x
			read_properties(load_x_coord,     load_y_coord + 1, load_z_coord,     &(surrounding_points[2]));	// Orgin + 1y
			read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord,     &(surrounding_points[3]));	// Orgin + x + y, forms top plane.
			read_properties(load_x_coord,     load_y_coord,     load_z_coord - 1, &(surrounding_points[4]));	// Bottom plane origin
			read_properties(load_x_coord + 1, load_y_coord,     load_z_coord - 1, &(surrounding_points[5]));	// +1x
			read_properties(load_x_coord,     load_y_coord + 1, load_z_coord - 1, &(surrounding_points[6]));	// +1y
			read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[7]));	// +x +y, forms bottom plane.

			trilinear_interpolation(x_percent, y_percent, z_percent, surrounding_points, &(data[i]));
		}

		// Calculate Qp and Qs.
		if (data[i].vs < 1500)
			data[i].qs = data[i].vs * 0.02;
		else
			data[i].qs = data[i].vs * 0.10;

		data[i].qp = data[i].qs * 1.5;
	}

	return SUCCESS;
}

/**
 * Retrieves the material properties (whatever is available) for the given data point, expressed
 * in x, y, and z co-ordinates.
 *
 * @param x The x coordinate of the data point.
 * @param y The y coordinate of the data point.
 * @param z The z coordinate of the data point.
 * @param data The properties struct to which the material properties will be written.
 */
void read_properties(int x, int y, int z, cs173_properties_t *data) {
	// Set everything to -1 to indicate not found.
	data->vp = -1;
	data->vs = -1;
	data->rho = -1;
	data->qp = -1;
	data->qs = -1;

	float *ptr = NULL;
	FILE *fp = NULL;
        long location = 0;

        // the z is inverted at line #145
        if ( strcmp(cs173_configuration->seek_axis, "fast-y") == 0 ||
                 strcmp(cs173_configuration->seek_axis, "fast-Y") == 0 ) { // fast-y,  cs173 
            if(strcmp(cs173_configuration->seek_direction, "bottom-up") == 0) { 
	        location = ((long) z * cs173_configuration->nx * cs173_configuration->ny) + (x * cs173_configuration->ny) + y;
                } else {
                    // nz starts from 0 up to nz-1
	            location = ((long)((cs173_configuration->nz -1) - z) * cs173_configuration->nx * cs173_configuration->ny) + (x * cs173_configuration->ny) + y;
            }
        } else {  // fast-X, cca data
            if ( strcmp(cs173_configuration->seek_axis, "fast-x") == 0 ||
                     strcmp(cs173_configuration->seek_axis, "fast-X") == 0 ) { // fast-y,  cs173 
                if(strcmp(cs173_configuration->seek_direction, "bottom-up") == 0) { 
 		    location = ((long)z * cs173_configuration->nx * cs173_configuration->ny) + (y * cs173_configuration->nx) + x;
                    } else { // bottom-up
 		        location = ((long)(cs173_configuration->nz - z) * cs173_configuration->nx * cs173_configuration->ny) + (y * cs173_configuration->nx) + x;
                }
            }
        }

	// Check our loaded components of the model.
	if (cs173_velocity_model->vs_status == 2) {
		// Read from memory.
		ptr = (float *)cs173_velocity_model->vs;
		data->vs = ptr[location];
	} else if (cs173_velocity_model->vs_status == 1) {
		// Read from file.
		fp = (FILE *)cs173_velocity_model->vs;
		fseek(fp, location * sizeof(float), SEEK_SET);
                float temp;
		fread(&(temp), sizeof(float), 1, fp);
                data->vs = temp;
	}

	// Check our loaded components of the model.
	if (cs173_velocity_model->vp_status == 2) {
		// Read from memory.
		ptr = (float *)cs173_velocity_model->vp;
		data->vp = ptr[location];
	} else if (cs173_velocity_model->vp_status == 1) {
		// Read from file.
		fp = (FILE *)cs173_velocity_model->vp;
		fseek(fp, location * sizeof(float), SEEK_SET);
                float temp;
		fread(&(temp), sizeof(float), 1, fp);
                data->vp=temp;
	}

	// Check our loaded components of the model.
	if (cs173_velocity_model->rho_status == 2) {
		// Read from memory.
		ptr = (float *)cs173_velocity_model->rho;
		data->rho = ptr[location];
	} else if (cs173_velocity_model->rho_status == 1) {
		// Read from file.
		fp = (FILE *)cs173_velocity_model->rho;
		fseek(fp, location * sizeof(float), SEEK_SET);
                float temp;
		fread(&(temp), sizeof(float), 1, fp);
                data->rho=temp;
	}
}

/**
 * Trilinearly interpolates given a x percentage, y percentage, z percentage and a cube of
 * data properties in top origin format (top plane first, bottom plane second).
 *
 * @param x_percent X percentage
 * @param y_percent Y percentage
 * @param z_percent Z percentage
 * @param eight_points Eight surrounding data properties
 * @param ret_properties Returned data properties
 */
void trilinear_interpolation(double x_percent, double y_percent, double z_percent,
							 cs173_properties_t *eight_points, cs173_properties_t *ret_properties) {
	cs173_properties_t *temp_array = calloc(2, sizeof(cs173_properties_t));
	cs173_properties_t *four_points = eight_points;

	bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[0]);

	// Now advance the pointer four "cvms5_properties_t" spaces.
	four_points += 4;

	// Another interpolation.
	bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[1]);

	// Now linearly interpolate between the two.
	linear_interpolation(z_percent, &temp_array[0], &temp_array[1], ret_properties);

	free(temp_array);
}

/**
 * Bilinearly interpolates given a x percentage, y percentage, and a plane of data properties in
 * origin, bottom-right, top-left, top-right format.
 *
 * @param x_percent X percentage.
 * @param y_percent Y percentage.
 * @param four_points Data property plane.
 * @param ret_properties Returned data properties.
 */
void bilinear_interpolation(double x_percent, double y_percent, cs173_properties_t *four_points, cs173_properties_t *ret_properties) {
	cs173_properties_t *temp_array = calloc(2, sizeof(cs173_properties_t));
	linear_interpolation(x_percent, &four_points[0], &four_points[1], &temp_array[0]);
	linear_interpolation(x_percent, &four_points[2], &four_points[3], &temp_array[1]);
	linear_interpolation(y_percent, &temp_array[0], &temp_array[1], ret_properties);
	free(temp_array);
}

/**
 * Linearly interpolates given a percentage from x0 to x1, a data point at x0, and a data point at x1.
 *
 * @param percent Percent of the way from x0 to x1 (from 0 to 1 interval).
 * @param x0 Data point at x0.
 * @param x1 Data point at x1.
 * @param ret_properties Resulting data properties.
 */
void linear_interpolation(double percent, cs173_properties_t *x0, cs173_properties_t *x1, cs173_properties_t *ret_properties) {
	ret_properties->vp  = (1 - percent) * x0->vp  + percent * x1->vp;
	ret_properties->vs  = (1 - percent) * x0->vs  + percent * x1->vs;
	ret_properties->rho = (1 - percent) * x0->rho + percent * x1->rho;
	ret_properties->qp  = (1 - percent) * x0->qp  + percent * x1->qp;
	ret_properties->qs  = (1 - percent) * x0->qs  + percent * x1->qs;
}

/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return SUCCESS
 */
int cs173_finalize() {
	pj_free(cs173_latlon);
	pj_free(cs173_utm);

	if (cs173_velocity_model) free(cs173_velocity_model);
	if (cs173_configuration) free(cs173_configuration);

	return SUCCESS;
}

/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int cs173_version(char *ver, int len)
{
  int verlen;
  verlen = strlen(version_string);
  if (verlen > len - 1) {
    verlen = len - 1;
  }
  memset(ver, 0, len);
  strncpy(ver, version_string, verlen);
  return 0;
}

/**
 * Reads the cs173_configuration file describing the various properties of CVM-S5 and populates
 * the cs173_configuration struct. This assumes cs173_configuration has been "calloc'ed" and validates
 * that each value is not zero at the end.
 *
 * @param file The cs173_configuration file location on disk to read.
 * @param config The cs173_configuration struct to which the data should be written.
 * @return Success or failure, depending on if file was read successfully.
 */
int read_cs173_configuration(char *file, cs173_configuration_t *config) {
	FILE *fp = fopen(file, "r");
	char key[40];
	char value[80];
	char line_holder[128];

	// If our file pointer is null, an error has occurred. Return fail.
	if (fp == NULL) {
		print_error("Could not open the cs173_configuration file.");
		return FAIL;
	}

	// Read the lines in the cs173_configuration file.
	while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
		if (line_holder[0] != '#' && line_holder[0] != ' ' && line_holder[0] != '\n') {
			sscanf(line_holder, "%s = %s", key, value);

			// Which variable are we editing?
			if (strcmp(key, "utm_zone") == 0) 				config->utm_zone = atoi(value);
			if (strcmp(key, "model_dir") == 0)				sprintf(config->model_dir, "%s", value);
			if (strcmp(key, "nx") == 0) 	  				config->nx = atoi(value);
			if (strcmp(key, "ny") == 0) 	  			 	config->ny = atoi(value);
			if (strcmp(key, "nz") == 0) 	  			 	config->nz = atoi(value);
			if (strcmp(key, "depth") == 0) 	  			 	config->depth = atof(value);
			if (strcmp(key, "top_left_corner_e") == 0) 		config->top_left_corner_e = atof(value);
			if (strcmp(key, "top_left_corner_n") == 0)		config->top_left_corner_n = atof(value);
			if (strcmp(key, "top_right_corner_e") == 0)		config->top_right_corner_e = atof(value);
			if (strcmp(key, "top_right_corner_n") == 0)		config->top_right_corner_n = atof(value);
			if (strcmp(key, "bottom_left_corner_e") == 0)	config->bottom_left_corner_e = atof(value);
			if (strcmp(key, "bottom_left_corner_n") == 0)	config->bottom_left_corner_n = atof(value);
			if (strcmp(key, "bottom_right_corner_e") == 0)	config->bottom_right_corner_e = atof(value);
			if (strcmp(key, "bottom_right_corner_n") == 0)	config->bottom_right_corner_n = atof(value);
			if (strcmp(key, "depth_interval") == 0)		config->depth_interval = atof(value);
			if (strcmp(key, "seek_axis") == 0)		sprintf(config->seek_axis, "%s", value);
			if (strcmp(key, "seek_direction") == 0)		sprintf(config->seek_direction, "%s", value);
		}
	}

	// Have we set up all cs173_configuration parameters?
	if (config->utm_zone == 0 || config->nx == 0 || config->ny == 0 || config->nz == 0 || config->model_dir[0] == '\0' || 
		config->seek_direction[0] == '\0' || config->seek_axis[0] == '\0' ||
		config->top_left_corner_e == 0 || config->top_left_corner_n == 0 || config->top_right_corner_e == 0 ||
		config->top_right_corner_n == 0 || config->bottom_left_corner_e == 0 || config->bottom_left_corner_n == 0 ||
		config->bottom_right_corner_e == 0 || config->bottom_right_corner_n == 0 || config->depth == 0 ||
		config->depth_interval == 0) {
		print_error("One cs173_configuration parameter not specified. Please check your cs173_configuration file.");
		return FAIL;
	}

	fclose(fp);

	return SUCCESS;
}

/**
 * Prints the error string provided.
 *
 * @param err The error string to print out to stderr.
 */
void print_error(char *err) {
	fprintf(stderr, "An error has occurred while executing CS173. The error was:\n\n");
	fprintf(stderr, "%s", err);
	fprintf(stderr, "\n\nPlease contact software@scec.org and describe both the error and a bit\n");
	fprintf(stderr, "about the computer you are running CS173 on (Linux, Mac, etc.).\n");
}

/**
 * Check if the data is too big to be loaded internally (exceed maximum
 * allowable by a INT variable)
 *
 */
int too_big() {
        long max_size= (long) (cs173_configuration->nx) * cs173_configuration->ny * cs173_configuration->nz;
        long delta= max_size - INT_MAX;

	if( delta > 0) {
	    return 1;
        } else {
	    return 0;
        }
}

/**
 * Tries to read the model into memory.
 *
 * @param model The model parameter struct which will hold the pointers to the data either on disk or in memory.
 * @return 2 if all files are read to memory, SUCCESS if file is found but at least 1
 * is not in memory, FAIL if no file found.
 */
int try_reading_model(cs173_model_t *model) {
	double base_malloc = cs173_configuration->nx * cs173_configuration->ny * cs173_configuration->nz * sizeof(float);
	int file_count = 0;
	int all_read_to_memory =0;
	char current_file[128];
	FILE *fp;

	// Let's see what data we actually have.
	sprintf(current_file, "%s/vp.dat", cs173_iteration_directory);
	if (access(current_file, R_OK) == 0) {
                if( !too_big() ) { // only if fit
		    model->vp = malloc(base_malloc);
		    if (model->vp != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->vp, 1, base_malloc, fp);
                        all_read_to_memory++;
			fclose(fp);
			model->vp_status = 2;
		    } else {
			model->vp = fopen(current_file, "rb");
			model->vp_status = 1;
                    }
		} else {
		    model->vp = fopen(current_file, "rb");
		    model->vp_status = 1;
                }
		file_count++;
	}

	sprintf(current_file, "%s/vs.dat", cs173_iteration_directory);
	if (access(current_file, R_OK) == 0) {
                if( !too_big( )) { // only if fit
		    model->vs = malloc(base_malloc);
		    if (model->vs != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->vs, 1, base_malloc, fp);
                        all_read_to_memory++;
			fclose(fp);
			model->vs_status = 2;
		    } else {
			model->vs = fopen(current_file, "rb");
			model->vs_status = 1;
		    }
		} else {
		    model->vs = fopen(current_file, "rb");
		    model->vs_status = 1;
                }
		file_count++;
	}

	sprintf(current_file, "%s/density.dat", cs173_iteration_directory);
	if (access(current_file, R_OK) == 0) {
                if(!too_big() ) { // only if fit
		    model->rho = malloc(base_malloc);
		    if (model->rho != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->rho, 1, base_malloc, fp);
                        all_read_to_memory++;
			fclose(fp);
			model->rho_status = 2;
		    } else {
			model->rho = fopen(current_file, "rb");
			model->rho_status = 1;
		    }
		} else {
		    model->rho = fopen(current_file, "rb");
		    model->rho_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/qp.dat", cs173_iteration_directory);
	if (access(current_file, R_OK) == 0) {
                if( !too_big() ) { // only if fit
		    model->qp = malloc(base_malloc);
		    if (model->qp != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->qp, 1, base_malloc, fp);
                        all_read_to_memory++;
			fclose(fp);
			model->qp_status = 2;
		    } else {
			model->qp = fopen(current_file, "rb");
			model->qp_status = 1;
		    }
		} else {
		    model->qp = fopen(current_file, "rb");
		    model->qp_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/qs.dat", cs173_iteration_directory);
	if (access(current_file, R_OK) == 0) {
                if( !too_big() ) { // only if fit
		    model->qs = malloc(base_malloc);
		    if (model->qs != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->qs, 1, base_malloc, fp);
                        all_read_to_memory++;
			model->qs_status = 2;
		    } else {
			model->qs = fopen(current_file, "rb");
			model->qs_status = 1;
		    }
		} else {
		    model->qs = fopen(current_file, "rb");
		    model->qs_status = 1;
		}
		file_count++;
	}

	if (file_count == 0)
		return FAIL;
	else if (file_count > 0 && all_read_to_memory != file_count)
		return SUCCESS;
	else
		return 2;
}

// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#ifdef DYNAMIC_LIBRARY

/**
 * Init function loaded and called by the UCVM library. Calls cs173_init.
 *
 * @param dir The directory in which UCVM is installed.
 * @return Success or failure.
 */
int model_init(const char *dir, const char *label) {
	return cs173_init(dir, label);
}

/**
 * Query function loaded and called by the UCVM library. Calls cs173_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return Success or fail.
 */
int model_query(cs173_point_t *points, cs173_properties_t *data, int numpoints) {
	return cs173_query(points, data, numpoints);
}

/**
 * Finalize function loaded and called by the UCVM library. Calls cs173_finalize.
 *
 * @return Success
 */
int model_finalize() {
	return cs173_finalize();
}

/**
 * Version function loaded and called by the UCVM library. Calls cs173_version.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_version(char *ver, int len) {
	return cs173_version(ver, len);
}


int (*get_model_init())(const char *, const char *) {
        return &cs173_init;
}
int (*get_model_query())(cs173_point_t *, cs173_properties_t *, int) {
         return &cs173_query;
}
int (*get_model_finalize())() {
         return &cs173_finalize;
}
int (*get_model_version())(char *, int) {
         return &cs173_version;
}



#endif
