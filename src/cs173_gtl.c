/**
 * @file cs173_gtl.c
 *
 * @section DESCRIPTION
 *
 */

#include "cs173.h"
#include "cs173_gtl.h"

/** Location of the ucvm.e e-tree file. */
char cs173_vs30_etree_file[128];

/** Holds the configuration parameters for the Vs30 map. */
cs173_vs30_map_config_t *cs173_vs30_map;

/** Proj.4 Vs30 map projection holder. */
projPJ cs173_aeqd;

/** The cosine of the Vs30 map's rotation. */
double cs173_cos_vs30_rotation_angle = 0;
/** The sine of the Vs30 map's rotation. */
double cs173_sin_vs30_rotation_angle = 0;


/**
 * Reads the format of the Vs30 data e-tree. This file location is typically specified
 * in the cs173_configuration file of the model.
 *
 * @param filename The e-tree's file location from which to read.
 * @param map The outputted map cs173_configuration structure.
 */
int cs173_read_vs30_map(char *filename, cs173_vs30_map_config_t *map) {
	char appmeta[512];
	char *token;
	int index = 0, retVal = 0;
	map->vs30_map = etree_open(filename, O_RDONLY, 64, 0, 3);
	retVal = snprintf(appmeta, sizeof(appmeta), "%s", etree_getappmeta(map->vs30_map));

	if (retVal >= 0 && retVal < 128) {
		return FAIL;
	}

	// Now we need to parse the map cs173_configuration.
	index = 0;
	token = strtok(appmeta, "|");

	while (token != NULL) {
		switch (index) {
	    case 0:
	    	snprintf(map->type, sizeof(map->type), "%s", token);
	    	break;
	    case 1:
	    	snprintf(map->description, sizeof(map->description), "%s", token);
	    	break;
	    case 2:
	    	snprintf(map->author, sizeof(map->author), "%s", token);
	    	break;
	    case 3:
	    	snprintf(map->date, sizeof(map->date), "%s", token);
	    	break;
	    case 4:
	    	sscanf(token, "%lf", &(map->spacing));
	    	break;
	    case 5:
	    	snprintf(map->schema, sizeof(map->schema), "%s", token);
	    	break;
	    case 6:
	    	snprintf(map->projection, sizeof(map->projection), "%s", token);
	    	break;
	    case 7:
	    	sscanf(token, "%lf,%lf,%lf", &(map->origin_point.longitude), &(map->origin_point.latitude),
	    			&(map->origin_point.depth));
	    	break;
	    case 8:
	    	sscanf(token, "%lf", &(map->rotation));
	    	break;
	    case 9:
	    	sscanf(token, "%lf,%lf,%lf", &(map->x_dimension), &(map->y_dimension), &(map->z_dimension));
	    	break;
	    case 10:
	    	sscanf(token, "%u,%u,%u", &(map->x_ticks), &(map->y_ticks), &(map->z_ticks));
	    	break;
	    default:
	    	fprintf(stderr, "Unexpected metadata. Please check your Vs30 e-tree within UCVM.\n");
	    	return FAIL;
	    	break;
		}
	    index++;
	    token = strtok(NULL, "|");
	}

	return SUCCESS;

}

/**
 * Given a latitude and longitude in WGS84 co-ordinates, we find the corresponding e-tree octant
 * in the Vs30 map e-tree and read the value as well as interpolate bilinearly.
 *
 * @param longitude The longitude in WGS84 format.
 * @param latitude The latitude in WGS84 format.
 * @param map The Vs30 map structure as defined during the initialization procedure.
 * @return The Vs30 value at that point, or -1 if outside the boundaries.
 */
double cs173_get_vs30_value(double longitude, double latitude, cs173_vs30_map_config_t *map) {
	// Convert both points to UTM.
	double longitude_utm_e = longitude * DEG_TO_RAD;
	double latitude_utm_n = latitude * DEG_TO_RAD;
	double vs30_long_utm_e = map->origin_point.longitude * DEG_TO_RAD;
	double vs30_lat_utm_n = map->origin_point.latitude * DEG_TO_RAD;
	double temp_rotated_point_n = 0.0, temp_rotated_point_e = 0.0;
	double rotated_point_n = 0.0, rotated_point_e = 0.0;
	double percent = 0.0;
	int loc_x = 0, loc_y = 0;
	etree_addr_t addr;
	cs173_vs30_mpayload_t vs30_payload[4];

	int max_level = ceil(log(map->x_dimension / map->spacing) / log(2.0));

	etree_tick_t edgetics = (etree_tick_t)1 << (ETREE_MAXLEVEL - max_level);
	double map_edgesize = map->x_dimension / (double)((etree_tick_t)1<<max_level);

	pj_transform(cs173_latlon, cs173_aeqd, 1, 1, &longitude_utm_e, &latitude_utm_n, NULL);
	pj_transform(cs173_latlon, cs173_aeqd, 1, 1, &vs30_long_utm_e, &vs30_lat_utm_n, NULL);

	// Now that both are in UTM, we can subtract and rotate.
	temp_rotated_point_e = longitude_utm_e - vs30_long_utm_e;
	temp_rotated_point_n = latitude_utm_n - vs30_lat_utm_n;

	rotated_point_e = cs173_cos_vs30_rotation_angle * temp_rotated_point_e - cs173_sin_vs30_rotation_angle * temp_rotated_point_n;
	rotated_point_n = cs173_sin_vs30_rotation_angle * temp_rotated_point_e + cs173_cos_vs30_rotation_angle * temp_rotated_point_n;

	// Are we within the box?
	if (rotated_point_e < 0 || rotated_point_n < 0 || rotated_point_e > map->x_dimension ||
		rotated_point_n > map->y_dimension) return -1;

	// Get the integer location of the grid point within the map.
	loc_x = floor(rotated_point_e / map_edgesize);
	loc_y = floor(rotated_point_n / map_edgesize);

	// We need the four surrounding points for bilinear interpolation.
	addr.level = ETREE_MAXLEVEL;
	addr.x = loc_x * edgetics; addr.y = loc_y * edgetics; addr.z = 0;
    /* Adjust addresses for edges of grid */
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[0]));
	addr.x = (loc_x + 1) * edgetics; addr.y = loc_y * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[1]));
	addr.x = loc_x * edgetics; addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[2]));
	addr.x = (loc_x + 1) * edgetics; addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[3]));

	percent = fmod(rotated_point_e / map->spacing, map->spacing) / map->spacing;
	vs30_payload[0].vs30 = percent * vs30_payload[0].vs30 + (1 - percent) * vs30_payload[1].vs30;
	vs30_payload[1].vs30 = percent * vs30_payload[2].vs30 + (1 - percent) * vs30_payload[3].vs30;

	return vs30_payload[0].vs30;
}

/**
 * Gets the GTL value using the Wills and Wald dataset, given a latitude, longitude and depth.
 *
 * @param point The point at which to retrieve the property. Note, depth is ignored.
 * @param data The material properties at the point specified, or -1 if not found.
 * @return Success or failure.
 */
int cs173_get_vs30_based_gtl(cs173_point_t *point, cs173_properties_t *data) {
        double a = 0.5, b = 0.6, c = 0.5;
	double percent_z = point->depth / cs173_configuration->depth_interval;
	double f = 0.0, g = 0.0;
	double vs30 = 0.0, vp30 = 0.0;

	// Double check that we're above the first layer.
	if (percent_z > 1) return FAIL;

	// Query for the point at depth_interval.
	cs173_point_t *pt = calloc(1, sizeof(cs173_point_t));
	cs173_properties_t *dt = calloc(1, sizeof(cs173_properties_t));

	pt->latitude = point->latitude;
	pt->longitude = point->longitude;
	pt->depth = cs173_configuration->depth_interval;

	if (cs173_query(pt, dt, 1) != SUCCESS) return FAIL;

	// Now we need the Vs30 data value.
	vs30 = cs173_get_vs30_value(point->longitude, point->latitude, cs173_vs30_map);

	if (vs30 == -1) {
		data->vp = -1;
		data->vs = -1;
	} else {
		// Get the point's material properties within the GTL.
		f = percent_z + b * (percent_z - pow(percent_z, 2.0f));
		g = a - a * percent_z + c * (pow(percent_z, 2.0f) + 2.0 * sqrt(percent_z) - 3.0 * percent_z);
		data->vs = f * dt->vs + g * vs30;
		vs30 = vs30 / 1000;
		vp30 = 0.9409 + 2.0947 * vs30 - 0.8206 * pow(vs30, 2.0f) + 0.2683 * pow(vs30, 3.0f) - 0.0251 * pow(vs30, 4.0f);
		vp30 = vp30 * 1000;
		data->vp = f * dt->vp + g * vp30;
	}

	free(pt);
	free(dt);

	return SUCCESS;
}

