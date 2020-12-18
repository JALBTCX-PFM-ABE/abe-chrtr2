
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

/***************************************************************************\
*                                                                           *
*   Programmer:         Jan C. Depner                                       *
*                                                                           *
*   Date Written:       April 2006                                          *
*                                                                           *
*   Purpose:            This program produces a uniform grid at a specified *
*                       grid spacing from irregularly spaced input data.    *
*                       This program is a combination of work by a number   *
*                       of people, most recently, myself and Dominic Avery. *
*                       In addition, quite a few people have added fixes    *
*                       to the routines that this program relies upon.      *
*                       These include Stephanie Lee (PSI), Steve Nosalik,   *
*                       Bill Rankin and a bunch of others.  If you'd like   *
*                       to see a bit more about the history of this program *
*                       and associated functions take a look at the         *
*                       comments in the misp_funcs.c file in the MISP       *
*                       library.                                            *
*                                                                           *
\***************************************************************************/

#include "nvutility.h"

#include "misp.h"
#include "chrtr2.h"

#include "version.h"


int32_t main (int32_t argc, char *argv[])
{
  FILE          *chp_fp;

  int32_t       i, j, k, m, error_control, gridcols, gridrows, reg_multfact, weight_factor, dn, up, bw, fw, chrtr2_hnd, row,
                numfiles, out_of_area, num_points, nibble, percent, old_percent, tmp_i;

  double        delta, y_griddeg, x_griddeg, center_x, center_y, maxvalue, minvalue, search_radius, tmp_pos, x, y,
                in_gridmin = 0.0, in_gridmeter = 0.0;

  float         *array;

  uint8_t       **val_array = NULL, input_file_flag, force_original_value = NVFalse, nominal = NVFalse, dateline, found;

  NV_F64_XYMBR  mbr;

  NV_F64_MBR    in_mbr = {0.0, 0.0, 0.0, 0.0};

  NV_F64_COORD3 xyz;

  char          chrtr2file[512], *input_filenames[4000], chp_file[512], varin[1024], info[1024];

  CHRTR2_HEADER chrtr2_header;

  CHRTR2_RECORD *chrtr2_array, chrtr2_record;


  int32_t reader (NV_F64_COORD3 *, int32_t, char *[], int32_t, uint8_t);
  void loadfiles (char *[], int32_t *);



  printf ("\n\n %s \n\n", VERSION);


  if (argc < 2)
    {
      fprintf (stderr, "\n\nUsage: %s CHRTRGUI_PARAMETER_FILE\n", argv[0]);
      fprintf (stderr, "Where:\n");
      fprintf (stderr, "\tCHRTRGUI_PARAMETER_FILE is a parameterfile created\n");
      fprintf (stderr, "\twith the chrtrGUI program (*.chp).\n\n");
      fflush (stderr);
      exit (-1);
    }


  /*  Set default values.  */

  error_control = 20;
  reg_multfact = 4;
  weight_factor = 2;
  delta = 0.05;
  maxvalue = 999999.0;
  minvalue = -999999.0;
  search_radius = 20.0;
  out_of_area = 0;
  num_points = 0;

  xyz.x = xyz.y = xyz.z = 0.0;


  strcpy (chp_file, argv[1]);

  if ((chp_fp = fopen (chp_file, "r")) == NULL)
    {
      perror (chp_file);
      exit (-1);
    }

  input_file_flag = NVFalse;
  numfiles = 0;
  while (ngets (varin, sizeof (varin), chp_fp) != NULL)
    {
      /*  Put everything to the right of the equals sign in 'info'.   */

      if (strchr (varin, '=') != NULL) strcpy (info, (strchr (varin, '=') + 1));

      if (strstr (varin, "[gridmin]")) sscanf (info, "%lf", &in_gridmin);
      if (strstr (varin, "[gridmeter]")) sscanf (info, "%lf", &in_gridmeter);
      if (strstr (varin, "[delta]")) sscanf (info, "%lf", &delta);
      if (strstr (varin, "[reg_mutfact]")) sscanf (info, "%d", &reg_multfact);
      if (strstr (varin, "[search_radius]")) sscanf (info, "%lf", &search_radius);
      if (strstr (varin, "[error_control]")) sscanf (info, "%d", &error_control);
      if (strstr (varin, "[weight_factor]")) sscanf (info, "%d", &weight_factor);
      if (strstr (varin, "[force_original_value]"))
        {
          sscanf (info, "%d", &tmp_i);
          force_original_value = (uint8_t) tmp_i;
        }
      if (strstr (varin, "[nibble_value]")) sscanf (info, "%d", &nibble);

      if (strstr (varin, "[nominal_depth]"))
        {
          sscanf (info, "%d", &tmp_i);
          nominal = (uint8_t) tmp_i;
        }
      if (strstr (varin, "[minvalue]")) sscanf (info, "%lf", &minvalue);
      if (strstr (varin, "[maxvalue]")) sscanf (info, "%lf", &maxvalue);
      if (strstr (varin, "[lat_south]")) 
        {
          posfix (info, &tmp_pos, POS_LAT);
          in_mbr.slat = tmp_pos;
        }
      if (strstr (varin, "[lat_north]"))
        {
          posfix (info, &tmp_pos, POS_LAT);
          in_mbr.nlat = tmp_pos;
        }
      if (strstr (varin, "[lon_west]"))
        {
          posfix (info, &tmp_pos, POS_LON);
          in_mbr.wlon = tmp_pos;
        }
      if (strstr (varin, "[lon_east]"))
        {
          posfix (info, &tmp_pos, POS_LON);
          in_mbr.elon = tmp_pos;
        }
      if (strstr (varin, "[output_file]")) get_string (varin, chrtr2file);

      if (input_file_flag)
        {
          if (strstr (varin, "**  End Input Files  **")) break;

          input_filenames[numfiles] = (char *) malloc (strlen (varin) + 1);
          strcpy (input_filenames[numfiles], varin);
          numfiles++;
        }

      if (strstr (varin, "**  Input Files  **")) input_file_flag = NVTrue;
    }

  if (force_original_value) weight_factor = -weight_factor;


  /*  Adjust the boundaries if the chart crosses over the dateline.  */

  dateline = NVFalse;
  if (in_mbr.wlon >= in_mbr.elon)
    {
      in_mbr.elon += 360.0;
      dateline = NVTrue;
    }


  /*  If we have meters input instead of minutes, Convert meters to the equivalent minutes of longitude at the
      center of the area.  */

  if (in_gridmeter > 0.00001)
    {
      /*  Convert from meters.    */

      center_x = in_mbr.wlon + (in_mbr.elon - in_mbr.wlon) / 2.0;
      center_y = in_mbr.slat + (in_mbr.nlat - in_mbr.slat) / 2.0;


      newgp (center_y, center_x, 90.0, in_gridmeter, &y, &x);
      x_griddeg = x - center_x;

      newgp (center_y, center_x, 0.0, in_gridmeter, &y, &x);
      y_griddeg = y - center_y;
    }
  else
    {
      x_griddeg = in_gridmin / 60.0;
      y_griddeg = in_gridmin / 60.0;
    }


  /*  Calculate grid rows of final file.  */

  /*SJ - adjust for change to grid orientation*/

  gridrows = (int32_t) ceil(((in_mbr.nlat - in_mbr.slat) / y_griddeg)-.5) + 1;


  /*  WARNING - the following code is non-functional (you can't get lats larger than 90.0).  We may use this in the
      future for doing proportional grids (by changing 90.0 to a reasonable value).  */

  /*****************************************************************************************************************/

  /*  Calculate proportional X grid size based on Y grid size at southern boundary.  If we're in the 
      southern hemisphere we'll use the Y grid size at the northern boundary.  In this way we'll have 
      a minimum bin size that is approximately square (spatially).  */
  /*
  if (in_mbr.slat >= 90.0 || in_mbr.nlat <= -90.0)
    {
      if (in_mbr.nlat < 0.0)
	{
	  invgp (NV_A0, NV_B0, in_mbr.nlat, in_mbr.wlon, in_mbr.nlat - y_griddeg, in_mbr.wlon, &dist, &az);
	}
      else
	{
	  invgp (NV_A0, NV_B0, in_mbr.slat, in_mbr.wlon, in_mbr.slat + y_griddeg, in_mbr.wlon, &dist, &az);
	}

      newgp (in_mbr.slat, in_mbr.wlon, 90.0, dist, &y, &x);

      x_griddeg = x - in_mbr.wlon;
    }
  */
  /*****************************************************************************************************************/


  /*SJ - adjust for change to grid orientation*/

  gridcols = (int32_t) ceil(((in_mbr.elon - in_mbr.wlon) / x_griddeg)-.5) + 1;


  /*  Add .ch2 extension to output file if it isn't already there.  */

  if (strstr (chrtr2file, ".ch2") == NULL) strcat (chrtr2file, ".ch2");


  /*  Populate the chrtr2 header prior to creating the file.  */

  memset (&chrtr2_header, 0, sizeof (CHRTR2_HEADER));

  strcpy (chrtr2_header.creation_software, VERSION);
  chrtr2_header.z_units = CHRTR2_METERS;
  chrtr2_header.mbr.wlon = in_mbr.wlon;
  chrtr2_header.mbr.slat = in_mbr.slat;
  chrtr2_header.width = gridcols;
  chrtr2_header.height = gridrows;
  chrtr2_header.lat_grid_size_degrees = y_griddeg;
  chrtr2_header.lon_grid_size_degrees = x_griddeg;
  chrtr2_header.min_observed_z = 0.0;
  chrtr2_header.max_observed_z = 0.0;
  chrtr2_header.min_z = minvalue;
  chrtr2_header.max_z = maxvalue;
  chrtr2_header.z_scale = 100.0;
  chrtr2_header.horizontal_uncertainty_scale = 0.0;
  chrtr2_header.vertical_uncertainty_scale = 0.0;
  chrtr2_header.uncertainty_scale = 0.0;
  chrtr2_header.max_number_of_points = 0;


  /*  Try to create and open the chrtr2 file.  */

  chrtr2_hnd = chrtr2_create_file (chrtr2file, &chrtr2_header);
  if (chrtr2_hnd < 0)
    {
      chrtr2_perror ();
      exit (-1);
    }


  printf ("\n\nOutput file: %s\n", chrtr2file);


  /*  If we had no input files we're just making an empty CHRTR2 file to be used with chrtr2_merge so we don't need to run
      MISP.  */

  if (numfiles)
    {
      /*  We're going to let MISP handle everything in zero based units of the bin size.  That is, we subtract off the
          west lon from longitudes then divide by the grid size in the X direction.  We do the same with the latitude using
          the south latitude.  This will give us values that range from 0.0 to gridcols in longitude and 0.0 to gridrows
          in latitude.  This allows us to use bins that are spatially nearly square but in terms of latitude and longitude
          would be stretched rectangles.  We'll use this to allow us to use different latitudinal and longitudinal bin sizes
          for areas in the extreme latitudes.  */

      mbr.min_x = 0.0;
      mbr.min_y = 0.0;
      mbr.max_x = (double) gridcols;
      mbr.max_y = (double) gridrows;


      /*  Initialize the MISP engine.  */


      if (misp_init (1.0, 1.0, (float) delta, reg_multfact, (float) search_radius, error_control, (float) maxvalue,
                     (float) minvalue, weight_factor, mbr)) return (-1);


      /*  Load all the data from the given input files.  */

      while (1)
        {
          if (reader (&xyz, dateline, input_filenames, numfiles, nominal)) break;


          /*  Move the lat and lon minutes into the grid domain.  */

          /*  IMPORTANT NOTE: Since MISP always wants to create a grid that has points at the corners of each cell and we want a grid
              with points at the center of each cell we're going to cheat a bit more here.  We have told MISP that the grid spacing
              is 1.0 in both directions (clever, no) so we are going to add .5 to the X and Y positions so that MISP will build a
              grid with points at the cell centers.  */


          /* we no longer need the half node shift since we moved chrtr2 to grid registration -SJ */

          xyz.x = (xyz.x - in_mbr.wlon) / x_griddeg;  /* + 0.5;*/
          xyz.y = (xyz.y - in_mbr.slat) / y_griddeg;  /* + 0.5;*/


          /*  Load data and check for out of area conditions.  */

          if (!misp_load (xyz))
            {
              out_of_area++;
            }
          else
            {
              num_points++;
            }
        }


      if (num_points == 0)
        {
          fprintf (stderr, "\n\nNo data points within specified bounds.\n");
          fprintf (stderr, "Check input area boundaries and/or min and max values.\n");
          fprintf (stderr, "Terminating!\n\n");
          fflush (stderr);
          exit (-1);
        }


      /*  Bump and grind!  */

      misp_proc ();


      /*  Process all the rows for the grid.  */

      array = (float *) calloc (gridcols + 1, sizeof (float));

      if (array == NULL)
        {
          perror ("Allocating retrieval array");
          exit (-1);
        }


      chrtr2_array = (CHRTR2_RECORD *) calloc (gridcols, sizeof (CHRTR2_RECORD));
      if (chrtr2_array ==NULL)
        {
          perror ("Allocating chrtr2_array");
          exit (-1);
        }


      if (nibble)
        {
          val_array = (uint8_t **) malloc (gridrows * sizeof (uint8_t *));

          if (val_array == NULL)
            {
              perror ("Allocating val_array");
              exit (-1);
            }


          for (i = 0 ; i < gridrows ; i++)
            {
              val_array[i] = (uint8_t *) malloc (gridcols * sizeof (uint8_t));

              if (val_array[i] == NULL)
                {
                  perror ("Allocating val_array[i]");
                  exit (-1);
                }
            }
        }


      chrtr2_header.min_observed_z = maxvalue + 1.0;
      chrtr2_header.max_observed_z = minvalue - 1.0;

      row = 0;
      while (misp_rtrv (array))
        {
          if (row >= gridrows) break;

          for (i = 0 ; i < gridcols ; i++) 
            {
              if (array[i] < chrtr2_header.min_observed_z) chrtr2_header.min_observed_z = array[i];
              if (array[i] > chrtr2_header.max_observed_z) chrtr2_header.max_observed_z = array[i];


              memset (&chrtr2_array[i], 0, sizeof (CHRTR2_RECORD));

              chrtr2_array[i].z = array[i];

              if (bit_test (array[i], 0))
                {
                  chrtr2_array[i].status = CHRTR2_REAL;

                  if (nibble) val_array[row][i] = NVTrue;
                }
              else
                {
                  chrtr2_array[i].status = CHRTR2_INTERPOLATED;

                  if (nibble) val_array[row][i] = NVFalse;
                }
            }


          /*  Note that we're using gridcols - 1 as the length.  This is beacuse MISP always places points on the corners of the cells
              and, thus, has one too many points for our purposes (we place points in the center of the cell, hence, why we added
              0.5 to the X and Y positions above).  */

          if (chrtr2_write_row (chrtr2_hnd, row, 0, gridcols - 1, chrtr2_array))
            {
              chrtr2_perror ();
              exit (-1);
            }

          row++;
        }



      if (!nibble) free (chrtr2_array);



      /*  Stop!  Nibble time!  */

      if (nibble)
        {
          /*  Nibble out the cells that aren't within our optional nibbling distance from a cell with valid data.  */

          fprintf (stderr, "\n\nNibbling                                                         \n\n");
          fflush (stderr);


          percent = 0;
          old_percent = -1;

          for (i = 0 ; i < gridrows ; i++)
            {
              dn = MAX (i - nibble, 0);
              up = MIN (i + nibble, gridrows - 1);

              for (j = 0 ; j < gridcols ; j++)
                {
                  if (!val_array[i][j])
                    {
                      bw = MAX (j - nibble, 0);
                      fw = MIN (j + nibble, gridcols - 1);

                      found = NVFalse;
                      for (k = dn ; k <= up ; k++)
                        {
                          for (m = bw ; m <= fw ; m++)
                            {
                              if (val_array[k][m])
                                {
                                  found = NVTrue;
                                  break;
                                }
                            }
                          if (found) break;
                        }

                      if (!found)
                        {
                          chrtr2_read_record_row_col (chrtr2_hnd, i, j, &chrtr2_record);

                          memset (&chrtr2_record, 0, sizeof (CHRTR2_RECORD));

                          chrtr2_write_record_row_col (chrtr2_hnd, i, j, chrtr2_record);
                        }
                    }
                }

              percent = ((float) (i) / (float) gridrows) * 100.0;
              if (old_percent != percent) 
                {
                  fprintf (stderr, "%03d%% nibbled             \r", percent);
                  fflush (stderr);
                  old_percent = percent;
                }
            }


          for (i = 0 ; i < gridrows ; i++) free (val_array[i]);
          free (val_array);


          percent = 0;
          old_percent = -1;
          chrtr2_header.min_observed_z = maxvalue + 1.0;
          chrtr2_header.max_observed_z = minvalue - 1.0;


          for (i = 0 ; i < gridrows ; i++)
            {
              chrtr2_read_row (chrtr2_hnd, i, 0, gridcols, chrtr2_array);


              /*  We need to recompute the mins and maxes.  */

              for (j = 0 ; j < gridcols ; j++)
                {
                  if (chrtr2_array[j].status)
                    {
                      if (chrtr2_array[j].z < chrtr2_header.min_observed_z) chrtr2_header.min_observed_z = chrtr2_array[j].z;
                      if (chrtr2_array[j].z > chrtr2_header.max_observed_z) chrtr2_header.max_observed_z = chrtr2_array[j].z;
                    }
                }


              percent = ((float) (i) / (float) gridrows) * 100.0;
              if (old_percent != percent) 
                {
                  fprintf (stderr, "%03d%% written             \r", percent);
                  fflush (stderr);
                  old_percent = percent;
                }
            }


          free (chrtr2_array);


          fprintf (stderr, "Completed                      \n");
          fflush (stderr);
        }
    }


  /*  Update the header.  */

  chrtr2_update_header (chrtr2_hnd, chrtr2_header);


  /*  Close the file.  */

  chrtr2_close_file (chrtr2_hnd);


  fprintf (stderr, "\n\nNorth latitude - %12.9f\n", chrtr2_header.mbr.nlat);
  fprintf (stderr, "South latitude - %12.9f\n", chrtr2_header.mbr.slat);
  fprintf (stderr, "West longitude - %12.9f\n", chrtr2_header.mbr.wlon);
  fprintf (stderr, "East longitude - %12.9f\n\n", chrtr2_header.mbr.elon);
  fflush (stderr);


  return (0);
}
