
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
*   Programmer(s):      Jan C. Depner, Jim Braud, Arnold F. Steed,          *
*                       & Gail Smith                                        *
*                                                                           *
*   Date Written:       July 1992                                           *
*                                                                           *
*   Last Modification:  April 2000                                          *
*                                                                           *
*   Module Name:        reader                                              *
*                                                                           *
*   Module Security                                                         *
*   Classification:     Unclassified                                        *
*                                                                           *
*   Data Security                                                           *
*   Classification:     Unknown                                             *
*                                                                           *
*   Purpose:            Input all the data from the given files in the      *
*                       input parameter list.                               *
*                                                                           *
*   Inputs:             lat             -   latitude value to pass to       *
*                                           caller                          *
*                       lon             -   longitude value to pass to      *
*                                           caller                          *
*                       zvalue          -   z value to pass to caller       *
*                       date_line       -   1 if area crosses date          *
*                                           line                            *
*                                                                           *
*   Outputs:            int             -   1 on end of last data file      *
*                                                                           *
*   Calling Routines:   main                                                *
*                                                                           *
*   Glossary:           fileptr         -   Pointer for the current file    *
*                                           being processed.                *
*                       firstfile       -   Indicate that the first file is *
*                                           being processed.                *
*                       filetype    -   Indicates the type of file to       *
*                                       be read.                            */



#define         LLZ_FILE        0
#define         DPG_FILE        1
#define         RDP_FILE        2
#define         HOF_FILE        3
#define         TOF_FILE        4
#define         GSF_FILE        5
#define         YXZ_FILE        6
#define         XYZ_FILE        7
#define         PFM_FILE        8


/*                      count           -   Index counter for the record    *
*                                           array.                          *
*                       recnum          -   Current record number.          *
*                       record          -   integer input buffer.           *
*                       byte_position   -   Current byte position of file.  *
*                       byte_order      -   4 byte array indicating int     *
*                                           byte order of rdp file.         *
*                       zvalue          -   Position in z-direction of an   *
*                                           xyz-coordinate system.          *
*                                                                           *
*									    *
*   Modifications:	Sam Mangin, PSI 6/01 - Added a swap byte check      *
*                       for dpg files.                                      *
*									    *
\***************************************************************************/

#include "FileHydroOutput.h"
#include "FileTopoOutput.h"

#include "nvutility.h"

#include "gsf.h"

#include "pfm.h"

#include "llz.h"

    

static PFM_OPEN_ARGS        open_args;
static LLZ_HEADER           llz_header;


int32_t big_endian ();


static void this_swap_int (int32_t *word)
{
  uint32_t      temp[4];

  temp[0] = *word & 0x000000ff;

  temp[1] = (*word & 0x0000ff00) >> 8;

  temp[2] = (*word & 0x00ff0000) >> 16;

  temp[3] = (*word & 0xff000000) >> 24;

  *word = (temp[0] << 24) + (temp[1] << 16) + (temp[2] << 8) + temp[3];
}


static void swap_rdp (int32_t rdp_record[])
{
  this_swap_int (&rdp_record[0]);
  this_swap_int (&rdp_record[1]);
  this_swap_int (&rdp_record[2]);
}


/***************************************************************************\
*                                                                           *
*   Programmer(s):      Jan C. Depner                                       *
*                                                                           *
*   Date Written:       July 1992                                           *
*                                                                           *
*   Module Name:        openfile                                            *
*                                                                           *
*   Module Security                                                         *
*   Classification:     Unclassified                                        *
*                                                                           *
*   Data Security                                                           *
*   Classification:     Unknown                                             *
*                                                                           *
*   Purpose:            Input all the files from the parameter list.        *
*                                                                           *
*   Inputs:             file        -   array of file names                 *
*                       numfiles    -   number of input files               *
*                       fileptr     -   input file pointers                 *
*                       filetype    -   Indicates the type of file          *
*                                       see #define's above                 *
*                       handle      -   handle, if this is a GSF, PFM, or   *
*                                       LLZ file                            *
*                       open_args   -   PFM_OPEN_ARGS                       *
*                       llz_header  -   LLZ header                          *
*                                                                           *
*   Outputs:            int32_t     -   1 on end of last file               *
*                                                                           *
*   Calling Routines:   read_files                                          *
*                                                                           *
\***************************************************************************/

static int32_t openfile (char *file[], int32_t numfiles, FILE **fileptr, int32_t *filetype, int32_t *handle)
{
  static int32_t       filecount = 0, prev_filetype = -1;


  int32_t open_out_file (char *, FILE **);



  /* If it is not the first input file close it                       */

  if (prev_filetype != -1)
    {
      if (prev_filetype == GSF_FILE)
        {
          gsfClose (*handle);
        }
      else if (prev_filetype == PFM_FILE)
        {
          close_pfm_file (*handle);
        }
      else if (prev_filetype == LLZ_FILE)
        {
          close_llz (*handle);
        }
      else
        {
          fclose (*fileptr);
        }
    }

  if (filecount == numfiles)
    {
      return (1);
    }


  fprintf (stderr, "\n\nData file %03d of %03d: %s\n\n", filecount + 1, numfiles, file[filecount]);
  fflush (stderr);

  if (strstr (file[filecount], ".llz") != NULL)
    {
      *handle = open_llz (file[filecount], &llz_header);

      if (*handle < 0)
        {
          perror (file[filecount]);
          exit (-1);
        }

      *filetype = LLZ_FILE;
    }
  else if (strstr (file[filecount], ".dpg") != NULL)
    {
      *filetype = DPG_FILE;
    }
  else if (strstr (file[filecount], ".rdp") != NULL)
    {
      *filetype = RDP_FILE;
    }
  else if (strstr (file[filecount], ".hof") != NULL)
    {
      *filetype = HOF_FILE;
    }
  else if (strstr (file[filecount], ".tof") != NULL)
    {
      *filetype = TOF_FILE;
    }
  else if (strstr (file[filecount], ".txt") != NULL || strstr (file[filecount], ".yxz") != NULL || 
           strstr (file[filecount], ".raw") != NULL)
    {
      *filetype = YXZ_FILE;
    }
  else if (strstr (file[filecount], ".xyz") != NULL)
    {
      *filetype = XYZ_FILE;
    }
  else if (strstr (file[filecount], ".pfm") != NULL)
    {
      open_args.checkpoint = 0;
      strcpy (open_args.list_path, file[filecount]);

      *handle = open_existing_pfm_file (&open_args);

      if (*handle < 0) pfm_error_exit (pfm_error);

      *filetype = PFM_FILE;
    }
  else
    {
      if (gsfOpen(file[filecount], GSF_READONLY, handle) == -1)
        {
          fprintf (stderr, "\n\nUnable to open file %s\n", file[filecount]);
          gsfPrintError (stderr);
          exit (-1);
        }
      *filetype = GSF_FILE;
    }


  if (*filetype != GSF_FILE && *filetype != PFM_FILE && *filetype != LLZ_FILE)
    {
      if (*filetype == HOF_FILE)
        {
          if ((*fileptr = open_hof_file (file[filecount])) == NULL)
            {
              perror(file[filecount]);
              exit(-1);
            }
        }
      else if (*filetype == TOF_FILE)
        {
          if ((*fileptr = open_tof_file (file[filecount])) == NULL)
            {
              perror(file[filecount]);
              exit(-1);
            }
        }
      else
        {
          /* Open the input file.  */

          if ((*fileptr = fopen (file[filecount], "r")) == NULL)
            {
              perror(file[filecount]);
              exit(-1);
            }
        }
    }


  filecount++;

  prev_filetype = *filetype;

  return (0);
}



int32_t reader (NV_F64_COORD3 *xyz, int32_t date_line, char *file[], int32_t numfiles, uint8_t nominal)
{
  static FILE          *fileptr = NULL;
  static int32_t       firstfile = 1, filetype = 0, recnum = 0, old_percent = -1, eof, 
                       handle, beam_num = -1, total_beams, row = 0, col = 0, rec = 0, numrecs = 0;
  char                 string[256], cut[50];
  int32_t              status, rdp_record[3], endian, year, day, hour, minute, percent, 
                       byte_position = 0, latdeg, londeg, latmin, lonmin;
  float                dpg_record[3], second, dep, dep2;
  double               lateral, lat1, lon1, lat2, lon2, latsec, lonsec;
  static double        nlat, nlon, ang1, ang2;
  static BIN_RECORD    bin;
  static DEPTH_RECORD  *depth_record = NULL;
  NV_I32_COORD2        coord;
  int32_t              llz_handle = 0;
  LLZ_REC              llz_rec;
  HYDRO_OUTPUT_T       hof;
  HOF_HEADER_T         hof_head;
  TOPO_OUTPUT_T        tof;
  TOF_HEADER_T         tof_head;
  static gsfDataID     gsf_data_id;
  static gsfRecords    gsf_records;
  static uint8_t       byte_swap = NVFalse, bad, just_opened = NVFalse;



  void newgp (double, double, double, double, double *, double *);
  uint8_t checkinput (FILE *dpgptr);




  bad = NVTrue;

  while (bad)
    {
      if (filetype != GSF_FILE && filetype != PFM_FILE && filetype != LLZ_FILE && fileptr != NULL) 
        byte_position = ftell (fileptr);


      /*  If byte_position is less than zero this is the first file, otherwise, check if current file is finished 
          processing and open a file.    */

      if (byte_position < 0 || byte_position >= eof || (filetype != GSF_FILE && filetype != PFM_FILE && 
                                                        filetype != LLZ_FILE && fileptr == NULL))
        {
          recnum = 0;
          firstfile = openfile (file, numfiles, &fileptr, &filetype, &handle);
          if (firstfile)
            {
              printf ("\n\n\n");
              return (1);
            }

          if (filetype != GSF_FILE && filetype != PFM_FILE && filetype != LLZ_FILE)
            {
              byte_position = ftell (fileptr);
              fseek (fileptr, 0, SEEK_END);
              eof = ftell (fileptr);
              fseek (fileptr, byte_position, 0);
            }

          switch (filetype)
            {
            case DPG_FILE:
              byte_swap = checkinput(fileptr);			/*SM-ADDED*/
              fseek(fileptr,0,SEEK_SET);
              break;
    
            case GSF_FILE:
            case PFM_FILE:
            case LLZ_FILE:
              eof = 1;
              break;
    
            case RDP_FILE:
              if (!fread (&endian, 4, 1, fileptr))
		{
		  fprintf (stderr, "Bad return in file %s, function %s at line %d.  This should never happen!", __FILE__, __FUNCTION__, __LINE__ - 2);
		  fflush (stderr);
		  exit (-1);
		}
              if (endian != 0x00010203)
                {
                  byte_swap = NVTrue;
                }
              else
                {
                  byte_swap = NVFalse;
                }
              break;
    
    
            case HOF_FILE:
              hof_read_header (fileptr, &hof_head);
              break;
                
    
            case TOF_FILE:
              tof_read_header (fileptr, &tof_head);
              break;
                

            case YXZ_FILE:
            case XYZ_FILE:
              break;
            }
          just_opened = NVTrue;
        }

      bad = NVFalse;


      /* Input a record from the current file being processed.        */

      switch (filetype)
        {
        case LLZ_FILE:

          if (read_llz (llz_handle, LLZ_NEXT_RECORD, &llz_rec))
            {
              byte_position = 0;

              xyz->y = llz_rec.xy.lat;
              xyz->x = llz_rec.xy.lon;
              xyz->z = llz_rec.depth;

              if (llz_rec.status & LLZ_INVAL) bad = NVTrue;
            }
          else
            {
              byte_position = -1;
              bad = NVTrue;
            }
          break;


        case DPG_FILE:

          status = fread (dpg_record, sizeof (dpg_record), 1, fileptr);

 
          if(byte_swap)				/*SM-ADDED*/
            {
              swap_float(&dpg_record[0]);
              swap_float(&dpg_record[1]);
              swap_float(&dpg_record[2]);
            }					/*SM-ADDED*/

          if ((dpg_record[0] == 0.0 && dpg_record[1] == 0.0 && dpg_record[2] == 0.0) || status <= 0)
            {
              bad = NVTrue;
            }
          else
            {
              xyz->y = dpg_record[0];
              xyz->x = dpg_record[1];
              xyz->z = dpg_record[2];
            }
          break;


        case RDP_FILE:

          status = fread (rdp_record, sizeof (rdp_record), 1, fileptr);

          if (byte_swap) swap_rdp (rdp_record);

          if ((rdp_record[0] == 0 && rdp_record[1] == 0 && rdp_record[2] == 0) || status <= 0)
            {
              bad = NVTrue;
            }
          else
            {
              xyz->y = rdp_record[0] / 10000000.0;
              xyz->x = rdp_record[1] / 10000000.0;
              xyz->z = rdp_record[2] / 10000.0;
            }
          break;


        case HOF_FILE:
          hof_read_record (fileptr, HOF_NEXT_RECORD, &hof);


          /*  HOF uses the lower three bits of the status field for status thusly :

          bit 0 = deleted    (1) 
          bit 1 = kept       (2) 
          bit 2 = swapped    (4)      */

          if ((hof.status & AU_STATUS_DELETED_BIT) || (hof.abdc < 70) || (hof.correct_depth == -998.0))
            {
              bad = NVTrue;
            }
          else
            {
              xyz->y = hof.latitude;
              xyz->x = hof.longitude;
              xyz->z = -hof.correct_depth;
            }
          break;


        case TOF_FILE:
          tof_read_record (fileptr, TOF_NEXT_RECORD, &tof);


          /*  TOF uses the lower two bits of the status field for status thusly :
              bit 0 = first deleted    (1) 
              bit 1 = second deleted   (2) */

          if (tof.elevation_last == -998.0 || tof.conf_last < 50)
            {
              bad = NVTrue;
            }
          else
            {
              xyz->y = tof.latitude_last;
              xyz->x = tof.longitude_last;
              xyz->z = -tof.elevation_last;
            }
          break;



        case YXZ_FILE:


          /********************************************************************
           *
           * Description : Reads a Hypack yxz file in the following format:
           *
           *                DDD-MM-SS.SSSs,DDD-MM-SS.SSSs,DD.DDD
           *
           *               Example:
           *
           *                013-27-18.771N,144-37-41.383E,53.230
           *
           *               Last field is not fixed size so we can only read these
           *               not update (like I care).
           *
           *               Or it will read files derived from Hypack RAW files that
           *               are in the following format:
           *
           * YEAR DAY HH:MM:SS.SSSS   DD.ddddddddd  DDD.ddddddddd   depth   depth   depth
           * 2001 017 03:25:10.0620   13.430910402  144.660955271    7.06    8.77    8.02
           *
           *                All fields that start with a # are comments.  Third depth is
           *                tide corrected.
           *
           *               Or it will read YXZ files in the following formats:
           *
           * [signed] lat degrees decimal  [signed] lon degrees decimal  depth
           *
           *               or
           *
           * [signed] lat degrees decimal,[signed] lon degrees decimal,depth
           *
           *               Example:
           *
           *                21.000278 -157.619722 223
           *
           *               or
           *
           *                21.000278,-157.619722,223
           *
           ********************************************************************/

          if (fgets (string, sizeof (string), fileptr) == NULL)
	    {
	      fprintf (stderr, "Bad return in file %s, function %s at line %d.  This should never happen!", __FILE__, __FUNCTION__, __LINE__ - 2);
	      fflush (stderr);
	      exit (-1);
	    }

          if (string[0] == '#')
            {
              bad = NVTrue;
            }
          else
            {
              if (strchr (string, ':'))
                {
                  sscanf (string, "%d %d %d:%d:%f %lf %lf %f %f %lf", &year, &day, &hour, &minute, &second, &xyz->y, 
                          &xyz->x, &dep, &dep2, &xyz->z);
                }
              else if (string[3] == '-' && string[7] == '-')
                {
                  strncpy (cut, &string[1], 13);
                  cut[13] = 0;
                  sscanf (cut, "%03d-%02d-%lf", &latdeg, &latmin, &latsec);

                  strncpy (cut, &string[16], 13);
                  cut[13] = 0;
                  sscanf (cut, "%03d-%02d-%lf", &londeg, &lonmin, &lonsec);

                  strcpy (cut, &string[31]);
                  sscanf (cut, "%lf", &xyz->z);


                  xyz->y = (double) latdeg + (double) latmin / 60.0 + latsec / 3600.0;
                  if (strchr (string, 'S')) xyz->y = -(xyz->y);
                  xyz->x = (double) londeg + (double) lonmin / 60.0 + lonsec / 3600.0;
                  if (strchr (string, 'W')) xyz->x = -(xyz->x);
                }
              else
                {
                  if (strchr (string, ','))
                    {
                      sscanf (string, "%lf,%lf,%lf", &xyz->y, &xyz->x, &xyz->z);
                    }
                  else
                    {
                      sscanf (string, "%lf %lf %lf", &xyz->y, &xyz->x, &xyz->z);
                    }
                }
            }
          break;


        case XYZ_FILE:


          /********************************************************************
           *
           *         All fields that start with a # are comments.
           *
           *         Will read XYZ files in the following formats:
           *
           * [signed] lon degrees decimal  [signed] lat degrees decimal  depth
           *
           *               or
           *
           * [signed] lon degrees decimal,[signed] lat degrees decimal,depth
           *
           *               Example:
           *
           *                -157.619722 21.000278 223.5
           *
           *               or
           *
           *                -157.619722,21.000278,223.5
           *
           ********************************************************************/

          if (fgets (string, sizeof (string), fileptr) == NULL)
	    {
	      fprintf (stderr, "Bad return in file %s, function %s at line %d.  This should never happen!", __FILE__, __FUNCTION__, __LINE__ - 2);
	      fflush (stderr);
	      exit (-1);
	    }

          if (string[0] == '#')
            {
              bad = NVTrue;
            }
          else
            {
              if (strchr (string, ','))
                {
                  sscanf (string, "%lf,%lf,%lf", &xyz->x, &xyz->y, &xyz->z);
                }
              else
                {
                  sscanf (string, "%lf %lf %lf", &xyz->x, &xyz->y, &xyz->z);
                }
            }
          break;


        case GSF_FILE:

          bad = NVTrue;
          byte_position = 0;

          if (beam_num == -1)
            {
              status = gsfRead (handle, GSF_RECORD_SWATH_BATHYMETRY_PING, &gsf_data_id, &gsf_records, NULL, 0);

              if (status == -1)
                {
                  byte_position = -1;
                  break;
                }


              nlat = gsf_records.mb_ping.latitude;
              nlon = gsf_records.mb_ping.longitude;
              ang1 = gsf_records.mb_ping.heading + 90.0;
              ang2 = gsf_records.mb_ping.heading;

              total_beams = gsf_records.mb_ping.number_beams;


              /*  If lat is 91 or lon is 181  */
    
              if (nlat > 90.0 || nlon > 180.0 || (gsf_records.mb_ping.ping_flags & GSF_IGNORE_PING)) break;

              beam_num = 0;
            }


          if (!(gsf_records.mb_ping.beam_flags[beam_num] & GSF_IGNORE_BEAM))
            {
              if (gsf_records.mb_ping.across_track != NULL)
                {
                  lateral = gsf_records.mb_ping.across_track[beam_num];
                  newgp (nlat, nlon, ang1, lateral, &lat1, &lon1);
     
                  xyz->y = lat1;
                  xyz->x = lon1;
                }
              else
                {
                  xyz->y = nlat;
                  xyz->x = nlon;
                }

    
              /* if the along track array is present, use it */

              if (gsf_records.mb_ping.along_track != NULL)
                {
                  lateral = gsf_records.mb_ping.along_track[beam_num];
                  newgp (lat1, lon1, ang2, lateral, &lat2, &lon2);
                  xyz->y = lat2;
                  xyz->x = lon2;
                }

              if (nominal)
                {
                  if (gsf_records.mb_ping.nominal_depth != NULL)
                    {
                      xyz->z = gsf_records.mb_ping.nominal_depth[beam_num];
                    }
                  else
                    {
                      if (just_opened)
                        {
                          fprintf (stderr, "Nominal depth requested but not available, using true depth               \n\n");
                          fflush (stderr);
                          just_opened = NVFalse;
                        }
                      xyz->z = gsf_records.mb_ping.depth[beam_num];
                    }
                }
              else
                {
                  if (gsf_records.mb_ping.depth != NULL)
                    {
                      xyz->z = gsf_records.mb_ping.depth[beam_num];
                    }
                  else
                    {
                      if (just_opened)
                        {
                          fprintf (stderr, "True depth requested but not available, using nominal depth\n\n");
                          fflush (stderr);
                          just_opened = NVFalse;
                        }
                      xyz->z = gsf_records.mb_ping.nominal_depth[beam_num];
                    }
                }

              bad = NVFalse;
            }
          beam_num++;
          if (beam_num == total_beams) beam_num = -1;
          break;


        case PFM_FILE:

          bad = NVTrue;
          byte_position = 0;

          if (rec == numrecs)
            {
              rec = 0;
              col++;

              if (col >= open_args.head.bin_width)
                {
                  col = 0;
                  row++;

                  if (row == open_args.head.bin_height)
                    {
                      byte_position = -1;
                      break;
                    }
                }

              coord.y = row;
              coord.x = col;

              read_bin_record_index (handle, coord, &bin);

              if (!bin.num_soundings)
                {
                  rec = 0;
                  numrecs = 0;
                  break;
                }

              if (depth_record) free (depth_record);
              if (read_depth_array_index (handle, coord, &depth_record, &numrecs))
                {
                  rec = 0;
                  numrecs = 0;
                  break;
                }
            }


          if (!(depth_record[rec].validity & (PFM_INVAL | PFM_DELETED | PFM_REFERENCE)))
            {
              xyz->y = depth_record[rec].xyz.y;
              xyz->x = depth_record[rec].xyz.x;
              xyz->z = depth_record[rec].xyz.z;
              bad = NVFalse;
            }

          rec++;
          break;
        }


      recnum++;


      if (filetype == GSF_FILE)
        {
          percent = gsfPercent (handle);
        }
      else if (filetype == PFM_FILE)
        {
          percent = ((float) row / (float) open_args.head.bin_height) * 100.0;
        }
      else if (filetype == LLZ_FILE)
        {
          percent = ((float) recnum / (float) llz_header.number_of_records) * 100.0;
        }
      else
        {
          byte_position = ftell (fileptr);
          percent = ((float) byte_position / (float) eof) * 100.0;
        }

      if (old_percent != percent)
        {
          fprintf (stderr, "%3d%% processed - %15f   \r", percent, xyz->z);
          fflush (stderr);
          old_percent = percent;
        }
    }


  /* Check if the chart crosses over the date line.              */

  if ((date_line) && (xyz->x < 0.0))
    {
      xyz->x += 360.0;
    }

  return (0);
}
