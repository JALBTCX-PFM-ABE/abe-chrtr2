
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
*   Programmer(s):      Stephanie Lee                                       *
*                                                                           *
*   Date Written:       March 21, 2000                                      *
*                                                                           *
*   Module Name:        checkinput                                          *
*                                                                           *
*   Module Security                                                         *
*   Classification:     Unclassified                                        *
*                                                                           *
*   Data Security                                                           *
*   Classification:     Unknown                                             *
*                                                                           *
*   Purpose:            Check the input file to determine if swapping       *
*                       bytes is necessary.                                 *
*                                                                           *
*   Inputs:             Pointer to file                                     *
*                                                                           *
*   Outputs:            None                                                *
*                                                                           *
*   Files Included:     dpgextract.h                                        *
*                       stdio.h                                             *
*                       stdlib.h                                            *
*                                                                           *
*   Calling Routines:   None                                                *
*                                                                           *
*   Routines Called:    None                                                *
*                                                                           *
*   Glossary:           depth       - Depth value measured in fathoms.      *
*                       dpgptr      - Pointer to the DPG file.              *
*                       latitude    - Latitude position measured in         *
*                                     decimal minutes.                      *
*                       longitude   - Longitude position measured in        *
*                                     decimal minutes.                      *
*                       status      - keeps track of end of file.           *
*                       swap        - keeps track of whether swap bytes     *
*                                     needs to be called and returns        *
*                                     0 (don't swap) or 1 (swap)            *
*                                                                           *
*   Method:             Read each individual value from a file.             *
*                       This method was necessary because dpg files         *
*                       have no headers and the first record in the file    *
*                       will not necessarily determine if the bytes need    *
*                       to be swapped.                                      *
*                                                                           *                           
*   Modifications:      Stephanie Lee, PSI 03/04/02                         *
*			Modified to cross the equator or international      *
*                       dateline, especially redpg output and rewind file.  *
*                                                                           *
*                       Jan Depner, 02/01/07                                *
*                       Using NAVO standard data types.                     *
*                                                                           *                           
\***************************************************************************/

#include "nvutility.h"

uint8_t checkinput (FILE *dpgptr)
{
  /* Variable declaration.   */

  int32_t status = 0, n = 0;
  uint8_t  swap = NVFalse;
    
  float latitude, longitude, depth;


  /* To save processing time, only go through the first ten enteries */

  while ((status != EOF) && (!swap) && (n < 10))
    {
      /* Input the latitude, longitude, and depth of the CHRTR file.   */

      status = fread (&latitude, sizeof (float), 1, dpgptr);
      status = fread (&longitude, sizeof (float), 1, dpgptr);
      status = fread (&depth, sizeof (float), 1, dpgptr);
      latitude /= 60.0;
      longitude /= 60.0;


      /*  if latitude, longitude, or depth is too small or too large 
          or all of lat, lon and depth equal zero.  Taking crossing
          the equator or international dateline into consideration */

      if ((latitude > 180.0) || (latitude < -180.0) || (longitude > 360.0) || (longitude < -360.0) ||  
          (depth > 12000.0)) swap = NVTrue;


      /* This statement was added to deal with output from redpg when
         crossing the international date line or equator */

      if ((latitude >= 0.0) && (latitude <= 0.0000001) && (longitude >= 0.0) && (longitude <= 0.0000001) && 
          (depth >= 0.0) && (depth <= 0.0000001) && status && !swap) swap = NVTrue;


      /* Check for zero in the latitude, longitude, and depth which indicates
         an end of file.                                                      */
    
      if (((latitude == 0.0) && (longitude == 0.0) && (depth == 0.0)) || (status <= 0)) status = EOF;

      n++;
    }

  rewind (dpgptr);
  return (swap);
}
