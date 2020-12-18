
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


/*********************************************************************************************

    This program is public domain software that was developed by 
    the U.S. Naval Oceanographic Office.

    This is a work of the US Government. In accordance with 17 USC 105,
    copyright protection is not available for any work of the US Government.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

*********************************************************************************************/

#ifndef VERSION

#define     VERSION     "PFM Software - chrtr2 V2.08 - 07/29/14"

#endif

/*

    Version 1.00
    Jan C. Depner
    07/20/10

    First version.


    Version 1.01
    Jan C. Depner
    03/18/11

    If no input files are specified this will create an empty CHRTR2 file for use in chrtr2_merge.


    Version 2.0
    Stacy Johnson
    08/20/2012

    Removed half node shift since chrtr2 was moved to grid registration


    Version 2.01
    Stacy Johnson
    10/01/2012

    added ceil call when calculating row/cols due to an issue going from float to int.


    Version 2.02
    Stacy Johnson
    11/14/2012

    Added subtraction on .5 to row/column calculation so that extra row/col is not added some cases


    Version 2.03
    Jan C. Depner (PFM Software)
    03/17/14

    Removed WLF support.  Top o' the mornin' to ye!


    Version 2.04
    Jan C. Depner (PFM Software)
    05/07/14

    Fixed freads and fgets without return check.


    Version 2.05
    Jan C. Depner (PFM Software)
    05/27/14

    Removed UNISIPS support.


    Version 2.06
    Jan C. Depner (PFM Software)
    07/21/14

    - Removed support for old SHOALS files.  It's been gone from pfmLoad for years
      so I don't know why it was left in here.


    Version 2.07
    Jan C. Depner (PFM Software)
    07/23/14

    - Switched from using the old NV_INT64 and NV_U_INT32 type definitions to the C99 standard stdint.h and
      inttypes.h sized data types (e.g. int64_t and uint32_t).


    Version 2.08
    Jan C. Depner (PFM Software)
    07/29/14

    - Fixed errors discovered by cppcheck.

*/
