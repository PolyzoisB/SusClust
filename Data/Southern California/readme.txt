Southern California Earthquake Catalogs (SCEDC)
==============================================

Source and citation (applies to both files):
  The relocated Southern California catalog of Hauksson et al. (2012) is obtained from SCEDC (2013).
  SCEDC (2013): Southern California Earthquake Center. Caltech.Dataset. doi:10.7909/C3WD3xH1.

1. sc_conv_1-1-1981_31-3-2022.txt
--------------------------------
  Converted catalogue with four whitespace-separated columns.

  Time span:
    1981-01-01 00:00:00 to 2022-03-31 23:59:59

  Region:
    Latitude between 32 and 37 degrees N.
    Longitude between -121 and -114 degrees (east-positive).
    Depth between 0 and 55 km in the source catalogue; depth is omitted here.

  Columns:
    1. Elapsed time (seconds since 1981-01-01 00:00:00)
    2. Latitude (degrees)
    3. Longitude (degrees; negative values denote west)
    4. Magnitude (best available magnitude estimate)

2. sc_1981_31_8_2023.txt
------------------------
  Catalogue with ten whitespace-separated numeric columns and no header.
  Calendar date/time and depth are retained in this file.

  Time span:
    1981-01-01 to 2023-08-31
    First record: 1981-01-01 01:49:29.506
    Last record:  2023-08-31 23:18:53.180

  Columns:
    1. Year
    2. Month
    3. Day
    4. Hour
    5. Minute
    6. Second (including fractional seconds)
    7. Latitude (degrees)
    8. Longitude (degrees; negative values denote west)
    9. Depth (km)
   10. Magnitude

  Example record:
    1981 01 01 01 49 29.506 33.73009 -118.81314 7.126 2.27
