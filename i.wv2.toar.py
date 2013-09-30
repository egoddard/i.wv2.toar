#!/usr/bin/env python

##############################################################################
#
# MODULE:	    i.wv2.toar
#
# AUTHOR(S):    Eric Goddard, egoddard@memphis.edu
#
# PURPOSE:	    Converts WorldView-2 pixel counts to Top of Atmosphere
#               Reflectance values. THIS HAS ONLY BEEN TESTED ON WV-2
#               PRODUCTLEVEL LV3D.
#
# COPYRIGHT:	(C) 2013 Eric Goddard
#
# LICENSE:      GPLv2
#
# REFERENCES:
#       2010  Updike Todd and Chris Comp. Radiometric Use of WorldView-2
#                  Imagery. DigitalGlobe. Electronic Document,
#                  http://www.digitalglobe.com/sites/default/files/Radiometric_Use_of_WorldView-2_Imagery%20%281%29.pdf.
#                  Accessed 12 July 2013.
##############################################################################

#%module
#%  description: Converts WorldView-2 counts to Top of Atmosphere Reflectances
#%  keywords: WorldView-2
#%  keywords: reflectance
#%  keywords: radiance
#%  keywords: atmosphere
#%  keywords: remote sensing
#%  keywords: radiometric correction
#%end
#%option G_OPT_R_INPUT
#%  key: map
#%  description: Select the WorldView-2 raster to be corrected.
#%  required: yes
#%end
#%option
#%  key: output_prefix
#%  description: Prefix for output reflectance maps.
#%  required:yes
#%end
#%option G_OPT_F_INPUT
#%  key: metadata_xml
#%  description: Path to the XML metadata file associated with the input map.
#%  required: yes
#%end
#%option
#%  key:band
#%  description: Choose the band to be corrected.
#%  required: yes
#%  options: Pan, Coastal, Blue, Green, Yellow, Red, Red Edge, NIR1, NIR2
#%  type: string
#%end
#%flag
#%  key: n
#%  description: Use band numbers instead of letters.
#%end
#%flag
#%  key: r
#%  description: Output radiance instead of reflectance

import sys
import os
from datetime import datetime
import math
import grass.script as grass
try:
    import xml.etree.cElementTree as etree
except ImportError:
    try:
        import xml.etree.ElementTree as etree
    except ImportError:
        try:
            import elementtree.ElementTree as etree
        except ImportError:
            grass.message("Error: unable to import ElementTree module to parse the \
                          XML metadata file.")
            sys.exit(1)


if "GISBASE" not in os.environ:
    grass.message("You bust be in GRASS GIS to run this program.")
    sys.exit(1)

def EarthSunDistance(date):
    """Take a string date as input and converts it to the Julian day and
    returns the Earth-Sun Distance for that day."""

    acqTime = datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%fZ")
    if acqTime.month == 1 or acqTime.month == 2:
        month = acqTime.month + 12
        year = acqTime.year - 1
    else:
        month = acqTime.month
        year = acqTime.year

    secondsMicro = acqTime.second + (acqTime.microsecond / 1000000.0)
    UT = acqTime.hour + (acqTime.minute / 60.0) + (secondsMicro / 3600.0)
    A = int(year / 100)
    B = 2 - A + int(A / 4)
    JD = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + \
            acqTime.day + (UT / 24.0) + B - 1524.5
    D = JD - 2451545.0
    g = math.radians(357.529 + 0.98560028 * D)

    return 1.00014 - 0.01671 * math.cos(g) - 0.00014 * math.cos(2 * g)


def ExtractVariablesFromMetadata(metadataFilePath, band):
    with open(metadataFilePath, 'rb') as metadataFile:
        root = etree.parse(metadataFile).getroot()

        try:
            absCalFactor = root.find("IMD/%s/ABSCALFACTOR" % band).text
            effectiveBandWidth = root.find("IMD/%s/EFFECTIVEBANDWIDTH" % band).text
            sunElevation = root.find("IMD/IMAGE/MEANSUNEL").text
            acqTime = root.find("IMD/MAP_PROJECTED_PRODUCT/EARLIESTACQTIME").text
            earthSunDist = EarthSunDistance(acqTime)

            #WV-2 Band-Averaged Solar Spectral Irradiance Table
            SSI = {"BAND_P": 1580.8140,
                   "BAND_C": 1758.2229,
                   "BAND_B": 1974.2416,
                   "BAND_G": 1856.4104,
                   "BAND_Y": 1738.4791,
                   "BAND_R": 1559.4555,
                   "BAND_RE": 1342.0695,
                   "BAND_N": 1069.7302,
                   "BAND_N2": 861.2866}

            return [float(absCalFactor), 90 - float(sunElevation),
                    float(earthSunDist), SSI[band], float(effectiveBandWidth)]

        except:
            grass.message("Error: Could not find required variables in the \
                          metadata file.")
            sys.exit(1)


def CalculateTOAR(rastItems, output, radiance):
    raster = rastItems[0]
    absCalFactor = rastItems[1]
    theta = rastItems[2]
    earthSunDist = rastItems[3]
    eSun = rastItems[4]
    effectiveBandwidth = rastItems[5]

    grass.use_temp_region()
    grass.run_command('g.region', rast=raster)

    if radiance:
        calc = "$output = ($absCal * $input_rast) / $efb"
        grass.message("Calculating Top of Atmosphere Radiance for {}".format(raster))
        grass.mapcalc(calc, output=output, absCal=absCalFactor, input_rast=raster,
                      efb=effectiveBandwidth)

    else:
        calc = "$output = ((($absCal * $input_rast) / $efb) * $esd^2 * $pi) / ($eSun * cos($theta))"
        grass.message("Calculating Top of Atmosphere Reflectance for %s" % raster)
        grass.mapcalc(calc, output=output, absCal=absCalFactor, input_rast=raster,
                      esd=earthSunDist, pi=math.pi, eSun=eSun, theta=theta,
                      efb=effectiveBandwidth)


def main():
    input_map = options['map']
    output_prefix = options['output_prefix']
    metadata_path = options['metadata_xml']
    band = options['band']
    use_band_names = flags['n']
    radiance = flags['r']
    #serialized = flags['s']

    #create lookup tables for band options to metadata band designations
    band_lookup = {'Pan': 'BAND_P',
                  'Coastal': 'BAND_C',
                  'Blue': 'BAND_B',
                  'Green': 'BAND_G',
                  'Yellow': 'BAND_Y',
                  'Red': 'BAND_R',
                  'Red Edge': 'BAND_RE',
                  'NIR1': 'BAND_N',
                  'NIR2': 'BAND_N2'}

    band_num = {'Pan': 'P',
               'Coastal': '1',
               'Blue': '2',
               'Green': '3',
               'Yellow': '4',
               'Red': '5',
               'Red Edge': '6',
               'NIR1': '7',
               'NIR2': '8'}

    if use_band_names:
        output = "%s.%s" % (output_prefix, band_lookup[band])
    else:
        output = "%s.%s" % (output_prefix, band_num[band])

    #Extract conversion attributes for input band.
    variables = [input_map] + ExtractVariablesFromMetadata(metadata_path, bandLookup[band])
    CalculateTOAR(variables, output, radiance)

if __name__ == '__main__':
    options, flags = grass.parser()
    main()
