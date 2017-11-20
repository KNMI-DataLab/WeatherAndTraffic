# sudo apt-get install python-pyshp
# http://opendata.ndw.nu/NDW_Shapefiles_20170926.zip
from audioop import reverse

import shapefile
import pandas
from scipy import stats
import datetime
import netCDF4
import numpy as np
import dateutil.parser
import time
from collections import OrderedDict, defaultdict
import itertools


def getSpeed(timeandtype, location, time):
    rows = timeandtype.loc[timeandtype.periodStart == time]

    indices = rows['index'].unique()
    firstIndex = indices[0]

    values = []
    modulo = -1
    if firstIndex.find('A') != -1:  # No category for speed: 2A, 4A, 6A and 8A
        modulo = 2
    if firstIndex.find('B') != -1:  # Three vehicle categories for speed: 8B, 16B, 24B and 32B
        modulo = 8
    if firstIndex.find('C') != -1:  # Six vehicle categories for speed: 12C, 24C, 36C and 48C
        modulo = 12

    for index in indices:
        value = int(index[:-1])
        if value < 200 and value % modulo == 0:
            speed = float(rows.loc[(rows['index'] == index)].avgVehicleSpeed)
            # print speed
            if speed >= 0:
                values.append(speed)
    return stats.hmean(values)


def getIntensity(timeandtype, location, time):
    rows = timeandtype.loc[(timeandtype.periodStart == time)]
    indices = rows['index'].unique()
    firstIndex = indices[0]

    values = []
    modulo = -1
    offset = 0
    if firstIndex.find('A') != -1:  # No category for flow: 1A, 3A, 5A and 7A
        modulo = 2
        offset = 1
    if firstIndex.find('B') != -1:  # Three vehicle categories for flow: 4B, 12B, 20B and 28B
        modulo = 8
        offset = 4
    if firstIndex.find('C') != -1:  # Six vehicle categories for flow: 6C, 18C, 30C, 42 and 54C
        modulo = 12
        offset = 6

    for index in indices:
        value = int(index[:-1])
        if value < 200 and (value + offset) % modulo == 0:
            flow = float(rows.loc[(rows['index'] == index)].avgVehicleFlow)
            # print flow
            if flow >= 0:
                values.append(flow)
    return sum(values)


def getRoadId(table,location,reverseDictRoad):
    rows = table.loc[table.measurementSiteReference == location]
    #locationID = \
    try:
        roadID = reverseDictRoad[str(rows.iloc[0].ROADNUMBER)+"-"+str(rows.iloc[0].measurementSide)] #.apply(lambda row: str(row.ROADNUMBER)+"-"+str(row.measurementSide), axis=1)
    except KeyError:
        roadID = 9999
    return roadID





# Make a map of locations with their coords
sf = shapefile.Reader("Telpunten_20170926.shp")
records = sf.records()
shapes = sf.shapes()
locationsByName = {}
for j in range(0, len(records)):
    name = records[j][2]
    lon = shapes[j].points[0][0]
    lat = shapes[j].points[0][1]
    locationsByName[name] = [lat, lon]

#Read in the files with data and create corresponing pandas tables
flowTable = pandas.read_csv('dataNDW/accident20151230_intensiteit_00001.csv')
speedTable = pandas.read_csv('dataNDW/accident20151230_snelheid_00001.csv')
#print getIntensity(flowTable, 'RWS01_MONIBAS_0020vwm0558ra', '2015-12-31 23:45:00');
#print getSpeed(speedTable, 'RWS01_MONIBAS_0020vwm0558ra', '2015-12-31 23:45:00');

speedTable.sort_values(['startLocatieForDisplayLat','startLocatieForDisplayLong'])

locations = speedTable.sort_values(['measurementSide','startLocatieForDisplayLat','startLocatieForDisplayLong']).measurementSiteReference.unique()
print(locations)
dates = speedTable.periodStart.unique()

#create dictionary with mapping road number and side (i.e., north/south bound) and flag
roadsTab = speedTable.drop_duplicates(['ROADNUMBER', 'measurementSide'])
roadsTab = roadsTab[['ROADNUMBER','measurementSide']].dropna()
roadsTab = roadsTab.apply(lambda x: '-'.join(x), axis=1)
roadsTab=roadsTab.reset_index()
roadsTab=roadsTab.drop("index",1)
seriesRoad=roadsTab.ix[:,0]
roadIdDict=seriesRoad.to_dict()
#print len(locations)
#print len(dates)
reverseDict = {v: k for k, v in roadIdDict.iteritems()}


##netcdf part initialization
fileOutName = "ndw_speed_intensityA2.nc"
latvar = []
lonvar = []
timevardata = []
roadIdvarData = []
flowvar = []
speedvar = []
station = np.array([], dtype='object')
###

#extract locations and filter the ones available in the measurements and the current measurement infrastructure (shape file)
locationsCorrected = locations




for loc in locations:
    timeandtypespeed = speedTable.loc[(speedTable.measurementSiteReference == loc)]
    try:
        latTemp = locationsByName[loc][0]
        lonTemp = locationsByName[loc][1]
    except KeyError, e:
        print 'Error: old station no more available for measurements', e
        #removing locations that are in the data set but NOT in the current shapefile (old locations that have been
        # removed from the current measurement locations)
        index = np.argwhere(locationsCorrected == loc)
        locationsCorrected = np.delete(locationsCorrected, index)
        continue
    roadIdTemp = getRoadId(timeandtypespeed,loc,reverseDict)
    if (roadIdTemp !=9999):
        latvar.append(latTemp)
        lonvar.append(lonTemp)
        roadIdvarData.append(roadIdTemp)
        station = np.append(station, loc)
    else:
        index = np.argwhere(locationsCorrected == loc)
        locationsCorrected = np.delete(locationsCorrected, index)

#extract flow and speed information and time and put into arrays to fill netcdf variables
for timeVariable in dates:  # Note localtime used , not UTC
    inputdate = timeVariable + 'CEST'
    date = dateutil.parser.parse(inputdate)  # ,"%Y-%m-%d %M:%H:%S")
    utctimestring = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.mktime(date.timetuple())))
    netcdfDate = netCDF4.date2num(dateutil.parser.parse(utctimestring), "seconds since 1970-01-01 00:00:00")
    timevardata.append(netcdfDate)
    for loc in locationsCorrected:
    #loc = 'RWS01_MONIBAS_0020vwm0558ra'
        timeandtypespeed = speedTable.loc[(speedTable.measurementSiteReference == loc)]
        timeandtypeflow = flowTable.loc[(flowTable.measurementSiteReference == loc)]
        speedvar.append(getSpeed(timeandtypespeed,loc,timeVariable))
        flowvar.append(getIntensity(timeandtypeflow, loc, timeVariable))

#create netcdf variables
numpoints = len(station)
timeSteps = len(timevardata)

ncfile = netCDF4.Dataset(fileOutName, 'w')
obs_dim = ncfile.createDimension('station', numpoints)  # latitude axis
time_dim = ncfile.createDimension('time', timeSteps)

lat = ncfile.createVariable('lat', 'd', ('station'))
lat.units = 'degrees_north'
lat.standard_name = 'latitude'
lon = ncfile.createVariable('lon', 'd', ('station'))
lon.units = 'degrees_east'
lon.standard_name = 'longitude'

timevar = ncfile.createVariable('time', 'd', ('time'))
timevar.units = "seconds since 1970-01-01 00:00:00"
timevar.standard_name = 'time'

roadIdvar = ncfile.createVariable('roadId', 'u2', ('station'))
roadIdMeaningsList = roadIdDict.values()
roadIdvar.flags_meanings = ', '.join(roadIdMeaningsList)
roadIdvar.flags_values = roadIdDict.keys()
minFlag = min(roadIdDict.keys())
maxFlag = max(roadIdDict.keys())
roadIdvar.valid_range = [minFlag, maxFlag]

floatVar = ncfile.createVariable('station', str, ('station'))
floatVar.units = '-'
floatVar.standard_name = 'station'

flowVar = ncfile.createVariable('flow', 'd', ('station', 'time'))


speedVar = ncfile.createVariable('speed', 'd', ('station', 'time'))

#fill netcdf variables
lat[:] = [latvar]
lon[:] = [lonvar]
timevar[:] = [timevardata]
floatVar[:] = station
roadIdvar[:] = [roadIdvarData]
flowVar[:] = [flowvar]
speedVar[:] = [speedvar]

#netcdf stuff for header
ncfile.featureType = "timeSeries";
ncfile.Conventions = "CF-1.4";
ncfile.close()


#Example netcdf file similar to the one written
"""
netcdf Amsterdam20170314 {
dimensions:
        station = 1110 ;
        time = 360 ;
        distance = 1110 ;
variables:
        double lat(station) ;
                lat:units = "degrees_north" ;
        double lon(station) ;
                lon:units = "degrees_east" ;
        int station(station) ;
        double time(time) ;
                time:units = "seconds since 1970-01-01 00:00:00" ;
                time:standard_name = "time" ;
        float flow(station, time) ;
        float speed(station, time) ;
        float distance(distance) ;
        short roadId(station) ;
                roadId:flag_meanings = "A1R A1L A10L A10R A4L A4R A9L A2R A2L " ;
                roadId:flag_values = 0s, 1s, 2s, 3s, 4s, 5s, 6s, 7s, 8s ;
                roadId:valid_range = 0s, 9s ;
// global attributes:
                :featureType = "timeSeries" ;
}
"""














