__author__ = 'ubriela'

import SwaggerMediaq as sm
import geojson
from geojson import Feature, Point, FeatureCollection
import requests
import json
import sys
import math


url = "http://local.eclipse.org/MediaQ_MVC_V3/api"
validate_endpoint = 'http://geojsonlint.com/validate'

# An abstract event
class Event(object):
    start, end = None, None

    def __str__(self):
        return str(self.start) + '\t' + str(self.end)

    def __init__(self, start, end):
        self.start = start
        self.end = end


# A directional event
class DirectionalEvent(Event):
    min, max, pos_min, pos_max = None, None, None, None

    def __str__(self):
        return str(Event.__str__(self)) + '\t' + str(self.min) + '\t' + str(self.max) + '\t' + str(self.pos_min) + '\t' + str(self.pos_max)

    def __init__(self, start, end, min, max, pos_min, pos_max):
        Event.__init__(self, start, end)
        self.min = min
        self.max = max
        self.pos_min = pos_min
        self.pos_max = pos_max

class SpatialEvent(Event):

    def __init__(self, start, end, min, max, pos_min, pos_max):
        Event.__init__(self, start, end)

def distance_km(lat1, lon1, lat2, lon2):
    """
    Distance between two geographical locations
    """
    R = 6371  # km
    dLat = math.radians(abs(lat2 - lat1))
    dLon = math.radians(abs(lon2 - lon1))
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)

    a = math.sin(dLat / 2) * math.sin(dLat / 2) + math.sin(dLon / 2) * math.sin(dLon / 2) * math.cos(lat1) * math.cos(
        lat2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = R * c
    return d

SPATIAL_DISTANCE = 0.001    # km
SPATIAL_MIN_LENGTH = 20

def events_by_locality(dirs):
    results = []
    for i in range(len(dirs)):
        curr = dirs[i]
        if i == 0:
            event = Event(0, 0)   # (start, end)
        elif distance_km(dirs[event.start][0], dirs[event.start][1], curr[0], curr[1]) <= SPATIAL_DISTANCE:    # within range
            event.end = i         # update end
        else:   # new range
            if event.end - event.start >= SPATIAL_MIN_LENGTH:
                results.append(event)
            # new event
            event = Event(i, i)
    return results


DIRECTIONAL_RANGE = 10
MIN_LENGTH = 10
def events_by_direction(dirs):
    results = []
    for i in range(len(dirs)):
        curr = dirs[i]
        if i == 0:
            sequence = DirectionalEvent(0, 0, curr, curr, 0, 0)   # (start, end, min, max, pos_min, pos_max)
        elif sequence.min <= curr <= sequence.max:    # within range
            sequence.end = i         # update end
        elif sequence.max - DIRECTIONAL_RANGE <= curr < sequence.min:  # left
            sequence.end, sequence.min, sequence.pos_min = i, curr, i  # update min, pos_min
        elif sequence.max < curr <= sequence.min + DIRECTIONAL_RANGE:
            sequence.end, sequence.max, sequence.pos_max = i, curr, i  # update max, pos_max
        elif sequence.min - DIRECTIONAL_RANGE < curr < sequence.max - DIRECTIONAL_RANGE:
            if sequence.end - sequence.start >= MIN_LENGTH - 1:  # direction event is detected
                results.append(sequence)
            # update max
            pos_max = i
            max = curr
            for j in reversed(range(sequence.start, i)):
                if dirs[j] > max and dirs[j] - curr <= DIRECTIONAL_RANGE:
                    max, pos_max = dirs[j], j
                if dirs[j] - curr <= DIRECTIONAL_RANGE:
                    continue
                else:
                    break
            sequence = DirectionalEvent(j + 1, i, curr, max, i, pos_max)
        elif sequence.min + DIRECTIONAL_RANGE < curr < sequence.max + DIRECTIONAL_RANGE:
            if sequence.end - sequence.start >= MIN_LENGTH - 1:  # direction event is detected
                results.append(sequence)
            # update min
            pos_min = i
            min = curr
            for j in reversed(range(sequence.start, i)):
                if dirs[j] < min and curr - dirs[j] <= DIRECTIONAL_RANGE:
                    min, pos_min = dirs[j], j
                if curr - dirs[j] <= DIRECTIONAL_RANGE:
                    continue
                else:
                    break
            sequence = DirectionalEvent(j + 1, i, min, curr, pos_min, i)

        elif curr <= sequence.min - DIRECTIONAL_RANGE or curr >= sequence.max + DIRECTIONAL_RANGE:    # out of range
            if sequence.end - sequence.start >= MIN_LENGTH - 1:  # direction event is detected
                results.append(sequence)
            sequence = DirectionalEvent(i, i, curr, curr, i, i)

    return results



# dirs = [120, 150, 160, 140, 150, 145, 135, 130, 150, 120]
# dirs = [181.94008, 182.37218, 171.61913, 160.90152, 157.55756, 153.66913, 156.21973, 175.0801, 159.62263, 153.34453, 151.51772, 151.01799, 152.5305, 155.17336, 154.8254, 150.19846, 166.72073, 158.32208, 148.25983, 154.22003, 162.4578, 167.7894, 181.1051, 170.45929, 175.29805, 175.39275, 186.90277, 188.76363, 200.17162, 200.49422, 194.32492, 197.19994, 202.31071, 204.23193, 197.63225, 198.6546, 205.44516, 198.59033, 200.10536, 200.5626, 202.47163, 204.49399, 204.32697, 205.05452, 209.90833, 210.01016, 220.22246, 240.85333, 247.69794, 254.44272, 265.39795, 270.58826, 274.80835, 270.61554, 266.99905, 194.38068, 161.73622, 168.23277, 176.72787, 187.61053, 204.8588]
# print dirs[31:43]
# results = events_by_direction(dirs)
# for s in results:
#     print s



# Create geoq client api
geoq = sm.GeoqApi(sm.ApiClient(url))

# Returns a set of video locations in GEOJSON format (small sample data)
#geoq.sample_videos()

# Returns a set of video frames (of a particular video) in GEOJSON format (small sample data)
#geoq.sample_fovs()

# Create geoq client with API key
# Replace KEY_VALUE by actual one
geoq = sm.GeoqApi(sm.ApiClient(url, "X-API-KEY", "8b51UFM2SlBltx3s6864eUO1zSoefeK5"))

# Returns a set of video locations
#geoq.rectangle_query(swlat=34.019972,swlng=-118.291588,nelat=34.021111,nelng=-118.287125)

# Returns a set of video locations that are captured within a time interval (startdate -> enddate)
# swlat=34.018212, swlng=-118.291716,nelat=34.025296, nelng=-118.279826,startdate="2014-04-13 00:00:00",enddate="2014-04-13 23:59:59"
# swlat=34.019972, swlng=-118.291588, nelat=34.021111, nelng=-118.287125
fc_videos = geoq.rectangle_query(swlat=34.018212, swlng=-118.291716, nelat=34.025296, nelng=-118.279826, startdate="2014-04-12 10:00:00",enddate="2014-04-13 23:59:59")
fc_videos = fc_videos.replace('u\'','\"').replace('\'','\"')

# validate GEOJSON
geojson_validation = requests.post(validate_endpoint, data=fc_videos)
if geojson_validation.json()['status'] != 'ok':
    print "Rectangle_query: Invalid geojson format"
    exit()

fc_videos = geojson.loads(fc_videos)
print "Number of videos: " + str(len(fc_videos.features))
for video in fc_videos.features:
    if video.properties['vid']:
        vid = str(video.properties['vid'])
	videoid = str(video.properties['videoid'])
        fovs = []
        # Returns a set of video frames
        try:
            fovs = geoq.video_metadata(vid).replace('u\'','\"').replace('\'','\"')
            geojson_validation = requests.post(validate_endpoint, data=fovs)
            if geojson_validation.json()['status'] != 'ok':
                print "Video_metadata: Invalid geojson format"
            fovs = geojson.loads(fovs)
        except Exception as inst:
            print vid
            print inst
            print "Unexpected error:", sys.exc_info()[0]

        if fovs:
            # print vid
            trajectory = [(fov.geometry.coordinates[0],fov.geometry.coordinates[0],float(fov.properties['theta_x'])) for fov in fovs.features]
            dirs = [dir[2] for dir in trajectory]
            locs = [(dir[0], dir[1]) for dir in trajectory]
	        #print dirs
            # results = events_by_direction(dirs)
            # for s in results:
            #     print len(dirs), s
                # print geoq.video_segment_url(vid, s.start, s.end)


            results = events_by_locality(locs)
            for s in results:
                print s
                print geoq.video_segment_url(vid, s.start, s.end)
