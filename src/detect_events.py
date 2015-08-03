__author__ = 'ubriela'

import SwaggerMediaq as sm
import geojson
from geojson import Feature, Point, FeatureCollection
import requests
import json
import sys


url = "http://mediaq.usc.edu/MediaQ_MVC_V3/api"
validate_endpoint = 'http://geojsonlint.com/validate'

class Sequence(object):
    start, end, min, max, pos_min, pos_max = None, None, None, None, None, None

    def __str__(self):
        return str(self.start) + '\t' + str(self.end) + '\t' + str(self.min) + '\t' + str(self.max) + '\t' + str(self.pos_min) + '\t' + str(self.pos_max)

    def __init__(self, start, end, min, max, pos_min, pos_max):
        self.start = start
        self.end = end
        self.min = min
        self.max = max
        self.pos_min = pos_min
        self.pos_max = pos_max



RANGE = 10
MIN_LENGTH = 50
def events_by_direction(dirs):
    results = []
    for i in range(len(dirs)):
        curr = dirs[i]
        if i == 0:
            sequence = Sequence(0, 0, curr, curr, 0, 0)   # (start, end, min, max, pos_min, pos_max)
        elif sequence.min <= curr <= sequence.max:    # within range
            sequence.end = i         # update end
        elif sequence.max - RANGE <= curr < sequence.min:  # left
            sequence.end, sequence.min, sequence.pos_min = i, curr, i  # update min, pos_min
        elif sequence.max < curr <= sequence.min + RANGE:
            sequence.end, sequence.max, sequence.pos_max = i, curr, i  # update max, pos_max
        elif sequence.min - RANGE < curr < sequence.max - RANGE:
            if sequence.end - sequence.start >= MIN_LENGTH - 1:  # direction event is detected
                results.append(sequence)
            # update max
            pos_max = i
            max = curr
            for j in reversed(range(sequence.start, i)):
                if dirs[j] > max and dirs[j] - curr <= RANGE:
                    max, pos_max = dirs[j], j
                if dirs[j] - curr <= RANGE:
                    continue
                else:
                    break
            sequence = Sequence(j + 1, i, curr, max, i, pos_max)
        elif sequence.min + RANGE < curr < sequence.max + RANGE:
            if sequence.end - sequence.start >= MIN_LENGTH - 1:  # direction event is detected
                results.append(sequence)
            # update min
            pos_min = i
            min = curr
            for j in reversed(range(sequence.start, i)):
                if dirs[j] < min and curr - dirs[j] <= RANGE:
                    min, pos_min = dirs[j], j
                if curr - dirs[j] <= RANGE:
                    continue
                else:
                    break
            sequence = Sequence(j + 1, i, min, curr, pos_min, i)

        elif curr <= sequence.min - RANGE or curr >= sequence.max + RANGE:    # out of range
            if sequence.end - sequence.start >= MIN_LENGTH - 1:  # direction event is detected
                results.append(sequence)
            sequence = Sequence(i, i, curr, curr, i, i)

    return results



# dirs = [120, 150, 160, 140, 150, 145, 135, 130, 150, 120]
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
            results = events_by_direction(dirs)

            for s in results:
                print vid, s
