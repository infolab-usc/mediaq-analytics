__author__ = 'ubriela'

import SwaggerMediaq as sm
import geojson
from geojson import Feature, Point, FeatureCollection
import requests
import json
import sys
import math
import re
import numpy as np
import time
import matplotlib.pyplot as plt
from math import atan2, degrees, pi


sys.path.append('/Users/ubriela/Dropbox/_USC/_Research/_Crowdsourcing/DisasterResponse/BDR/src/bdr/')
from FOV import FOV
from Params import Params
from UtilsBDR import mbr_to_cellids, cell_coord, distance_km


url = "http://mediaq.usc.edu/MediaQ_MVC_V3/api"
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
MIN_LENGTH = 50
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

"""
1.5GB --> 1500
1.5MB --> 1.5
250K --> 0.25
"""
def size_in_mb(size):
    if size:
        value = re.sub("[^0-9.]", "", size)
        unit = re.sub("[^a-zA-z]", "", size)

        if unit == "B":
            return float(value) * 0.000001
        if unit == "K":
            return float(value)  * 0.001
        elif unit == "M":
            return float(value)
        elif unit == "G":
            return float(value) * 1000
        elif unit == "T":
            return float(value) * 1000000
        else:
            return None
    else:
        return None
# print size_in_mb("100.5M")

# Create geoq client api
geoq = sm.GeoqApi(sm.ApiClient(url))

# Returns a set of video locations in GEOJSON format (small sample data)
#geoq.sample_videos()

# Returns a set of video frames (of a particular video) in GEOJSON format (small sample data)
#geoq.sample_fovs()

# Create geoq client with API key
# Replace KEY_VALUE by actual one
geoq = sm.GeoqApi(sm.ApiClient(url, "X-API-KEY", "8b51UFM2SlBltx3s6864eUO1zSoefeK5"))

def getVideos(swlat=34.018212, swlng=-118.291716, nelat=34.025296, nelng=-118.279826):

    # Returns a set of video locations
    # http://mediaq.usc.edu/MediaQ_MVC_V3/api/geoq/rectangle_query?swlat=34.018212&swlng=-118.291716&nelat=34.025296&nelng=-118.279826&X-API-KEY=REAL_KEY&X-API-KEY=8b51UFM2SlBltx3s6864eUO1zSoefeK5

    # Returns a set of video locations that are captured within a time interval (startdate -> enddate)
    # swlat=34.018212, swlng=-118.291716,nelat=34.025296, nelng=-118.279826,startdate="2014-04-13 00:00:00",enddate="2014-04-13 23:59:59"
    # swlat=34.019972, swlng=-118.291588, nelat=34.021111, nelng=-118.287125
    fc_videos = geoq.rectangle_query(swlat, swlng, nelat, nelng, startdate="2014-04-12 00:00:00",enddate="2014-04-13 23:59:59")
    fc_videos = fc_videos.replace('None','null').replace('u\'','\"').replace('\'','\"')

    # print fc_videos

    # validate GEOJSON
    geojson_validation = requests.post(validate_endpoint, data=fc_videos)
    if geojson_validation.json()['status'] != 'ok':
        print "Rectangle_query: Invalid geojson format"
        exit()

    fc_videos = geojson.loads(fc_videos)
    print "Number of videos: " + str(len(fc_videos.features))
    fov_counts = []
    video_sizes = []
    valid_fc_videos = [] # with video of size > 0
    for video in fc_videos.features:
        fov_counts.append(video.properties['fov_count'])
        size = size_in_mb(video.properties['size'])
        if size and size > 0:
            valid_fc_videos.append(video)
            video_sizes.append(size)

    print "Total/Mean video size (M): " + str(sum(video_sizes)) + "/" + str(np.mean(video_sizes))
    bins = range(0,100,10)
    hist_video_counts = np.histogram(video_sizes, bins)[0]
    print "Histogram:" + str(hist_video_counts)

    # the histogram of the data
    plt.hist(video_sizes, bins=bins)
    plt.xlabel('Range (MB)')
    plt.ylabel('Video Count')
    plt.title(r'$\mathrm{Histogram\ of\ Video\ Size}$')
    # plt.axis([0, 100, 0, 500])
    plt.grid(True)
    # plt.show()

    print "Total/Mean FOV count: " + str(sum(fov_counts)) + "/" + str(np.mean(fov_counts))
    bins = range(0,200,20)
    hist_fov_counts = np.histogram(fov_counts, bins)[0]
    print "Histogram:" + str(hist_fov_counts)

    # the histogram of the data
    plt.hist(fov_counts, bins=bins)
    plt.xlabel('Range (Number of Frames)')
    plt.ylabel('Frame Count')
    plt.title(r'$\mathrm{Histogram\ of\ Frame\ Count}$')
    # plt.axis([0, 100, 0, 500])
    plt.grid(True)
    # plt.show()

    return valid_fc_videos

def getFOVs(vid):
    fovs = []
    # Returns a set of video frames
    try:
        fovs = geoq.video_metadata(vid).replace('None','null').replace('u\'','\"').replace('\'','\"')
        geojson_validation = requests.post(validate_endpoint, data=fovs)
        if geojson_validation.json()['status'] != 'ok':
            print "Video_metadata: Invalid geojson format"
        fovs = geojson.loads(fovs)
    except Exception as inst:
        print vid
        print inst
        print "Unexpected error:", sys.exc_info()[0]

    return fovs

def detect_events(video):
    if video.properties['vid']:
        vid = str(video.properties['vid'])
        videoid = str(video.properties['videoid'])
        fovs = getFOVs(vid)

        if fovs and len(fovs) > 0:
            # print vid
            trajectory = [(fov.geometry.coordinates[0],fov.geometry.coordinates[1],float(fov.properties['theta_x'])) for fov in fovs.features]
            dirs = [dir[2] for dir in trajectory]
            locs = [(dir[0], dir[1]) for dir in trajectory]
            #print dirs
            results = events_by_direction(dirs)
            for s in results:
                print len(dirs), videoid, s, dirs[s.start:s.end]
                print geoq.video_segment_url(vid, s.start, s.end)


            # results = events_by_locality(locs)
            # for s in results:
            #     print s
            #     print geoq.video_segment_url(vid, s.start, s.end)


def test_detect_events(fc_videos):
    for video in fc_videos:
        detect_events(video)


def dump_metadata_dataset(filename = "mediaq_fovs.txt"):

    file = open(filename,"w")

    videos = getVideos()
    videoidx = 1
    for video in videos:
        # vid,videoid,fovnum,Plat,Plng,Px,Py,prevX,prevY,speed,dir,prevDir,R,alpha,timestamp
        vid = str(video.properties['vid'])

        fovs = getFOVs(vid)
        # print vid, fovs
        if fovs and len(fovs) > 0:
            for fov in fovs.features:
                plat = fov.geometry.coordinates[1]
                plng = fov.geometry.coordinates[0]
                dir = float(fov.properties['theta_x'])
                R = float(fov.properties['r']) * 200
                alpha = float(fov.properties['alpha'])
                fov_num = int(fov.properties['fov_num'])
                timestamp = None
                line = " ".join(map(str, [vid, videoidx, fov_num, plat, plng, None, None, None, None, None, dir, None, R, alpha, timestamp]))
                file.write(line + "\n")
        videoidx = videoidx + 1

    file.close()

def compute_coverage_map(grid_size = 200):
    swlat=34.018212
    swlng=-118.291716
    nelat=34.025296
    nelng=-118.279826
    videos = getVideos(swlat, swlng, nelat, nelng)

    map = np.ndarray(shape=(grid_size, grid_size), dtype=int)
    for video in videos:
        if video.properties['vid']:
            vid = str(video.properties['vid'])
            fovs = getFOVs(vid)

            if fovs and len(fovs) > 0:
                for fov in fovs.features:
                    f = FOV(fov)
                    param = Params(1000, swlat, swlng, nelat, nelng)
                    param.GRID_SIZE = grid_size
                    for cid in f.cellids(param):
                        cell_lat, cell_lng = cell_coord(cid, param)
                        if f.cover(cell_lat, cell_lng):
                            y_idx = cid/param.GRID_SIZE
                            x_idx = cid - y_idx*param.GRID_SIZE
                            # print x_idx, y_idx, map[x_idx][y_idx]
                            map[x_idx][y_idx] = map[x_idx][y_idx] + 1


    fig, ax = plt.subplots()
    heatmap = ax.pcolor(map, cmap=plt.cm.Blues)
    plt.show()
    np.savetxt("mediaq_coverage_heatmap.txt" , map, fmt='%i\t')




# dump_metadata_dataset()

# start = time.time()
videos = getVideos()
# print time.time() - start

for video in videos:
    detect_events(video)

# compute_coverage_map()


