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
import itertools
import csv
from collections import defaultdict
from math import atan2, degrees, pi

sys.path.append('/Users/ubriela/git/BDR/src/bdr/')
sys.path.append('/Users/ubriela/git/BDR/src/plot/code')
sys.path.append('/Users/ubriela/git/privategeocrowddynamic/src/common')
sys.path.append('/Users/ubriela/git/privategeocrowddynamic/src/minball')

from FOV import FOV
from Video import Video
from Params import Params
from UtilsBDR import mbr_to_cellids, cell_coord, distance_km
from smallestenclosingcircle import make_circle

#url = "http://mediaq.usc.edu/MediaQ_MVC_V3/api"
url = "http://mediaq1.cloudapp.net/MediaQ_MVC_V3/api"
validate_endpoint = 'http://geojsonlint.com/validate'

VIDEO_FILE = "FOVMetadata.txt"
IS_FROM_MEDIAQ = True # true --> obtain from mediaq; otherwise, retrieve from file

MAX_RADIUS = 0.001    # km

ONE_KM = 0.0089982311916  # convert km to degree

# An abstract scene
class Scene(object):
    start, end = None, None

    def __str__(self):
        return str(self.start) + '\t' + str(self.end)

    def __init__(self, start, end):
        self.start = start
        self.end = end


# A track scene
class TrackScene(Scene):
    min, max, min_fov, max_fov = None, None, None, None

    def __str__(self):
        return str(Scene.__str__(self)) + '\t' + str(self.min) + '\t' + str(self.max) + '\t' + str(self.min_fov) + '\t' + str(self.max_fov)

    def __init__(self, start, end, min, max, min_fov, max_fov):
        Scene.__init__(self, start, end)
        self.min = min
        self.max = max
        self.min_fov = min_fov
        self.max_fov = max_fov

# A arch scene
class ArchScene(Scene):

    def __str__(self):
        return str(Scene.__str__(self))

    def __init__(self, start, end):
        Scene.__init__(self, start, end)

# A zooming scene
class ZoomScene(Scene):

    def __str__(self):
        return str(Scene.__str__(self))

    def __init__(self, start, end):
        Scene.__init__(self, start, end)

# A panning scene
class PanScene(Scene):

    def __str__(self):
        return str(Scene.__str__(self))

    def __init__(self, start, end):
        Scene.__init__(self, start, end)


def within(low, high, mid):
    if high < low:  # swap if high < low
        tmp = low
        low = high
        high = tmp
    if high - low < 180:
        if low <= mid <= high:
            return True
    else:
        if mid >= high or mid <= low:
            return True

    return False

def less_than(low, high):
    if abs(high - low) < 180:
        return low < high
    else:
        return low > high

def within_range(low, high, range):
    if abs(high - low) < 180:
        # if high < low:
        #     print "high shough be larger than low", low, high, range, high - low
        return abs(high - low) <= range
    else:
        # if high > low:
        #     print "high shough be smaller than low", low, high, range, high - low
        return 360 - abs(low - high) <= range

"""
correct the angle
"""
def angle(theta):
    if 0 <= theta <= 360:
        return theta
    elif theta < 0:
        return 360 + theta
    else:
        return theta - 360

# print within(90, 170, 100)
# print within(170, 90, 100)
# print within(90, 271, 10)
# print within(90, 100, 100)
# print less_than(70, 80)
# print within_range(50, 100, 50)
# print within_range(350, 10, 15)


"""
correct mediaq compass value to typical range
"""
def correct_dir(dir):
    if dir <= 90:
        return 90 - dir
    else:
        return 360 - (dir - 90)

# print correct_dir(270)

"""
calculate the end point of a vector

(x2-x1)= r cos
(y2-y1)= r sin

lat/y: 0.00088
lng/x: 0.00107

"""
R_100 = 1
def vector_end_point(lat1, lon1, dir):
    c_dir = correct_dir(dir)
    lat2 = lat1 + 0.00088 * math.sin(math.radians(c_dir)) * R_100
    lon2 = lon1 + 0.00107 * math.cos(math.radians(c_dir)) * R_100
    return (lat2, lon2)


# end = vector_end_point(-118.291561, 34.018326, 30)
# print end
# print end[1] - 34.019106, end[0] + 118.291021
# print 34.019106 - 34.018326, -118.291021 + 118.291561

# 34.019106, -118.291021


# def line(l1, l2):
#     A = (l1[1] - l2[1])
#     B = (l2[0] - l1[0])
#     C = (l1[0]*l2[1] - l2[0]*l1[1])
#     return A, B, -C
# def intersection(L1, L2):
#     D  = L1[0] * L2[1] - L1[1] * L2[0]
#     Dx = L1[2] * L2[1] - L1[1] * L2[2]
#     Dy = L1[0] * L2[2] - L1[2] * L2[0]
#     if D != 0:
#         x = Dx / D
#         y = Dy / D
#         return x,y
#     else:
#         return False


def line(l1, l2):
    A = (l1[0] - l2[0])
    B = (l2[1] - l1[1])
    C = (l1[1]*l2[0] - l2[1]*l1[0])
    return A, B, -C



"""
find intersection point of two vectors
"""
def interesection_point(l1, l2, d1, d2):
    end1 = vector_end_point(l1[0], l1[1], d1)
    end2 = vector_end_point(l2[0], l2[1], d2)

    L1 = line(l1, end1)
    L2 = line(l2, end2)

    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return y, x
    else:
        return False

# locs=[(34.018349, -118.291559),(34.018918, -118.291502),(34.018349, -118.290193), (34.018391, -118.291175)]
# dirs=[60,120,300,0]
# print interesection_point((34.018326, -118.291561), (34.020060, -118.291486), 90, 210)


def angle_bwn_points(l1, l2):
    dx = l2[1] - l1[1]
    dy = l2[0] - l1[0]
    rads = math.atan2(dy,dx)
    rads %= 2*pi
    return math.degrees(rads)

# print angle_bwn_points((34.018364, -118.291552),(34.016487, -118.291037))

# def angle_bwn_vector_dir(l1, l2, dir):
    alpha1 = angle_bwn_points(l1, l2)
    alpha2 = correct_dir(dir)



def is_in_clockwise_range(dir, start, end):
    if start <= end:
        return start <= dir <= end
    else:   #  start > end
        return (start <= dir <= 360) or (0 <= dir <= end)

# print is_in_clockwise_range(45, 0, 90)
# print is_in_clockwise_range(120, 270, 90)

def combine_clockwise_ranges(start1, end1, start2, end2):
    # if start1, end1
    if is_in_clockwise_range(start2, start1, end1) and is_in_clockwise_range(end2, start1, end1):
        return start1, end1
    elif is_in_clockwise_range(end1, start1, end2) and is_in_clockwise_range(start2, start1, end2):
        return start1, end2
    elif is_in_clockwise_range(start1, start2, end1) and is_in_clockwise_range(end2, start2, end1):
        return start2, end1
    elif is_in_clockwise_range(start1, start2, end2) and is_in_clockwise_range(end1, start2, end2):
        return start2, end2
    else:
        print "combine_clockwise_ranges ", start1, end1, start2, end2

# print combine_clockwise_ranges(270, 45, 0, 90)
# print combine_clockwise_ranges(270, 45, 180, 300)
# print combine_clockwise_ranges(270, 0, 180, 45)
# print combine_clockwise_ranges(0, 180, 45, 270)

def angle_clockwise(start, end):
    if start <= end:
        return end - start + 1
    else:
        return 360 - start + end + 1

# print angle_clockwise(0, 0)
# print angle_clockwise(45, 45)

"""
Assumption: two consecutive camera directions forms a angle that is <= 180
"""
def correct_clockwise_pan(start, end):
    if angle_clockwise(start, end) <= 180:
        return start, end
    else:
        return end, start

# print correct_clockwise_pan(300, 45)
# print correct_clockwise_pan(45, 300)
# print correct_clockwise_pan(45, 180)

# if the camera directions from scene.start --> scene.end satisfy track scene
def is_tracking_scene(dirs, scene, MAX_ANGLE):
    low, high = 360, 0
    for i in range(scene.start, scene.end + 1, 1):
        low = min(low, dirs[i])
        high = max(high, dirs[i])

    if abs(high - low) < 180:
        return abs(high - low) <= MAX_ANGLE
    else:
        max_angle = 0
        for i in range(scene.start, scene.end + 1, 1):
            max_angle = max(max_angle, 360 - abs(high - dirs[i]))
        return max_angle <= MAX_ANGLE

# test
# MAX_ANGLE = 10
# dirs = [1,10,3,355]
# scene = Scene(0,3)
# print is_tracking_scene(dirs, scene)


"""
a sequence of FOVs satisfy arching scene iff exist any non-parallel pairs of camera directions that intersects at O,
for the rest of FOVs f, the angle between the two vectors (l_fO, d_f) is less than a threshold
"""
def is_arching_scene(locs, dirs, scene, MAX_ANGLE):
    # for i in range(scene.start, scene.end, 1):
    #     print locs[i], dirs[i]
    if (scene.end - scene.start) < 2:
        return False

    # print scene.end - scene.start

    comb = itertools.combinations(range(scene.start, scene.end + 1, 1), 2)
    for (i, j) in comb:
        if j != i and dirs[i] != dirs[j]:   # non-parallel pairs
            convergent = True
            # only consider pair of locations whose distance is less than 1 meter
            if distance_km(locs[i][0], locs[i][1], locs[j][0], locs[j][1]) < 0.001:
                continue
            l_o = interesection_point(locs[i], locs[j], dirs[i], dirs[j])   # find intersection point
            # print l_o
            if not l_o:
                continue
            for k in range(scene.start, scene.end + 1, 1):
                if k != i and k != j:   # the rest of FOVs
                    if distance_km(l_o[0], l_o[1], locs[k][0], locs[k][1]) > 0.01 * R_100:
                        # print "distance ", distance_km(l_o[0], l_o[1], locs[k][0], locs[k][1])
                        convergent = False
                        break
                    alpha1 = angle_bwn_points(locs[k], l_o)
                    alpha2 = correct_dir(dirs[k])
                    # print alpha1, alpha2, MAX_ANGLE
                    if not within_range(alpha1, alpha2, MAX_ANGLE):
                        convergent = False
                        break

            if convergent:  # exist one convergent center --> return True
                return True
    return False


# locs=[(34.018349, -118.291559),(34.018918, -118.291502),(34.018349, -118.290193)]
# dirs=[30,120,300]
# scene = ArchScene(0,2)
# print is_arching_scene(locs, dirs, scene, 25)


"""
A zooming scene if it satisfy both tracking scene and the camera locations are within a circle
"""
def is_zooming_scene(locs, dirs, scene, MAX_ANGLE, MAX_RADIUS):
    if (scene.end - scene.start) < 1:
        return False

    # if it is tracking scene
    if not is_tracking_scene(dirs, scene, MAX_ANGLE):
        return False

    points = [(locs[i][0], locs[i][1]) for i in range(scene.start, scene.end + 1, 1)]
    # if radius of the smallest bounding circle of the camera locations are small
    x = make_circle(points)

    if x is not None:
        cx, cy, r = x
        # print cx, cy, r
    else:
        print "xxx"
        cx, cy, r = 0, 0, 0

    # print r/ONE_KM
    if r/ONE_KM <= MAX_RADIUS:
        return True
    return False

# locs = [(34.239320, -116.951589),(34.239322, -116.951588), (34.239323, -116.951581)]
# dirs = [10,20,15]
# scene = ZoomScene (0,2)
# print is_zooming_scene(locs, dirs, scene, 15, 0.001)

"""
if the coverage angle is larger than a threshold
"""
def is_panning_angle(dirs, scene, MIN_ANGLE):
    if scene.end - scene.start < 2: # a valid panning angle needs at least three camera directions
        return False
    start, end = dirs[scene.start], dirs[scene.start + 1]
    for i in range(scene.start + 2, scene.end + 1, 1):
        if is_in_clockwise_range(dirs[i], start, end):
            continue    # do nothing
        else:
            curr_pan = correct_clockwise_pan(dirs[i - 1], dirs[i])
            start, end = combine_clockwise_ranges(start, end, curr_pan[0], curr_pan[1])

    if angle_clockwise(start, end) >= MIN_ANGLE:
        return True

    return False

# dirs = [270, 45, 0, 90, 120]
# scene = PanScene(0,4)
# print is_panning_angle(dirs, scene, 220)


"""
A panning scene if the camera locations are
 a circle and the wide angle is larger than a threshold
"""
def is_panning_scene(locs, dirs, scene, MIN_ANGLE, MAX_RADIUS):
    if (scene.end - scene.start) < 1:
        return False

    # if it is panning scene
    if not is_panning_angle(dirs, scene, MIN_ANGLE):
        return False

    points = [(locs[i][0], locs[i][1]) for i in range(scene.start, scene.end + 1, 1)]
    # if radius of the smallest bounding circle of the camera locations are small
    x = make_circle(points)

    if x is not None:
        cx, cy, r = x
        # print cx, cy, r
    else:
        print "xxx"
        cx, cy, r = 0, 0, 0

    # print r/ONE_KM
    if r/ONE_KM <= MAX_RADIUS:
        return True
    return False

def search_pan_scenes(locs, dirs, MIN_DURATION, MIN_ANGLE, MAX_RADIUS):
    if len(locs) < 2:
        return []

    results = []

    # init sequence
    sequence = PanScene(0, 0)
    i = 0
    while i < len(locs):
        # jump if possible
        if sequence.start + MIN_DURATION - 1 < len(locs) and sequence.end - sequence.start < MIN_DURATION - 1:
            i = sequence.start + MIN_DURATION - 1
        tmp_seq = PanScene(sequence.start, i)
        if is_panning_scene(locs, dirs, tmp_seq, MIN_ANGLE, MAX_RADIUS):
            sequence.end = i
        else:
            # print sequence.end - sequence.start
            if sequence.end - sequence.start >= MIN_DURATION - 1:
                results.append(sequence)

            # update start - find the largest index that satisfies
            sequence = PanScene(sequence.start + 1, i)
            j = sequence.start + 1
            for j in range(sequence.start + 1, i):
                tmp_seq = PanScene(j, i)
                if is_panning_scene(locs, dirs, tmp_seq, MIN_ANGLE, MAX_RADIUS):
                    break
            sequence.start = j

        i = i + 1
    if sequence.end - sequence.start >= MIN_DURATION - 1:
        results.append(sequence)

    return results

def search_zoom_scenes(locs, dirs, MIN_DURATION, MAX_ANGLE, MAX_RADIUS):
    if len(locs) < 2:
        return []

    results = []

    # init sequence
    sequence = ZoomScene(0, 0)
    i = 0
    while i < len(locs):
        # jump if possible
        if sequence.start + MIN_DURATION - 1 < len(locs) and sequence.end - sequence.start < MIN_DURATION - 1:
            i = sequence.start + MIN_DURATION - 1
        tmp_seq = ZoomScene(sequence.start, i)
        if is_zooming_scene(locs, dirs, tmp_seq, MAX_ANGLE, MAX_RADIUS):
            sequence.end = i
        else:
            # print sequence.end - sequence.start
            if sequence.end - sequence.start >= MIN_DURATION - 1:
                results.append(sequence)

            # update start - find the smallest index that satisfies
            sequence = ZoomScene(sequence.start + 1, i)
            for j in reversed(range(sequence.start + 1, i )):
                tmp_seq = ZoomScene(j, i)
                if not is_zooming_scene(locs, dirs, tmp_seq, MAX_ANGLE, MAX_RADIUS):
                    sequence.start= j - 1
                    break
        i = i + 1

    if sequence.end - sequence.start >= MIN_DURATION - 1:
        results.append(sequence)

    return results

# locs = [(34.239320, -116.951589),(34.239322, -116.951588), (34.239323, -116.951581)]
# dirs = [10,20,15]
# scene = ZoomScene (0,2)
# print len(search_zooming_scenes(locs, dirs, 3, 15, 0.001))


"""
Search track scenes
min-dir | max-dir | min | max | min+dir | max+dir
"""
def search_track_scenes(dirs, MIN_DURATION, MAX_ANGLE):
    results = []
    for i in range(len(dirs)):
        curr_dir = dirs[i]
        if i == 0:
            sequence = TrackScene(0, 0, curr_dir, curr_dir, 0, 0)   # (start, end, min, max, min_fov, max_fov)
        elif within(sequence.min, sequence.max, curr_dir):    # within range
            sequence.end = i         # update end
        elif within(angle(sequence.max - MAX_ANGLE), sequence.min, curr_dir):  # left
            sequence.end, sequence.min, sequence.min_fov = i, curr_dir, i  # update min, min_fov
        elif within(sequence.max, angle(sequence.min + MAX_ANGLE), curr_dir):   # right
            sequence.end, sequence.max, sequence.max_fov = i, curr_dir, i  # update max, max_fov
        elif within(angle(sequence.min - MAX_ANGLE), angle(sequence.max - MAX_ANGLE), curr_dir):
            if sequence.end - sequence.start >= MIN_DURATION - 1:  # direction scene is detected
                results.append(sequence)
            # update max - find the most recent dir that violates the range
            max_fov = i
            max = curr_dir
            for j in reversed(range(sequence.start, i)):
                # if less_than(max, dirs[j]) and within_range(curr_dir, dirs[j], MAX_ANGLE - 1):
                #     max, max_fov = dirs[j], j
                if within_range(curr_dir, dirs[j], MAX_ANGLE - 1):
                    if less_than(max, dirs[j]):
                        max, max_fov = dirs[j], j
                    continue
                else:
                    break
            sequence = TrackScene(j + 1, i, curr_dir, max, i, max_fov)
        elif within(angle(sequence.min + MAX_ANGLE), angle(sequence.max + MAX_ANGLE), curr_dir):
            if sequence.end - sequence.start >= MIN_DURATION - 1:  # direction scene is detected
                results.append(sequence)
            # update min - find the most recent dir that violates the range
            min_fov = i
            min = curr_dir
            for j in reversed(range(sequence.start, i)):
                # if less_than(dirs[j], min) and within_range(dirs[j], curr_dir, MAX_ANGLE):
                #     min, min_fov = dirs[j], j
                if within_range(dirs[j], curr_dir, MAX_ANGLE):
                    if less_than(dirs[j], min):
                        min, min_fov = dirs[j], j
                    continue
                else:
                    break
            sequence = TrackScene(j + 1, i, min, curr_dir, min_fov, i)

        elif not within_range(curr_dir, sequence.min, MAX_ANGLE) or not within_range(sequence.max, curr_dir, MAX_ANGLE):    # out of range
            if sequence.end - sequence.start >= MIN_DURATION - 1:  # direction scene is detected
                results.append(sequence)
            sequence = TrackScene(i, i, curr_dir, curr_dir, i, i) # create new scene

    return results




"""
Search arching scenes
min-dir | max-dir | min | max | min+dir | max+dir
"""
def search_arch_scenes(locs, dirs, MIN_DURATION, MAX_ANGLE):
    if len(locs) < 3:
        return []

    results = []

    # init sequence with the first two
    sequence = ArchScene(0, 1)
    i = 0
    while i < len(locs):
        # jump if possible
        if sequence.start + MIN_DURATION - 1 < len(locs) and sequence.end - sequence.start < MIN_DURATION - 1:
            i = sequence.start + MIN_DURATION - 1
        tmp_seq = ArchScene(sequence.start, i)
        if is_arching_scene(locs, dirs, tmp_seq, MAX_ANGLE):
            sequence.end = i
        else:
            # print sequence.end - sequence.start
            if sequence.end - sequence.start >= MIN_DURATION - 1:
                results.append(sequence)

            # update start - find the smallest index that satisfy arching scene
            sequence = ArchScene(sequence.start + 1, i)
            for j in reversed(range(sequence.start + 1, i - 1)):
                tmp_seq = ArchScene(j, i)
                if not is_arching_scene(locs, dirs, tmp_seq, MAX_ANGLE):
                    sequence.start= j - 1
                    break
        i = i + 1

    if sequence.end - sequence.start >= MIN_DURATION - 1:
        results.append(sequence)

    return results


# locs=[(34.018918, -118.291502), (34.018349, -118.290193),(34.018349, -118.291559),(34.018918, -118.291502),(34.018349, -118.290193), (34.018349, -118.291559)]
# dirs=[90,0,30,120,300,60]
# print search_arch_scenes(locs, dirs, 3, 25)[0]


# MAX_ANGLE = 15      # degree
# MIN_DURATION = 4
# dirs = [120, 150, 160, 140, 150, 145, 135, 130, 150, 120]
# dirs = [20, 5, 0, 350, 355, 10, 0, 5, 20, 15]
# dirs = [181.94008, 182.37218, 171.61913, 160.90152, 157.55756, 153.66913, 156.21973, 175.0801, 159.62263, 153.34453, 151.51772, 151.01799, 152.5305, 155.17336, 154.8254, 150.19846, 166.72073, 158.32208, 148.25983, 154.22003, 162.4578, 167.7894, 181.1051, 170.45929, 175.29805, 175.39275, 186.90277, 188.76363, 200.17162, 200.49422, 194.32492, 197.19994, 202.31071, 204.23193, 197.63225, 198.6546, 205.44516, 198.59033, 200.10536, 200.5626, 202.47163, 204.49399, 204.32697, 205.05452, 209.90833, 210.01016, 220.22246, 240.85333, 247.69794, 254.44272, 265.39795, 270.58826, 274.80835, 270.61554, 266.99905, 194.38068, 161.73622, 168.23277, 176.72787, 187.61053, 204.8588]
#  dirs[31:43]
# results = search_tracking_scenes(dirs)
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


"""
videoid 0
fovid 1
lat 2
lon 3
dir 4
r 7
theta 8
timestamp 9
"""

def read_fovs(file):
    data = np.genfromtxt(file, dtype=None, unpack=True, delimiter='\t')

    prev_vid = data[0][0]
    fovs = []
    idx = 0
    videos = []
    for d in data:
        vid = d[0]
        # print str(vid), str(prev_vid)
        if vid == prev_vid:
            # lat, lon, compass, R, alpha
            fov = FOV(data[idx][2],data[idx][3],data[idx][4],data[idx][7],data[idx][8])
            fovs.append(fov)
        else:
            # new video
            v = Video(fovs)
            v.id = prev_vid
            videos.append(v)
            # print v.to_str()

            # new fovs
            fovs = []
            fov = FOV(data[idx][2],data[idx][3],data[idx][4],data[idx][7],data[idx][8])
            fovs.append(fov)

        idx = idx + 1
        prev_vid = vid

    print "number of videos", len(videos)
    return videos



# Create geoq client api
geoq = sm.GeoqApi(sm.ApiClient(url))

# Returns a set of video locations in GEOJSON format (small sample data)
#geoq.sample_videos()

# Returns a set of video frames (of a particular video) in GEOJSON format (small sample data)
#geoq.sample_fovs()

# Create geoq client with API key
# Replace KEY_VALUE by actual one

geoq = sm.GeoqApi(sm.ApiClient(url, "X-API-KEY", "8b51UFM2SlBltx3s6864eUO1zSoefeK5"))



"""
Returns a set of video locations that are captured within a time interval (startdate -> enddate)
"""
# def get_videos(swlat=34.018212, swlng=-118.291716, nelat=34.025296, nelng=-118.279826, startdate="2014-04-13 00:00:00", enddate="2014-04-13 23:59:59"):

def get_videos(swlat=-90, swlng=-180, nelat=90, nelng=+180, startdate="2010-04-13 00:00:00", enddate="2017-04-13 23:59:59"):

    valid_fc_videos = [] # with video of size > 0

    if IS_FROM_MEDIAQ:
        fc_videos = geoq.rectangle_query(swlat, swlng, nelat, nelng)
        fc_videos = fc_videos.replace('None','null').replace('u\'','\"').replace('\'','\"')

        # print fc_videos

        # validate GEOJSON
        geojson_validation = requests.post(validate_endpoint, data=fc_videos)
        if geojson_validation.json()['status'] != 'ok':
            print "Rectangle_query: Invalid geojson format"
            exit()

        fc_videos = geojson.loads(fc_videos)
        print "Number of videos: " + str(len(fc_videos.features))

        # fov_counts = []
        video_sizes = []
        for video in fc_videos.features:
            # fov_counts.append(video.properties['fov_count'])
            size = size_in_mb(video.properties['size'])
            if size and size > 0:
                print video
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

        # print "Total/Mean FOV count: " + str(sum(fov_counts)) + "/" + str(np.mean(fov_counts))
        # bins = range(0,200,20)
        # hist_fov_counts = np.histogram(fov_counts, bins)[0]
        # print "Histogram:" + str(hist_fov_counts)

        # the histogram of the data
        # plt.hist(fov_counts, bins=bins)
        # plt.xlabel('Range (Number of Frames)')
        # plt.ylabel('Frame Count')
        # plt.title(r'$\mathrm{Histogram\ of\ Frame\ Count}$')
        # plt.axis([0, 100, 0, 500])
        plt.grid(True)
        # plt.show()

    else:
        valid_fc_videos = read_fovs(VIDEO_FILE)

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

"""
Search for various kinds of video scenes,
"""
def search_scenes(video, type, MIN_DURATION, ANGLE_THRESHOLD):

    if IS_FROM_MEDIAQ:
        if video.properties['vid']:
            vid = str(video.properties['vid'])
            videoid = str(video.properties['videoid'])
            fovs = getFOVs(vid)

            if fovs and len(fovs) > 0:
                # print vid
                trajectory = [(fov.geometry.coordinates[0],fov.geometry.coordinates[1],float(fov.properties['theta_x'])) for fov in fovs.features]
                dirs = [dir[2] for dir in trajectory]
                locs = [(dir[0], dir[1]) for dir in trajectory]



    else:
        dirs = [fov.compass for fov in video.fovs]
        locs = [(fov.lat, fov.lon) for fov in video.fovs]

    if type == "track":
        results = search_track_scenes(dirs, MIN_DURATION, ANGLE_THRESHOLD)
    elif type == "arch":
        results = search_arch_scenes(locs, dirs, MIN_DURATION, ANGLE_THRESHOLD)
    elif type == "zoom":
        results = search_zoom_scenes(locs, dirs, MIN_DURATION, ANGLE_THRESHOLD, MAX_RADIUS)
    elif type == "pan":
        results = search_pan_scenes(locs, dirs, MIN_DURATION, ANGLE_THRESHOLD, MAX_RADIUS)

    return len(results)
    
    # for s in results:
    #     print geoq.video_segment_url(vid, s.start, s.end)

    # results = search_zooming_scenes(locs, dirs)
    # for s in results:
    #     print geoq.video_segment_url(vid, s.start, s.end)



"""
Dump mediaq's FOV data in CSV format
"""
def dump_mediaq_metadata(filename = "mediaq_fovs.txt"):

    file = open(filename,"w")

    videos = get_videos()
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

"""
Compute coverage map (2D histogram)
"""
def compute_coverage_map(grid_size = 20):
    swlat=34.018212
    swlng=-118.291716
    nelat=34.025296
    nelng=-118.279826
    videos = get_videos(swlat, swlng, nelat, nelng)

    map = np.ndarray(shape=(grid_size, grid_size), dtype=int)
    for video in videos:
        if IS_FROM_MEDIAQ:
            if video.properties['vid']:
                vid = str(video.properties['vid'])
                fovs = getFOVs(vid)

                if fovs and len(fovs) > 0:
                    for fov in fovs.features:
                        f = FOV(fov)
                        param = Params(200, swlat, swlng, nelat, nelng)
                        param.GRID_SIZE = grid_size
                        for cid in f.cellids(param):
                            cell_lat, cell_lng = cell_coord(cid, param)
                            if f.cover(cell_lat, cell_lng):
                                y_idx = cid/param.GRID_SIZE
                                x_idx = cid - y_idx*param.GRID_SIZE
                                # print x_idx, y_idx, map[x_idx][y_idx]
                                map[x_idx][y_idx] = map[x_idx][y_idx] + 1
        else:
            for f in video.fovs:
                param = Params(200, swlat, swlng, nelat, nelng)
                param.GRID_SIZE = grid_size
                for cid in f.cellids(param):
                    cell_lat, cell_lng = cell_coord(cid, param)
                    if f.cover(cell_lat, cell_lng):
                        y_idx = cid/param.GRID_SIZE
                        x_idx = cid - y_idx*param.GRID_SIZE
                        # print x_idx, y_idx, map[x_idx][y_idx]
                        map[x_idx][y_idx] = map[x_idx][y_idx] + 1

    fig, ax = plt.subplots()
    heatmap = ax.pcolor(map, cmap=plt.cm.Reds)
    plt.show()
    np.savetxt("mediaq_coverage_heatmap.txt" , map, fmt='%i\t')


# dump_mediaq_metadata()

start = time.time()
videos = get_videos()
print time.time() - start


# min_duration = [5,10,15,20,25,30,35,40,45,50,55,60]
min_duration = [15]
# angle_threshold = [5,10,15,20,25,30,35,40,45]
angle_threshold = [15]
# angle_threshold = [60, 90,120,150,180,210,240,270,300]

# tracking scenes
TEST = False

start = time.time()
type = "arch"
if TEST:
    scenes = np.zeros((len(min_duration), len(angle_threshold)))

    for i in range(len(min_duration)):
        for j in range(len(angle_threshold)):
            found = 0
            num = 0
            for video in videos:
                num = num + 1

                # check validity of the data
                if any(v.compass < 0 or v.compass > 360 for v in video.fovs): # do nothing if exist camera direction is negative
                    # print video.id, [v.compass for v in video.fovs]
                    continue

                if len(video.fovs) > 2000:  # do not consider too long videos
                    # print "do not consider too long videos ", len(video.fovs)
                    continue

                # print num, len(video.fovs)
                found = found + search_scenes(video, type, min_duration[i], angle_threshold[j])

                if num % 100 == 0:
                    print "done ", num, " videos"

            scenes[i,j] = found

    np.savetxt(type + "_scenes.txt" , scenes, fmt='%.4f\t')
print time.time() - start

compute_coverage_map()

