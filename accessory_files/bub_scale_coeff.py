import numpy as np
import os
import cv2
import matplotlib.pyplot as plt

with open('data/single_pos_24.txt', 'r') as f:
    data = f.readlines()

for i in range(len(data)):
    data[i] = data[i].split()

position_dict = {}
for bubble in data:
    position_dict.update({(bubble[0], int(bubble[1])) : [float(bubble[2])/10, float(bubble[3])/10, float(bubble[4])/10]})

cam_positions = np.array(
                [[-19.0500-0.8351, -32.9956-1.4465, 40.8681],
                [-19.0500-0.7743, -32.9956-1.3411, 13.0119],
                [19.0500+1.2247, -32.9956-2.1212, 40.9293],
                [19.0500+1.6113, -32.9956-2.7908, 12.8003]]
                )
dist_to_camera = lambda p, cam: np.linalg.norm(cam_positions[cam] - p)

#def get_area_sum(bubble):
#    contours, _ = cv2.findContours(bubble, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
#    A = 0
#    for cont in contours:
#        A += cv2.contourArea(cont)
#    return A
#
#def get_area_max(bubble):
#    contours, _ = cv2.findContours(bubble, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
#    biggest_cont = max(contours, key=cv2.contourArea)
#    return cv2.contourArea(biggest_cont)
#
#
#data_dir = 'data/bub_single_24new/'
#
#area_dict = {}
#for bubble in os.listdir(data_dir):
#    img = cv2.imread(data_dir + bubble, cv2.IMREAD_GRAYSCALE)
#    area_dict.update({(bubble[:10], int(bubble[11:-8]), int(bubble[-7])) : (get_area_sum(img), get_area_max(img))})
#
#
#
#area_sum, area_max, distances = [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]
#for key in area_dict.keys():
#    cam = key[-1]
#    area_sum[cam].append(area_dict[key][0])
#    area_max[cam].append(area_dict[key][1])
#    distances[cam].append(dist_to_camera(position_dict[key[:2]], cam))
#
#all_areas_sum = area_sum[0] + area_sum[1] + area_sum[2] + area_sum[3]
#all_areas_max = area_max[0] + area_max[1] + area_max[2] + area_max[3]
#all_distances = distances[0] + distances[1] + distances[2] + distances[3]
#
#plt.figure(figsize=(12,9))
#plt.scatter(all_distances, all_areas_max, s=2, c='r')
#plt.scatter(all_distances, all_areas_sum, s=2, c='b')
#
#fit_coeff = np.polyfit(all_distances, all_areas_sum, 1)
#print(fit_coeff)
#
#x = np.linspace(25, 65, 200)
#y = fit_coeff[0]*x + fit_coeff[1]
#
#plt.plot(x, y, c='g')
#
#plt.show()



def get_radius(im):
    contours, _ = cv2.findContours(im.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # Get the coordinates of the center of the bubble
    centers = []
    for cont in contours:
        M = cv2.moments(cont)
        try:
            center = (M['m10']/M['m00'], M['m01']/M['m00'])
        except ZeroDivisionError:
            continue
        centers.append(center)
    
    if len(centers) > 1:
        center = np.array([0,0], dtype=float)
        for i in range(len(centers)):
            center += np.array(centers[i])
        center /= len(centers)
    
    if len(centers) == 1:
        center = np.array(centers[0], dtype=float)

    if len(centers) == 0:
        center = np.array([-100., -100.])
        
        
    radius = 0
    if sum(center) > 0:
        dists = []
        for cont in contours:
            for point in cont:
                dists.append(np.linalg.norm(center - point))
        radius = np.median(dists)
    
    return radius

data_dir = 'data/bub_single_24new/'

radius_dict = {}
for bubble in os.listdir(data_dir):
    img = cv2.imread(data_dir + bubble, cv2.IMREAD_GRAYSCALE)
    radius_dict.update({(bubble[:10], int(bubble[11:-8]), int(bubble[-7])) : get_radius(img)})

rad = []
dist = []
for key in radius_dict:
    rad.append(radius_dict[key])
    dist.append(dist_to_camera(position_dict[key[:2]], key[2]))

plt.figure(figsize=(12,9))
plt.scatter(dist, rad, s=2)
plt.xlabel('Distance to camera [cm]')
plt.ylabel('Bubble radius [pixels]')

fit_coeff = np.polyfit(dist, rad, 1)
print(fit_coeff)

x = np.linspace(24, 69, 200)
y = fit_coeff[0]*x + fit_coeff[1]

plt.plot(x, y, c='g')

plt.show()

