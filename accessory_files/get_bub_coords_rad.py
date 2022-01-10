import cv2
import os
import numpy as np
import imutils

imdir = './data/full_single_33/'
write_file = 'single_coord_rad_33.txt'

im_names = os.listdir(imdir)

coord_dict = {}
run_ev = []
for name in im_names:
    im = cv2.imread(imdir + name, cv2.IMREAD_GRAYSCALE)
    contours, hierarchy = cv2.findContours(im.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
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
        
    name = name.split('-')
    
    coord_dict.update( {(name[0], int(name[1]), int(name[2][0])) : (center, radius)} )
    if (name[0], int(name[1])) not in run_ev:
        run_ev.append( (name[0], int(name[1])) )

coord_rad_dict = {}
for ev_key in run_ev:
    cam_info = [[], [], [], []]
    radii = []
    for key in coord_dict:
        if key[:2] != ev_key:
            continue
        
        cam_info[key[2]].append(coord_dict[key][0])
        radii.append(coord_dict[key][1])
    radius = np.mean(radii)
    '''
    RADIUS IS INCORRECT.  CANNOT AVERAGE OVER ALL CAMERAS SINCE BUBBLE RADIUS IS MEASURED IN PIXELS.
    '''
    for i in range(4):
        if cam_info[i] == []:
            cam_info[i].append(np.array([-100., -100.]))
    coord_rad_dict.update( {ev_key : [cam_info[0][0], cam_info[1][0], cam_info[2][0], cam_info[3][0], radius]} )


with open('./data/' + write_file, 'w') as f:
 for key in coord_rad_dict:
    line = key[0] + ' ' + str(key[1]) + ' '
    for item in coord_rad_dict[key]:
        if isinstance(item, np.ndarray):
            line += str(item[0]) + ' ' + str(item[1]) + ' '
        else:
            line += str(item) + '\n'
    f.write(line)
    
