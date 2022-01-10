import numpy as np
import cv2
import pickle
import os

def get_radius(im):
    # Returns the radius of a bubble in pixels
    contours, _ = cv2.findContours(im.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Get the coordinates of the center of the bubble contour
    centers = []
    for cont in contours:
        M = cv2.moments(cont)
        try:
            center = (M['m10']/M['m00'], M['m01']/M['m00'])
        except ZeroDivisionError:
            continue
        centers.append(center)

    if len(centers) > 1:
        center = np.array([0.,0.], dtype=float)
        for i in range(len(centers)):
            center += np.array(centers[i])
        center /= len(centers)
    elif len(centers) == 1:
        center = np.array(centers[0], dtype=float)
    else:
        center = None
        
    # Radius is defined as the median distance of the center of the bubble to the edge of the bubble contour
    radius = 0
    if center is not None:
        dists = []
        for cont in contours:
            for point in cont:
                dists.append(np.linalg.norm(center - point))
        radius = np.median(dists)

    return radius


im_dir = './singles_new_scale_mar21/'

with open(im_dir + 'data_dict.pkl', 'rb') as f:
    data_dict = pickle.load(f)

rad_dict = {}
im_names = [name for name in os.listdir(im_dir) if '.png' in name]
for i, im_name in enumerate(im_names):
    im = cv2.imread(im_dir + im_name, cv2.IMREAD_GRAYSCALE)
    radius = get_radius(im)
    
    name = im_name.split('-')
    data_key = (int(name[0]), int(name[2][0])-1)
    
    coord_pos = list(data_dict[data_key])
    
    coord_pos.append(radius)
    coord_pos_rad = tuple(coord_pos)
    
    rad_dict.update( {(int(name[0]), int(name[2][0])-1, int(name[1])) : coord_pos_rad} )
    
    if i % 1000 == 0:
        print('Done %i of %i images' % (i, len(im_names)))

with open(im_dir + 'data_with_rad_dict.pkl', 'wb') as f:
    pickle.dump(rad_dict, f)
