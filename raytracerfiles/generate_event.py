import matlab.engine
import pickle
import os
import numpy as np
import cv2
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--n_events", type=int, default=1, help="Number of events to generate. Default: 1")
parser.add_argument("--mult", type=int, default=1, help="Multiplicity of generated events. Default: 1")
parser.add_argument("--separate_events", type=int, default=0, help="Set to 1 to separate events in different folders, each with their own dictionary and text file.  Set to 0 to generate all events in one directory. Default: 0")
parser.add_argument("--omit_txt", type=float, default=0, help="Set to 1 to omit the text file containing event data. Default: 0")
parser.add_argument("--omit_dict", type=float, default=0, help="Set to 1 to omit the pickle file containing a dictionary of event data. Default: 0")
parser.add_argument("--omit_imgs", type=float, default=0, help="Set to 1 to omit the generation of images. Default: 0")
parser.add_argument("--image_dir", type=str, default="../dcgan_images/", help="Directory in which the 32x32 pixel images of bubbles are saved. Default: ../dcgan_images/")
parser.add_argument("--write_dir", type=str, default="../artificial_events/", help="Directory in which the artificial events are saved. Default: ../artificial_events/")
opt = parser.parse_args()

events = opt.n_events

# int or list of ints of length events
#multiplicity = [np.random.randint(1, 11) for _ in range(events)]
multiplicity = opt.mult

# Separate events in different folders, each with their own dictionary and text file.
separate_events = opt.separate_events

omit_txt = opt.omit_txt
omit_dict = opt.omit_dict
omit_imgs = opt.omit_imgs

# Directory of the 32x32 artificial images of bubbles
image_dir = opt.image_dir

# Directory to write atrificial image files to
write_dir = opt.write_dir
os.makedirs(write_dir, exist_ok=True)

# If fit parameters are saved, load them. Otherwise, compute them.
new_params = False
try:
    with open('fitparams.pkl', 'rb') as f:
        allparams, cam1params, cam2params, cam3params, cam4params = pickle.load(f)
    print('Read file fitparams.pkl.')
except:
    print('Computing fitparams.')
    eng = matlab.engine.start_matlab()
    new_params = True
    allparams, cam1params, cam2params, cam3params, cam4params = eng.FitPICO60Geometry(nargout=5)
    with open('fitparams.pkl', 'wb') as f:
        pickle.dump([allparams, cam1params, cam2params, cam3params, cam4params], f)
        
try:
    assert not new_params
    with open('rays.pkl', 'rb') as f:
        _, _, pixels, gen_params = pickle.load(f)
    del _
    print('Read file rays.pkl.')
except:
    print('Computing rays.')
    if not 'eng' in locals() and not 'eng' in globals():
        eng = matlab.engine.start_matlab()
    _, _, pixels, gen_params = eng.GetRaysAndPixels(allparams, [1700,1088], nargout=4)
    del _
    with open('rays.pkl', 'wb') as f:
        pickle.dump([scatters, ray_startingpoints, pixels, gen_params], f)

if 'eng' in locals() or 'eng' in globals():
    eng.quit()


r1, r2, r3, liquid_level, cyl_axis = gen_params[0][1], gen_params[0][3], gen_params[0][5], gen_params[0][6], gen_params[0][7:]
# r1 is main inner cylinder radius
# r2 is minor radius of torus defining the inner surface of the knuckle of the jar
# r3 is the radius of the sphere defining the inner surface of the dome of the jar
# liquid_level is the level of the water-freon interface
s = r3 * (r1 - r2)/(r3 - r2)
z = -r2 * np.sqrt(1 - (s/r3)**2)
d = r3 * z * (1/r3 - 1/r2)
cyl_axis = np.array(cyl_axis)
resolution = (1700, 1088)  # Camera resolution

# Camera positions
# Taken from GetRaysAndPixels.m
cam_positions = np.array(
                [[-19.0500-0.8351, -32.9956-1.4465, 40.8681],
                [-19.0500-0.7743, -32.9956-1.3411, 13.0119],
                [19.0500+1.2247, -32.9956-2.1212, 40.9293],
                [19.0500+1.6113, -32.9956-2.7908, 12.8003]]
                )

# Size of the bubble images
# Should be a positive even integer
bs = 32

# Size of the partition boxes in cm
ps = 3

# These are the dimensions of a box bounding the volume of the detector in which we wish to generate points
box = (30., 30., 96.)  # Each dimension should be an integer multiple of the partition size 'ps'
# Coordinate of the center of the box
box_center = (0., 0., 38.)
# These should be the same as specified in partition_rays.py
# These need to be tuples of floats, not ints

def injar(p, r1, r2, r3, cyl_axis, z, d):
    # Returns true if point p is inside the jar, false otherwise
    if p[2] > liquid_level:
        return False
    if p[0]**2 + p[1]**2 < r1**2 and p[2] >= 0:
        return True
    if np.linalg.norm(p - cyl_axis*d) < r3 and np.dot(cyl_axis, p) < z:
        return True
    if np.dot(cyl_axis, p) < 0 and np.dot(cyl_axis, p) >= z:
        if p[0]**2 + p[1]**2 < r3**2 - (np.abs(z) + cyl_axis[2]*d)**2:
            return True
        if (np.dot(p,p) - (r1 - r2)**2 - r2**2)**2 - 4*(r1 - r2)**2*(r2**2 - p[2]**2) < 0:
            return True
    return False


def generate_points(n, r1, r2, r3, cyl_axis, z, d):
    # Returns an n-by-3 array of random points within the detector geometry
    
    if n == 1:
        p = [box[0]*np.random.rand() -  box[0]/2 + box_center[0],
             box[1]*np.random.rand() -  box[1]/2 + box_center[1],
             box[2]*np.random.rand() -  box[2]/2 + box_center[2]]
        while not injar(p, r1, r2, r3, cyl_axis, z, d):
            p = [box[0]*np.random.rand() -  box[0]/2 + box_center[0],
                 box[1]*np.random.rand() -  box[1]/2 + box_center[1],
                 box[2]*np.random.rand() -  box[2]/2 + box_center[2]]
        return np.array(p)
    
    pts = np.zeros([n,3])
    for i in range(n):
        p = [box[0]*np.random.rand() -  box[0]/2 + box_center[0],
             box[1]*np.random.rand() -  box[1]/2 + box_center[1],
             box[2]*np.random.rand() -  box[2]/2 + box_center[2]]
        while not injar(p, r1, r2, r3, cyl_axis, z, d):
            p = [box[0]*np.random.rand() -  box[0]/2 + box_center[0],
                 box[1]*np.random.rand() -  box[1]/2 + box_center[1],
                 box[2]*np.random.rand() -  box[2]/2 + box_center[2]]
        pts[i] = p
    return pts
    
    
partition_coords = [[i,j,k] for i in range(int(box[0]//ps)) for j in range(int(box[1]//ps)) for k in range(int(box[2]//ps))]
partition_ray_dicts = [{} for _ in range(len(partition_coords))]

# List of pairs of coordinates, minimal and maximal, defining each box in the partition
box_coords = []
for coords in partition_coords:
    B1 = (int((coords[0]*ps - box[0]/2 + box_center[0])*100),
          int((coords[1]*ps - box[1]/2 + box_center[1])*100),
          int((coords[2]*ps - box[2]/2 + box_center[2])*100))
    box_coords.append(B1)

box_file_dict = {}
for cam in range(4):
    for i, B1 in enumerate(box_coords):
        box_file_dict.update( {(cam, B1) : 'cam%ibox%i.pkl' % (cam, i)} )

def get_ray_box(p, cam):
    # Returns the ray file associated with the partition containing p
    p = ((p[0] + box[0]/2 - box_center[0])//ps,
         (p[1] + box[1]/2 - box_center[1])//ps,
         (p[2] + box[2]/2 - box_center[2])//ps)
    B1 = (int((p[0]*ps - box[0]/2 + box_center[0])*100),
          int((p[1]*ps - box[1]/2 + box_center[1])*100),
          int((p[2]*ps - box[2]/2 + box_center[2])*100))
    file_name = box_file_dict[(cam, B1)]
    with open('partitions/' + file_name, 'rb') as f:
        return pickle.load(f)


def ray_to_point_dist(ray, point):
    # Return the squared distance from the point to the ray where the ray is defined by two points on either end.
    numer = np.cross(ray[0] - point, ray[1] - point)
    denom = ray[1] - ray[0]
    return np.dot(numer, numer)/np.dot(denom, denom)

def closest_ray_indices(p):
    ray_indices = []
    for cam in range(4):
        rays = get_ray_box(p, cam)
        if not rays:
            ray_indices.append(-1)
            continue
        min_ray_index = min(rays.keys(), key=lambda ray_key: ray_to_point_dist(rays[ray_key], p))
        if ray_to_point_dist(np.array(rays[min_ray_index]), p) > 0.0016:  # Determined empirically
            ray_indices.append(-1)
            continue
        ray_indices.append(min_ray_index)
    return ray_indices

print('Generating positions and computing rays.')
# Compute the ray closest to a random point for each camera.
t1 = time.time()
if isinstance(multiplicity, int):
    n = events * multiplicity
elif isinstance(multiplicity, list) and len(multiplicity) == events:
    n = sum(multiplicity)
else:
    print('\nImproper multiplicity specified.\n')
    raise ValueError
pts, ray_indices = [], []
while len(ray_indices) < n:
    p = generate_points(1, r1, r2, r3, cyl_axis, z, d)
    rays = closest_ray_indices(p)
    if sum(rays) > 0:
        pts.append(p)
        ray_indices.append(rays)
    if len(ray_indices) % 100 == 0:
        print('   Computed %i of %i bubble positions' % (len(ray_indices), n))
pts = np.array(pts)

print('Computing pixel coordinates.')
# Get the pixel coordinates of the ray passing through the point
pixel_coords = []
for indices in ray_indices:
    event_pixels = []
    for index in indices:
        if index == -1:
            pix = [-1,-1]
        else:
            pix = np.array(pixels[index][1:], dtype=np.int)
        event_pixels.append(pix)
    pixel_coords.append(event_pixels)
pixel_coords = np.array(pixel_coords)
print('%i pixel coordinates computed in %f seconds' % (4*n, time.time()-t1))

del ray_indices

'''
pts is an nx3 array of random 3D points chosen within the jar.
pixel_coords is an nx4x2 array of pixel coordinates corresponding to the n points.  The second dimension (4) refers to the number of cameras.  If the point is not visible in a camera, the pixel coordinate will be reported as [-1,-1].
'''

def dist_to_camera(p, cam):
    return np.linalg.norm(cam_positions[cam] - p)

## Old method of scaling bubble.  Uses bubble area in pixels.  New method uses radius.
#def scale_size(dist):
#    # Return the area of a bubble from the fit of the data
#    # dist is in cm and area is in pixels
#    return -0.1467404*dist + 23.95305973
#
#def scale_bubble(bubble, p, cam):
#    # Scale a bubble according to its position p and distance from camera cam
#    size = np.sum(bubble/255)
#    dist = dist_to_camera(p, cam)
#    scaled_size = scale_size(dist)
#    scale_factor = np.sqrt(scaled_size/size)
#    if scale_factor < np.inf:
#        new_dim = int(bs*scale_factor)
#        resized_bub = cv2.resize(bubble.copy(), (new_dim, new_dim))
#        _, resized_bub = cv2.threshold(resized_bub, 75, 255, cv2.THRESH_BINARY)
#    else:
#        new_dim = bs
#        resized_bub = np.zeros((bs,bs))
#
#    center = new_dim//2
#    if center < bs//2:
#        bub_img = np.zeros((bs,bs))
#        bub_img[bs//2-center:bs//2-center+new_dim, bs//2-center:bs//2-center+new_dim] += resized_bub
#        return bub_img
#    return resized_bub[center-bs//2:center+bs//2, center-bs//2:center+bs//2]


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

# Emperically determined coefficients for each camera for the (assumed) linear relationship between bubble radius and distance from camera
m = [-0.02120990766348932, -0.023366570364067447, -0.013947505326630556, -0.017606693132123528]
b = [3.196533075363447, 3.0463249066649655, 3.579675796810233, 3.069843681175793]
def scale_rad(dist):
    # Return the radius of a bubble from the fit of the data
    # dist is in cm and radius is in pixels
    return m[cam]*dist + b[cam]

def scale_bubble(bubble, p, cam):
    # Scale a bubble according to its position p and distance from camera cam
    rad = get_radius(bubble)
    dist = dist_to_camera(p, cam)
    scaled_rad = scale_rad(dist)
    scale_factor = scaled_rad/rad
    if scale_factor < np.inf:
        new_dim = int(bs*scale_factor*0.9)
        resized_bub = cv2.resize(bubble.copy(), (new_dim, new_dim))
        _, resized_bub = cv2.threshold(resized_bub, 50, 255, cv2.THRESH_BINARY)
    else:
        new_dim = bs
        resized_bub = np.zeros((bs,bs))
        
    center = new_dim//2
    if center < bs//2:
        bub_img = np.zeros((bs, bs))
        bub_img[bs//2-center:bs//2-center+new_dim, bs//2-center:bs//2-center+new_dim] += resized_bub
        return bub_img
    return resized_bub[center-bs//2:center+bs//2, center-bs//2:center+bs//2]
    

print('Generating events.')
all_bubs = os.listdir(image_dir)
full_data_dict = {}
for event in range(events):
    event_ims = np.zeros((4, *resolution))
    event_data_dict = {}
    if isinstance(multiplicity, int):
        for bub_num in range(multiplicity):
            index = event*multiplicity + bub_num
            position = pts[index]
            pix = pixel_coords[index]
            
            # Update a dictionary containing (event number, bubble number) as keys and (3D bubble position, pixel coordinates in all cameras) as values
            event_data_dict.update( {(event, bub_num) : (position, pix)} )
            full_data_dict.update( {(event, bub_num) : (position, pix)} )
            
            for cam in range(4):
                cam_pix = pix[cam].copy()
                if np.sum(cam_pix) < 0:
                    continue
                    
                # This accounts for the fact that the image is rotated prior to being saved (to match the data)
                cam_pix[1] = resolution[1] - cam_pix[1]
                
                bubble = np.random.choice(all_bubs)
                bubble = cv2.imread(image_dir + bubble, cv2.IMREAD_GRAYSCALE)
                
                bub_scaled = scale_bubble(bubble, position, cam)
                if cam_pix[0] < bs//2:
                    bub_scaled = bub_scaled[bs//2-cam_pix[0]:, :]
                    event_ims[cam, :cam_pix[0]+bs//2, cam_pix[1]-bs//2:cam_pix[1]+bs//2] += bub_scaled
                    continue
                if cam_pix[1] < bs//2:
                    bub_scaled = bub_scaled[:, bs//2-cam_pix[0]:]
                    event_ims[cam, cam_pix[0]-bs//2:cam_pix[0]+bs//2, :cam_pix[1]+bs//2] += bub_scaled
                    continue
                if cam_pix[0] + bs//2 > resolution[0]:
                    bub_scaled = bub_scaled[:resolution[0]+bs//2-cam_pix[0], :]
                    event_ims[cam, cam_pix[0]-bs//2:, cam_pix[1]-bs//2:cam_pix[1]+bs//2] += bub_scaled
                    continue
                if cam_pix[1] + bs//2 > resolution[1]:
                    bub_scaled = bub_scaled[:, :resolution[1]+bs//2-cam_pix[0]]
                    event_ims[cam, cam_pix[0]-bs//2:cam_pix[0]+bs//2, cam_pix[1]-bs//2:] += bub_scaled
                    continue
                    
                event_ims[cam, cam_pix[0]-bs//2:cam_pix[0]+bs//2, cam_pix[1]-bs//2:cam_pix[1]+bs//2] += bub_scaled
                
        if separate_events:
            # Create an event directory to save images and event data
            if not omit_imgs or not omit_dict or not omit_txt:
                os.makedirs(write_dir + '%i-%i' % (event, multiplicity), exist_ok=True)
            
            # Write the images to png files.  Image names have the format event-camera-multiplicity.png
            if not omit_imgs:
                for cam in range(4):
                    cv2.imwrite(write_dir + '%i-%i/%i-%i-%i.png' % (event, multiplicity, event, cam, multiplicity), np.rot90(event_ims[cam]))
            
            # Write the event data dictionary to a pickle file for easy access in future python programs
            if not omit_dict:
                with open(write_dir + '%i-%i/data_dict.pkl' % (event, multiplicity), 'wb') as f:
                    pickle.dump(event_data_dict, f)
            
            # Write the bubble 3D positions and pixel coordinates in the event to a text file saved in the event directory
            if not omit_txt:
                with open(write_dir + '%i-%i/data_text.txt' % (event, multiplicity), 'w') as f:
                    for key in event_data_dict:
                        line = r''
                        for i in range(2):
                            for num in event_data_dict[key][i]:
                                if i == 0:
                                    line += str(num) + ' '
                                else:
                                    for j in range(2):
                                        line += str(num[j]) + ' '
                        line += '\n'
                        f.write(line)
                        
        else:
            # Write the images to png files.  Image names have the format event-camera-multiplicity.png
            if not omit_imgs:
                for cam in range(4):
                    cv2.imwrite(write_dir + '%i-%i-%i.png' % (event, cam, multiplicity), np.rot90(event_ims[cam]))
    
    else:
        for bub_num in range(multiplicity[event]):
            index = sum(multiplicity[:event]) + bub_num
            position = pts[index]
            pix = pixel_coords[index]
            
            # Update a dictionary containing (event number, bubble number) as keys and (3D bubble position, pixel coordinates in all cameras) as values
            event_data_dict.update( {(event, bub_num) : (position, pix)} )
            full_data_dict.update( {(event, bub_num) : (position, pix)} )
            
            for cam in range(4):
                cam_pix = pix[cam].copy()
                if np.sum(cam_pix) < 0:
                    continue
                
                # This accounts for the fact that the image is rotated prior to being saved (to match the data)
                cam_pix[1] = resolution[1] - cam_pix[1]
                
                bubble = np.random.choice(all_bubs)
                bubble = cv2.imread(image_dir + bubble, cv2.IMREAD_GRAYSCALE)
                
                bub_scaled = scale_bubble(bubble, position, cam)
                if cam_pix[0] < bs//2:
                    bub_scaled = bub_scaled[bs//2-cam_pix[0]:, :]
                    event_ims[cam, :cam_pix[0]+bs//2, cam_pix[1]-bs//2:cam_pix[1]+bs//2] += bub_scaled
                    continue
                if cam_pix[1] < bs//2:
                    bub_scaled = bub_scaled[:, bs//2-cam_pix[0]:]
                    event_ims[cam, cam_pix[0]-bs//2:cam_pix[0]+bs//2, :cam_pix[1]+bs//2] += bub_scaled
                    continue
                if cam_pix[0] + bs//2 > resolution[0]:
                    bub_scaled = bub_scaled[:resolution[0]+bs//2-cam_pix[0], :]
                    event_ims[cam, cam_pix[0]-bs//2:, cam_pix[1]-bs//2:cam_pix[1]+bs//2] += bub_scaled
                    continue
                if cam_pix[1] + bs//2 > resolution[1]:
                    bub_scaled = bub_scaled[:, :resolution[1]+bs//2-cam_pix[0]]
                    event_ims[cam, cam_pix[0]-bs//2:cam_pix[0]+bs//2, cam_pix[1]-bs//2:] += bub_scaled
                    continue
                    
                event_ims[cam, cam_pix[0]-bs//2:cam_pix[0]+bs//2, cam_pix[1]-bs//2:cam_pix[1]+bs//2] += bub_scaled
        
        if separate_events:
            # Create an event directory to save images and event data
            if not omit_imgs or not omit_dict or not omit_txt:
                os.makedirs(write_dir + '%i-%i' % (event, multiplicity[event]), exist_ok=True)
            
            # Write the images to png files.  Image names have the format event-camera-multiplicity.png
            if not omit_imgs:
                for cam in range(4):
                    cv2.imwrite(write_dir + '%i-%i/%i-%i-%i.png' % (event, multiplicity[event], event, cam, multiplicity[event]), np.rot90(event_ims[cam]))
            
            # Write the event data dictionary to a pickle file for easy access in future python programs
            if not omit_dict:
                with open(write_dir + '%i-%i/data_dict.pkl' % (event, multiplicity[event]), 'wb') as f:
                    pickle.dump(event_data_dict, f)
            
            # Write the bubble 3D positions and pixel coordinates in the event to a text file saved in the event directory
            if not omit_txt:
                with open(write_dir + '%i-%i/data_text.txt' % (event, multiplicity[event]), 'w') as f:
                    for key in event_data_dict:
                        line = r''
                        for i in range(2):
                            for num in event_data_dict[key][i]:
                                if i == 0:
                                    line += str(num) + ' '
                                else:
                                    for j in range(2):
                                        line += str(num[j]) + ' '
                        line += '\n'
                        f.write(line)
                        
        else:
            # Write the images to png files.  Image names have the format event-camera-multiplicity.png
            if not omit_imgs:
                for cam in range(4):
                    cv2.imwrite(write_dir + '%i-%i-%i.png' % (event, cam, multiplicity[event]), np.rot90(event_ims[cam]))


if not separate_events:
    # Write the data dictionary to a pickle file for easy access in future python programs
    if not omit_dict:
        with open(write_dir + 'data_dict.pkl', 'wb') as f:
            pickle.dump(full_data_dict, f)
    
    # Write the bubble 3D positions and pixel coordinates in all events to a text file
    if not omit_txt:
        with open(write_dir + 'data_text.txt', 'w') as f:
            for key in full_data_dict:
                line = str(key[0]) + ' ' + str(key[1]) + ' '
                for i in range(2):
                    for num in full_data_dict[key][i]:
                        if i == 0:
                            line += str(num) + ' '
                        else:
                            for j in range(2):
                                line += str(num[j]) + ' '
                line += '\n'
                f.write(line)


with open(write_dir + 'info.txt', 'w') as f:
    if separate_events:
        text = 'Events are separated in separate directories named "e-m" where e is the event number \n' + \
               'and m is the bubble multiplicity of the event.  The images in the directories are named \n' + \
               '"e-c-m.png" where e and m are as before and c is the camera number.  Each directory also \n' + \
               'contains a text file "data_text.txt".  Each line in this text file corresponds to a \n' + \
               'bubble in the event and contains the 3D coordinates of the bubble in the detector \n' + \
               'followed by the pixel coordinates in camera 0, 1, 2, and 3 in that order.  The pixel \n' + \
               'coordinates are listed as [-1,-1] if the bubble does not appear in the image.  Each \n' + \
               'directory also contains a pickle file "data_dict.pkl" which can be loaded in python via \n' + \
               'the command pickle.load(f) where f is the opened file.  The dictionary contains the \n' + \
               'tuple (event number, bubble number) as a key and the tuple (3D bubble position, camera \n' + \
               'pixel coordinates) as a value.'
        
        f.write(text)

    else:
        text = 'Images are named "e-c-m.png" where e is the event number, c is the camera number of the \n' + \
               'image, and m is the bubble multiplicity of the event.  The text file "data_text.txt" \n' + \
               'contains one line corresponding to each bubble.  Each line gives the event number, bubble \n' + \
               'number, 3D bubble position, and pixel positions in cameras 0, 1, 2, and 3 in that order. \n' + \
               'Pixel coordinates are listed as [-1,-1] if the bubble does not appear in that image.  The \n' + \
               'pickle file "data_dict.pkl" contains a dictionary which can be loaded with the command \n' + \
               'pickle.load(f) where f is the opened file.  The dictionary has the tuple (event number, \n' + \
               'bubble number) as a key and corresponding value (3D bubble position, camera pixel \n' + \
               'coordinates).'
        
        f.write(text)
