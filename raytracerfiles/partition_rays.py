import matlab.engine
import pickle
import numpy as np
import os

'''
This code propagates rays through the PICO-60 geometry. The jar volume is partitioned into 3x3x3 cm boxes, 
and files containing all rays passing through each box are saved in the folder './partitions'.  The 
files can be read using pickle.load, and each contain a dictionary with the ray index as the key and 
the start and end points of the ray as values.

This should be run before running generate_event.py to create the necessary files.  Note that this 
code runs the matlab files FitPICO60Geometry.m and GetRaysAndPixels.m which requires the Matlab Engine 
API to be installed.
'''

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
        scatters, ray_startingpoints, pixels, gen_params = pickle.load(f)
    print('Read file rays.pkl.')
except:
    print('Computing rays.')
    if not 'eng' in locals() and not 'eng' in globals():
        eng = matlab.engine.start_matlab()
    scatters, ray_startingpoints, pixels, gen_params = eng.GetRaysAndPixels(allparams, [170,108], nargout=4)
    with open('rays.pkl', 'wb') as f:
        pickle.dump([scatters, ray_startingpoints, pixels, gen_params], f)

if 'eng' in locals() or 'eng' in globals():
    eng.quit()


r1, r2, r3, cyl_axis = gen_params[0][1], gen_params[0][3], gen_params[0][5], gen_params[0][6:]
# r1 is main inner cylinder radius
# r2 is minor radius of torus defining the inner surface of the knuckle of the jar
# r3 is the radius of the sphere defining the inner surface of the dome of the jar
s = r3 * (r1 - r2)/(r3 - r2)
z = -r2 * np.sqrt(1 - (s/r3)**2)
d = r3 * z * (1/r3 - 1/r2)
cyl_axis = np.array(cyl_axis)

# These are the dimensions of a box bounding the volume of the detector in which we wish to generate points
box = (30., 30., 96.)  # Each dimension should be an integer multiple of 3
# Coordinate of the center of the box
box_center = (0., 0., 38.)
# These should be the same as specified in partition_rays.py
# These need to be tuples of floats, not ints

# Size of the partition boxes in cm
ps = 3

def injar(p, r1, r2, r3, cyl_axis, z, d):
    # Returns true if point p is inside the jar, false otherwise
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


# Creates a dict of all rays which scattered at least once
print('Creating scatter dictionary.')
scatter_dicts = [{}, {}, {}, {}]
for cam in range(len(scatter_dicts)):
    for scatter in scatters:
        index = int(scatter[3])-1
        if int(pixels[index][0])-1 != cam or index < 0:
            continue
        if index in scatter_dicts[cam]:
            scatter_dicts[cam][index].append(list(scatter[:3]))
        else:
            scatter_dicts[cam].update( {index : [list(scatter[:3])]} )
    for i, starting_point in enumerate(ray_startingpoints):
        if i in scatter_dicts[cam]:
            scatter_dicts[cam][i].insert(0, list(starting_point))

del scatters
del ray_startingpoints


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
    

# Dictionaries of the ray segments for each camera that are inside the jar
# The dictionary key is the ray index and the value is the two scatter points on either side of the ray segment
print('Getting rays which intersect the jar.')
rays_in_jar_dict = [{}, {}, {}, {}]
for cam, s_dict in enumerate(scatter_dicts):
    for index in s_dict:
        scats = s_dict[index]
        npts = len(scats)
        ray_seg_in_jar_index = 0
        for s in range(1,npts):
            pt1 = np.array(scats[-s])
            pt2 = np.array(scats[-s-1])
            if injar(0.5*(pt1 + pt2), r1, r2, r3, cyl_axis, z, d):
                ray_seg_in_jar_index = s
                break
        if ray_seg_in_jar_index != 0:
            rays_in_jar_dict[cam].update( {index : [pt2, pt1]} )

del scatter_dicts


def get_intersection(dist1, dist2, p1, p2):
    # Helper function for LBintersection
    if dist1*dist2 >= 0 or dist1 == dist2:
        return None
    return p1 - (p2-p1)*dist1/(dist2-dist1)

def in_box(hit, B1, B2, axis):
    # Helper function for LBintersection
    if axis == 1 and hit[1] > B1[1] and hit[1] < B2[1] and hit[2] > B1[2] and hit[2] < B2[2]:
        return True
    if axis == 2 and hit[0] > B1[0] and hit[0] < B2[0] and hit[2] > B1[2] and hit[2] < B2[2]:
        return True
    if axis == 3 and hit[0] > B1[0] and hit[0] < B2[0] and hit[1] > B1[1] and hit[1] < B2[1]:
        return True
    return False

def LBintersection(L1, L2, B1, B2):
    # Determines if a line, defined by points L1 and L2 intersects a box, defined by minimal and maximal corners B1 and B2.
    if (L1[0] < B1[0] and L2[0] < B1[0]) or (L1[0] > B2[0] and L2[0] > B2[0]):
        return False
    if (L1[1] < B1[1] and L2[1] < B1[1]) or (L1[1] > B2[1] and L2[1] > B2[1]):
        return False
    if (L1[2] < B1[2] and L2[2] < B1[2]) or (L1[2] > B2[2] and L2[2] > B2[2]):
        return False
    
    hit = get_intersection(L1[0]-B1[0], L2[0]-B1[0], L1, L2)
    if hit is not None and in_box(hit, B1, B2, 1):
        return True
    hit = get_intersection(L1[1]-B1[1], L2[1]-B1[1], L1, L2)
    if hit is not None and in_box(hit, B1, B2, 2):
        return True
    hit = get_intersection(L1[2]-B1[2], L2[2]-B1[2], L1, L2)
    if hit is not None and in_box(hit, B1, B2, 3):
        return True
    hit = get_intersection(L1[0]-B2[0], L2[0]-B2[0], L1, L2)
    if hit is not None and in_box(hit, B1, B2, 1):
        return True
    hit = get_intersection(L1[1]-B2[1], L2[1]-B2[1], L1, L2)
    if hit is not None and in_box(hit, B1, B2, 2):
        return True
    hit = get_intersection(L1[2]-B2[2], L2[2]-B2[2], L1, L2)
    if hit is not None and in_box(hit, B1, B2, 3):
        return True
    return False
        

partition_coords = [[i,j,k] for i in range(int(box[0]//ps)) for j in range(int(box[1]//ps)) for k in range(int(box[2]//ps))]
partition_ray_dicts = [{} for _ in range(len(partition_coords))]

# List of pairs of coordinates, minimal and maximal, defining each box in the partition
box_coords = []
for coords in partition_coords:
    B1 = [float(coords[0]*ps - box[0]/2 + box_center[0]),
          float(coords[1]*ps - box[1]/2 + box_center[1]),
          float(coords[2]*ps - box[2]/2 + box_center[2])]
    B2 = [float((coords[0] + 1)*ps - box[0]/2 + box_center[0]),
          float((coords[1] + 1)*ps - box[1]/2 + box_center[1]),
          float((coords[2] + 1)*ps - box[2]/2 + box_center[2])]
    box_coords.append([B1, B2])

# For each camera and each box in the partition, make a dictionary containing all rays passing through that box
print('Getting rays in volume partitions.')
os.makedirs('partitions', exist_ok=True)
for cam, ray_dict in enumerate(rays_in_jar_dict):
    for i, [B1, B2] in enumerate(box_coords):
        rays_in_box_dict = {}
        for ray_index in ray_dict:
            L1, L2 = ray_dict[ray_index]
            if LBintersection(L1, L2, B1, B2):
                rays_in_box_dict.update( {ray_index : [L1,L2]} )
        with open('partitions/cam%ibox%i.pkl' % (cam,i), 'wb') as f:
            pickle.dump(rays_in_box_dict, f)
        print('Wrote file partitions/cam%ibox%i.pkl.' % (cam,i))
