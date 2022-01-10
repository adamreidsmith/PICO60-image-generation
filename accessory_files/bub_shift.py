import os
import cv2
import numpy as np

im_dir = 'bub_single_24/'
dest_dir = 'bub_single_24_shifted/'
os.makedirs(dest_dir, exist_ok=True)

def top_ind(im):
    if np.sum(im[0]) > 0:
        return 0
        
    for i in range(len(im)):
        if np.sum(im[i]) > 0:
            return i

def shift_up(im, top):
    shift = np.random.randint(low=-1, high=top)

    if shift > -1:
        for i in range(shift+1):
            im = np.delete(im, 0, 0)
        for i in range(shift+1):
            im = np.append(im, np.zeros((1, len(im[0]))), 0)
            
    return im

for im_name in os.listdir(im_dir):
    im = cv2.imread(im_dir + im_name, cv2.IMREAD_GRAYSCALE)
    
    if im.shape != (32, 32):
        continue

    if np.random.random() >= 0.5:
        # shift up
        top = top_ind(im)
        im = shift_up(im, top)
        
    else:
        # shift down
        im = np.flip(im, 0)
        top = top_ind(im)
        im = shift_up(im, top)
        im = np.flip(im, 0)
        
    if np.random.random() >= 0.5:
        # shift left
        im = im.transpose()
        top = top_ind(im)
        im = shift_up(im, top)
        im = im.transpose()
    
    else:
        # shift right
        im = im.transpose()
        im = np.flip(im, 0)
        top = top_ind(im)
        im = shift_up(im, top)
        im = np.flip(im, 0)
        im = im.transpose()

    assert im.shape == (32, 32)
    
    cv2.imwrite(dest_dir + im_name, im)
    
    
                
