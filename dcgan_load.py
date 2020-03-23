import torch as t
from torch import nn
from torchvision.utils import save_image
import cv2
import os
import argparse

from dcgan import *

'''
Load and run the model created by dcgan.py to generate artificial images.  The model is trained
on Gaussian inputs with standard deviation 1, but it seems best to run the model with data with
standard deviation ~0.8.
'''

parser = argparse.ArgumentParser()
parser.add_argument("--n_images", type=int, default=100, help="number of images to generate")
parser.add_argument("--stdev", type=float, default=0.8, help="standard deviation of Gaussian input data")
parser.add_argument("--thresholding", type=int, default=1, help="turn on or off thresholding")
parser.add_argument("--threshold", type=float, default=0.95, help="threshold in [0,1] to apply")
opt = parser.parse_args()

# Directory to write genereated images to
write_dir = './dcgan_images/'

model = 'dcgan_model.pt'

# ---------------
# Generate images
# ---------------
        
os.makedirs(write_dir, exist_ok=True)

G = t.load(model)
G.eval()

for i in range(opt.n_images):
    input = t.normal(mean=0, std=opt.stdev, size=(1, latent_dim))
    img = G(input)
    if opt.thresholding:
        img = t.where(img > opt.threshold, t.ones_like(img), t.zeros_like(img))
    save_image(img*255, 'dcgan_images/%i.png' % i)
