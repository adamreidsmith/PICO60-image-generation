import torch as t
from torch import nn
from torch.utils.data import Dataset, DataLoader
from torch.optim import Adam

import numpy as np
import cv2
import os

batch_size = 16
img_size = 32
latent_dim = 64
channels = 1
lr = 0.00001 #0.0002
b1, b2 = 0.5, 0.999
n_epochs = 50
n_ch = 64  # Number of channels in initial convolution layers
sample_interval = 200  # Save a generated image every sample_interval number of batches

write_dir = '../images/gan_conv_test/'
data_dir = '../data/bub_single_24new/'
os.makedirs(write_dir, exist_ok=True)

Zeros = t.zeros((batch_size, 1, img_size, img_size))
Ones = t.ones((batch_size, 1, img_size, img_size))

class Data(Dataset):
    def __init__(self):
        self.data = [cv2.imread(data_dir + im, cv2.IMREAD_GRAYSCALE)//255 for im in os.listdir(data_dir) if '-0.png' in im]
        self.data = t.Tensor([im for im in self.data if im.shape == (img_size, img_size)])
        
        self.len = self.data.shape[0]

    def __len__(self):
        return self.len

    def __getitem__(self, index):
        return self.data[index]

dataset = Data()

dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

img_shape = (channels, img_size, img_size)

class Generator(nn.Module):
    # Takes as input a random tensor of size (batch_size, latent_dim) generated by noise()
    '''def __init__(self):
        super(Generator, self).__init__()
        
        self.lin = nn.Sequential(nn.Linear(latent_dim, n_ch * (img_size//2)**2))
        # Reshape to (N, n_ch, im_size//2, im_size//2)
        self.conv_block = nn.Sequential(
            nn.BatchNorm2d(n_ch),
            nn.Conv2d(n_ch, n_ch//2, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(n_ch//2),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Upsample(scale_factor=2),
            nn.Conv2d(n_ch//2, n_ch//4, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(n_ch//4),
            nn.LeakyReLU(0.2, inplace=True),
            #nn.Upsample(scale_factor=2),
            nn.Conv2d(n_ch//4, channels, kernel_size=3, stride=1, padding=1),
            nn.Sigmoid()
            )
    
    def forward(self, im):
        im = self.lin(im)
        im = im.view(im.shape[0], n_ch, img_size//2, img_size//2)
        return self.conv_block(im)
            
    # (N, latent_dim) ->[linear] (N, n * 16**2) ->[reshape] (N, n, 16, 16) ->[conv] (N, n/2, 16, 16) ->[upsampling] (N, n/2, 32, 32) ->[conv] (N, n/4, 32, 32) ->[conv] (N, 1, 32, 32)'''
    
    
    def __init__(self):
        super(Generator, self).__init__()
        
        self.init_size = img_size//4
        self.linear = nn.Sequential(nn.Linear(latent_dim, n_ch * self.init_size**2))  #Outsize=(8x8xn) with n=128
        
        self.conv_blocks = nn.Sequential(
            nn.BatchNorm2d(n_ch),
            nn.Upsample(scale_factor=2),
            nn.Conv2d(n_ch, n_ch, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(n_ch, 0.8),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Upsample(scale_factor=2),
            nn.Conv2d(n_ch, n_ch//2, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(n_ch//2, 0.8),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Conv2d(n_ch//2, channels, kernel_size=3, stride=1, padding=1),
            nn.Sigmoid()
            )
    
    def forward(self, noise):
        generated_im = self.linear(noise)
        generated_im = generated_im.view(generated_im.shape[0], n_ch, self.init_size, self.init_size)
        generated_im = self.conv_blocks(generated_im)
        #generated_im = t.where(generated_im > 0.5, Ones, generated_im)
        return generated_im
        
    # (N, latent_dim) ->[linear] (N, n * 8**2) ->[reshape] (N, n, 8, 8) ->[upsampling] (N, n, 16, 16) ->[conv] (N, n, 16, 16) ->[upsampling] (N, n, 32, 32) ->[conv] (N, n/2, 32, 32) ->[conv] (N, 1, 32, 32)
            

class Discriminator(nn.Module):
    # Takes as input a data image or image generated by Generator
    def __init__(self):
        super(Discriminator, self).__init__()
        
        self.model = nn.Sequential(
            nn.Linear(int(np.prod(img_shape)), 8*latent_dim),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Linear(8*latent_dim, 4*latent_dim),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Linear(4*latent_dim, 1),
            nn.Sigmoid(),
        )
    
    def forward(self, img):
        img_flat = img.view(img.size(0), -1)
        validity = self.model(img_flat)
        return validity
        

def generate_noise():
    return t.rand((batch_size, latent_dim))

# Loss function
adversarial_loss = nn.BCELoss()

G = Generator()
D = Discriminator()

# Optimizers
optimizer_G = Adam(G.parameters(), lr=lr, betas=(b1, b2))
optimizer_D = Adam(D.parameters(), lr=lr, betas=(b1, b2))

for epoch in range(n_epochs):
    for i, data_imgs in enumerate(dataloader):
        
        # Configure generator input
        data_shape = data_imgs.shape
        data_imgs = data_imgs.view(data_imgs.shape[0], channels, data_imgs.shape[1], data_imgs.shape[2])
        
        # Adversarial ground truths
        valid = t.ones((data_shape[0], 1), requires_grad=False)
        fake = t.zeros((data_shape[0], 1), requires_grad=False)
        
        # ----------------
        # Train Generator
        # ----------------
        
        optimizer_G.zero_grad()
        
        # Sample noise as generator input
        noise = generate_noise()
        
        # Generate a batch of images
        gen_imgs = G(noise)
        gen_imgs = gen_imgs[:data_shape[0]]
        
        # Loss measures generator's ability to fool the discriminator
        G_loss = adversarial_loss(D(gen_imgs), valid)

        G_loss.backward()
        optimizer_G.step()
        
        # --------------------
        # Train Discriminator
        # --------------------
        
        optimizer_D.zero_grad()
        
        # Measure discriminator's ability to classify real from generated samples
        real_loss = adversarial_loss(D(data_imgs), valid)
        fake_loss = adversarial_loss(D(gen_imgs.detach()), fake)
        D_loss = (real_loss + fake_loss) / 2
        
        D_loss.backward()
        optimizer_D.step()
        
        # ----------------
        # Log Progress
        # ----------------
        
        print(
            "[Epoch %d/%d] [Batch %d/%d] [D loss: %f] [G loss: %f]"
            % (epoch, n_epochs, i, len(dataloader), D_loss.item(), G_loss.item())
        )

        batches_done = epoch*len(dataloader) + i
        if batches_done % sample_interval == 0:
            os.makedirs(write_dir + str(batches_done), exist_ok=True)
            for j, gen_im in enumerate(gen_imgs):
                cv2.imwrite('%s%d/%d.png' % (write_dir, batches_done, j), gen_im[0].detach().numpy()*255)
