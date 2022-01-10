import torch as t
from torch import nn
from torch.utils.data import Dataset, DataLoader, random_split
from torch.optim import Adam

import numpy as np
import cv2
import os


batch_size = 16
n_epochs = 1000
img_size = np.array((1088, 1700))  # Size of the incoming images
channels = 4  # Number of images in each event
lr = 0.00001 #0.0001
b1, b2 = 0.9, 0.999  # Betas for Adam optimization function
n_ch = 128  # Number of channels in initial convolution layers
sample_interval = 30  # Save a generated image every sample_interval number of batches
normalization_vals = t.Tensor([150, 150, 550])  # Values to normalize the positions
train_frac = 0.9  # Train on this fraction of the data, evaluate on the rest

#write_dir = './write_dir'
data_dir = './data/'
im_folder = 'full_single_24/'
pos_file = 'single_pos_24.txt'
#model_path = './model_path'
#os.makedirs(write_dir, exist_ok=True)

class Data(Dataset):
    def __init__(self):
        '''
        self:
            run_ev      List of tuple's (run, event)
            data        Dictionary of all 4 images in an event with tuple's (run, event) as keys
            pos_dict    Dictionary of 3d positions of bubble in event with tuple's (run, event) as keys
            len         Total nimber of events
            
        Images must be named:  run-event-camera.png
            run is 10 character run number
            event is one ot two digit run number
            camera is one digit camera number
        '''
        
        im_dir = data_dir + im_folder
    
        # --------------------------------
        # Load event images in groups of 4
        # --------------------------------

        def sort_key(im_name):
            return int(im_name[-5])
        
        # Create a list containing sublists of the image names for each event
        ims_by_event = []
        im_list = os.listdir(im_dir)
        run_ev_set = set([im[:-6] for im in im_list])  # Using a set makes sure elements are unique
        for run_ev in run_ev_set:
            event_ims = []
            for im in im_list:
                if im[:-6] == run_ev:
                    event_ims.append(im)
            if 0 < len(event_ims) < 5:
                event_ims = sorted(event_ims, key=sort_key)
                ims_by_event.append(event_ims)
        
        self.len = len(ims_by_event)
        
        # For each sublist in ims_by_event, append 'blank' if the camera did not image the bubble in that event
        for j, event_ims in enumerate(ims_by_event):
            if len(event_ims) == 4:
                continue
            cams_with_ims = [int(im_name[-5]) for im_name in event_ims]
            for k in range(4):
                if not k in cams_with_ims:
                    ims_by_event[j].insert(k, 'blank')
        
        # Get a list of image names corresponding to events in ims_by_event
        self.run_ev = []
        for event in ims_by_event:
            for im_name in event:
                if im_name != 'blank':
                    run = im_name[:10]
                    ev = int(im_name[11:-6])
                    self.run_ev.append((run,ev))
                    break

        # Load images. Load a blank image where 'blank' is found
        self.data = t.empty((self.len, 4, *img_size))
        blank = t.zeros([*img_size])
        for i, event in enumerate(ims_by_event):
            loaded_ims_for_current_event = t.empty((4, *img_size))
            for j, im_name in enumerate(event):
                # Load a blank image if the image name is 'blank'
                if im_name == 'blank':
                    loaded_ims_for_current_event[j] = blank
                # Otherwise, load the image
                else:
                    im = t.from_numpy(cv2.imread(im_dir + im_name, cv2.IMREAD_GRAYSCALE)//255)
                    # Stop if image is not the required dimension
                    loaded_ims_for_current_event[j] = im
            # Update the data dictionary
            self.data[i] = loaded_ims_for_current_event
            if i % 100 == 0:
                print('Loaded images from event %i of %i' % (i, len(ims_by_event)))
            
        # ---------------------
        # Load bubble 3D positions
        # ---------------------
        
        f = open(data_dir + pos_file, 'r')
        self.pos_dict = {}
        for line in f:
            ld = line[:-2].split(' ')
            self.pos_dict.update( { (ld[0], int(ld[1])) : t.Tensor((float(ld[2]), float(ld[3]), float(ld[4]))) } )
        f.close()
                
        for key in self.pos_dict:
            self.pos_dict[key] /= normalization_vals

    __len__ = lambda self: self.len

    def __getitem__(self, index):
        key = self.run_ev[index]
        return self.data[index], self.pos_dict[key]
        
dataset = Data()

train_len = int(train_frac * dataset.len)
valid_len = dataset.len - train_len

train_data, valid_data = random_split(dataset, (train_len, valid_len))

train_loader = DataLoader(dataset=train_data, batch_size=batch_size, shuffle=True)
valid_loader = DataLoader(dataset=valid_data, batch_size=batch_size, shuffle=True)

class Model(nn.Module):

    def __init__(self):
        super(Model, self).__init__()
        
        mp_dim = lambda ds, ks, s=1, p=0, d=1: int((ds + 2*p - d*(ks - 1) - 1 + s)//s)
        
        def convolution_block(in_filters, out_filters, bn=True):
            block = [
                nn.Conv2d(in_filters, out_filters, kernel_size=3, stride=1, padding=1),
                nn.LeakyReLU(0.2, inplace=True),
                nn.MaxPool2d(kernel_size=2, stride=2),
                nn.Dropout2d(0.25)
            ]
            if bn:
                block.append(nn.BatchNorm2d(out_filters, eps=0.1))
            return block
        
        self.conv_blocks = nn.Sequential(
            *convolution_block(channels, 12, bn=False),
            *convolution_block(12, 24),
            *convolution_block(24, 48)
            # Maybe add some linear layers here
        )
        d_size = img_size // 2**3
        
        self.lin_layers = nn.Sequential(
            nn.Linear(48*np.prod(d_size), 100),
            nn.ReLU(inplace=True),
            nn.Linear(100, 3),
            nn.Tanh()
        )
    
    def forward(self, imgs):
        out = self.conv_blocks(imgs)
        out = out.view(batch_size, -1)
        pos = self.lin_layers(out)
        return pos

model = Model()

loss_func = nn.MSELoss()

optimizer = Adam(model.parameters(), lr=lr, betas=(b1, b2))

def train():
    model.train()
    train_loss = 0
    for data in train_loader:
        ims, pos = data

        # Reset gradients to zero
        optimizer.zero_grad()

        # Generate a predicted 3d position
        prediction = model(ims)

        # Compute the loss
        loss = loss_func(prediction, pos)

        # Update weights
        loss.backward()
        optimizer.step()

        train_loss += loss.detach().item()

        prediction_corrected = prediction.detach() * normalization_vals
        pos_corrected = pos.detach() * normalization_vals
        dist = t.norm(prediction_corrected - pos_corrected, dim=1)
        mean_dist = t.mean(dist)
 
    return train_loss/(len(train_loader)*batch_size), mean_dist

def evaluate():
    model.eval()
    valid_loss = 0
    for data in valid_loader:
        ims, pos = data

        # Generate a predicted 3d position
        prediction = model(ims)

        # Compute the loss
        loss = loss_func(prediction, pos)

        valid_loss += loss.detach().item()

        prediction_corrected = prediction.detach() * normalization_vals
        pos_corrected = pos.detach() * normalization_vals
        dist = t.norm(prediction_corrected - pos_corrected, dim=1)
        mean_dist = t.mean(dist)
 
    return valid_loss/(len(valid_loader)*batch_size), mean_dist
 
 
for epoch in range(n_epochs):

    # Train the network
    train_loss, train_dist = train()

    # Evaluate the network
    valid_loss, valid_dist = evaluate()

    print(
     '[Epoch %d/%d] [Train Loss: %f] [Valid Loss: %f] [Train Distance: %f] [Valid Distance: %f]'
     % (epoch, n_epochs, train_loss, valid_loss, train_dist, valid_dist)
    )
    
        
        
        
        
'''
class Model(nn.Module):

    def __init__(self):
        super(Model, self).__init__()
        
        def convolution_block(in_filters, out_filters, bn=True):
            block = [nn.Conv2d(in_filters, out_filters, kernel_size=3, stride=2, padding=1),
                     nn.LeakyReLU(0.2, inplace=True),
                     nn.Dropout2d(0.25)]
            if bn:
                block.append(nn.BatchNorm2d(out_filters, 0.8))
            return block

        self.model = nn.Sequential(
            *convolution_block(channels, n_ch//8, bn=False),
            *convolution_block(n_ch//8, n_ch//4),
            *convolution_block(n_ch//4, n_ch//2),
            *convolution_block(n_ch//2, n_ch)
        )
        
        # The height and width of downsampled image
        ds_size = img_size // 2 ** 4
        self.lin_layer = nn.Sequential(nn.Linear(n_ch * ds_size ** 2, 3), nn.Tanh())
        
    def forward(self, img):
        out = self.model(img)
        out = out.view(out.shape[0], -1)
        pos = self.lin_layer(out)
        return pos


model = Model()

loss_func = nn.MSELoss()

optimizer = Adam(model.parameters(), lr=lr, betas=(b1, b2))

def train():
    model.train()
    train_loss = 0
    for data in train_loader:
        ims, pos = data
        
        # Reset gradients to zero
        optimizer.zero_grad()
        
        # Generate a predicted 3d position
        prediction = model(ims)
        
        # Compute the loss
        loss = loss_func(prediction, pos)
        
        # Update weights
        loss.backward()
        optimizer.step()
        
        train_loss += loss.detach().item()
        
        prediction_corrected = prediction.detach() * normalization_vals
        pos_corrected = pos.detach() * normalization_vals
        dist = t.norm(prediction_corrected - pos_corrected, dim=1)
        mean_dist = t.mean(dist)
    
    return train_loss/(len(train_loader)*batch_size), mean_dist

def evaluate():
    model.eval()
    valid_loss = 0
    for data in valid_loader:
        ims, pos = data
        
        # Generate a predicted 3d position
        prediction = model(ims)
        
        # Compute the loss
        loss = loss_func(prediction, pos)
        
        valid_loss += loss.detach().item()
        
        prediction_corrected = prediction.detach() * normalization_vals
        pos_corrected = pos.detach() * normalization_vals
        dist = t.norm(prediction_corrected - pos_corrected, dim=1)
        mean_dist = t.mean(dist)
    
    return valid_loss/(len(valid_loader)*batch_size), mean_dist
    
    
for epoch in range(n_epochs):

    # Train the network
    train_loss, train_dist = train()
    
    # Evaluate the network
    valid_loss, valid_dist = evaluate()
    
    print(
        '[Epoch %d/%d] [Train Loss: %f] [Valid Loss: %f] [Train Distance: %f] [Valid Distance: %f]'
        % (epoch, n_epochs, train_loss, valid_loss, train_dist, valid_dist)
    )
    
'''
   
'''
def generate_noise(name):
    if name == 'Gaussian':
        return t.normal(mean=0, std=1, size=(batch_size, latent_dim))
    return t.rand((batch_size, latent_dim))
    
def save_model():
    print('Saving Saving generator model to ' + model_path)
    t.save(G, model_path)

# Loss function
adversarial_loss = nn.BCELoss()

G = Generator()
D = Discriminator()

# Initialize weights
#G.apply(weights_init_normal)
#D.apply(weights_init_normal)

# Optimizers
optimizer_G = Adam(G.parameters(), lr=lr, betas=(b1, b2))
optimizer_D = Adam(D.parameters(), lr=lr, betas=(b1, b2))
try:
    for epoch in range(n_epochs):
        for i, data_imgs in enumerate(dataloader):
            
            # Configure generator input
            data_shape = data_imgs.shape
            data_imgs = data_imgs.view(data_imgs.shape[0], channels, data_imgs.shape[-2], data_imgs.shape[-1])
            
            # Adversarial ground truths
            valid = t.ones((data_shape[0], 1), requires_grad=False)
            fake = t.zeros((data_shape[0], 1), requires_grad=False)
            
            # ----------------
            # Train Generator
            # ----------------
            
            optimizer_G.zero_grad()
            
            # Sample noise as generator input
            noise = generate_noise('Gaussian')
            
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
                    for k, im in enumerate(gen_im):
                        cv2.imwrite('%s%d/%d_%d.png' % (write_dir, batches_done, j, k), im.detach().numpy()*255)
    
    save_model()
            
except KeyboardInterrupt:
    save_model()
'''
