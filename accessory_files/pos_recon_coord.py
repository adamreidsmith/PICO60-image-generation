import torch as t
from torch import nn
from torch.utils.data import Dataset, DataLoader, random_split, SubsetRandomSampler
from torch.optim import Adam

import numpy as np
import pickle


batch_size = 16
n_epochs = 300
lr = 0.0002 #0.0001
b1, b2 = 0.9, 0.999  # Betas for Adam optimization function
train_frac = 0.8  # Train on this fraction of the data, evaluate on the rest

real_data_file = './data/single_coord_pos_24.txt'
fake_data_file = './artificial_events/singles/data_dict.pkl'


def normalize_pos(pos):
    # Normalize position from [-15,15]x[-15,15]x[-10,85] into [-1,1]^3.
    if len(pos.size()) == 1:
        return t.Tensor([pos[0]/15, pos[1]/15, (pos[2] - 37.5)/47.5])
    d = len(pos)
    return t.cat(((pos[:,0]/15).view(d, 1), (pos[:,1]/15).view(d, 1), ((pos[:,2] - 37.5)/47.5).view(d, 1)), 1)
    
def denormalize_pos(pos):
    # Denormalize position from [-1,1]^3 into [-15,15]x[-15,15]x[-10,85]
    if len(pos.size()) == 1:
        return t.Tensor([pos[0]*15, pos[1]*15, pos[2]*47.5 + 37.5])
    d = len(pos)
    return t.cat(((pos[:,0]*15).view(d, 1), (pos[:,1]*15).view(d, 1), (pos[:,2]*47.5 + 37.5).view(d, 1)), 1)


class Data(Dataset):
    def __init__(self):
        
        # Load the real images as data
        with open(real_data_file, 'r') as f:
            lines = f.readlines()
        
        # 1124 real data points
        self.real_data, self.real_pos = {}, {}
        for line in lines:
            items = line.split()
            coords = t.Tensor( ((float(items[2]), float(items[3])),
                                (float(items[4]), float(items[5])),
                                (float(items[6]), float(items[7])),
                                (float(items[8]), float(items[9]))) )
            self.real_data.update( {(items[0], int(items[1])) : coords} )
            
            pos_3d = t.Tensor( (float(items[10])/10, float(items[11])/10, float(items[12])/10) )
            self.real_pos.update( {(items[0], int(items[1])) : normalize_pos(pos_3d)} )
        
#        for key in self.real_data:
#            for i, coord in enumerate(self.real_data[key]):
#                if t.sum(coord) < 0:
#                    self.real_data[key][i] = t.Tensor([-1,-1])
        
        # Load the artificial images as data
        with open(fake_data_file, 'rb') as f:
            fake_data_dict = pickle.load(f)
        
        self.fake_data, self.fake_pos = {}, {}
        for key in fake_data_dict:
            coords = t.Tensor(fake_data_dict[key][1])
            self.fake_data.update( {key : coords} )
            
            pos_3d = t.Tensor(fake_data_dict[key][0])
            self.fake_pos.update( {key : pos_3d} )
        
        for key in self.fake_data:
            for i, coord in enumerate(self.fake_data[key]):
                if t.sum(coord) < 0:
                    self.fake_data[key][i] = t.Tensor([-100, -100])
        
        self.all_data = {**self.real_data, **self.fake_data}
        self.all_pos = {**self.real_pos, **self.fake_pos}
        
        self.len = len(self.all_data)
        self.run_ev = list(self.all_data.keys())
        
        '''
        Real data has [-100, -100] as pixel coordiates for no bubble in an image, but artificial data has [-1, -1].
        '''
        
        
    def __len__(self):
        return self.len

    def __getitem__(self, index):
        key = self.run_ev[index]
        return self.all_data[key].flatten(), self.all_pos[key]

## Random split of data
#dataset = Data()
#
#train_len = int(train_frac * dataset.len)
#valid_len = dataset.len - train_len
#
#train_data, valid_data = random_split(dataset, (train_len, valid_len))
#
#train_loader = DataLoader(dataset=train_data, batch_size=batch_size, shuffle=True)
#valid_loader = DataLoader(dataset=valid_data, batch_size=batch_size, shuffle=True)

# Custom split of data
# Validate on half the real data, train on the rest
dataset = Data()
indices = list(range(dataset.len))
#split = len(dataset.real_data)
#valid_indices, train_indices = indices[:split], indices[split:]

split = 2124
valid_indices, train_indices = indices[1124:split], indices[split:4124]

valid_sampler = SubsetRandomSampler(valid_indices)
train_sampler = SubsetRandomSampler(train_indices)

valid_loader = DataLoader(dataset, batch_size=batch_size, sampler=valid_sampler)
train_loader = DataLoader(dataset, batch_size=batch_size, sampler=train_sampler)


class Swish(nn.Module):
    def __init__(self):
        super(Swish, self).__init__()
    
    def forward(self, x):
        return t.mul(t.sigmoid(x), x)


class Model(nn.Module):

    def __init__(self):
        super(Model, self).__init__()
        
        '''
        Original paper claims batch norm should be implemented after activation, but sources say it works better after.
        See: https://stackoverflow.com/questions/39691902/ordering-of-batch-normalization-and-dropout
        
        -> FC -> ReLu -> Dropout -> BatchNorm -> FC
        
        OR
        
        -> FC -> BatchNorm -> ReLu -> Dropout -> FC
        '''
        
        a, b, c, d = 8, 100, 30, 3
        self.layers = nn.Sequential(nn.Linear(a, b),
                                    nn.BatchNorm1d(b),
                                    #nn.LeakyReLU(0.1, inplace=True),
                                    Swish(),
                                    #nn.Dropout(p=0.4),
                                    
                                    nn.Linear(b, c),
                                    nn.BatchNorm1d(c),
                                    #nn.LeakyReLU(0.1, inplace=True),
                                    Swish(),
                                    #nn.Dropout(p=0.4),
                                    
                                    nn.Linear(c, d),
                                    nn.Tanh()  # Doesn't do operation inplace, but this is fine since it occurs last
                                   )
                                   
        a, b, c, d = 8, 100, 50, 3
        self.layers2 = nn.Sequential(nn.Linear(a, b),
                                    nn.ReLU(inplace=True),
                                    nn.Dropout(p=0.4),
                                    nn.BatchNorm1d(b),
                                    
                                    nn.Linear(b, c),
                                    nn.ReLU(inplace=True),
                                    nn.Dropout(p=0.4),
                                    nn.BatchNorm1d(c),
                                    
                                    nn.Linear(c, d),
                                    nn.Tanh()  # Doesn't do operation inplace, but this is fine since it occurs last
                                   )

    def forward(self, coords):
        return self.layers(coords)


model = Model()

loss_func = nn.MSELoss()

optimizer = Adam(model.parameters(), lr=lr, betas=(b1, b2))


def train():
    model.train()
    train_loss = 0
    dists = t.empty(0)
    for data in train_loader:
        coords, pos = data

        # Reset gradients to zero
        optimizer.zero_grad()

        # Generate a predicted 3d position
        prediction = model(coords)

        # Compute the loss
        loss = loss_func(denormalize_pos(prediction), pos)

        # Update weights
        loss.backward()
        optimizer.step()

        train_loss += loss.detach().item()
                
        diff = denormalize_pos(prediction) - pos
        dist = t.norm(diff, dim=1)
        dists = t.cat((dists, dist))

    return train_loss/len(train_indices), dists.mean(), dists.std()


def evaluate():
    model.eval()
    valid_loss = 0
    dists = t.empty(0)
    for data in valid_loader:
        coords, pos = data

        # Generate a predicted 3d position
        prediction = model(coords)

        # Compute the loss
        loss = loss_func(denormalize_pos(prediction), pos)

        valid_loss += loss.detach().item()
        
        diff = denormalize_pos(prediction) - pos
        dist = t.norm(diff, dim=1)
        dists = t.cat((dists, dist))
        
    return valid_loss/len(valid_indices), dists.mean(), dists.std()
 
 
for epoch in range(n_epochs):

    # Train the network
    train_loss, train_dist, train_std = train()

    # Evaluate the network
    valid_loss, valid_dist, valid_std = evaluate()

    print(
          '[Epoch %d/%d] [Train Loss: %f] [Valid Loss: %f] [Train Dist: %f ± %f] [Valid Dist: %f ± %f]'
          % (epoch, n_epochs, train_loss, valid_loss, train_dist, train_std, valid_dist, valid_std)
         )
