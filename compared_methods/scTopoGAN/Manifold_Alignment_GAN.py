"""Generative Adversarial Network while projecting source into target space"""

# Import modules
import numpy as np
import torch
from torch import nn, optim

from torch.autograd.variable import Variable
import matplotlib.pyplot as plt

import torch
if torch.cuda.is_available():
    device = torch.device("cuda:1")
else:
    device ="cpu"

# Define function to initialise weights
def weights_init(m):
    if isinstance(m, nn.Linear):
        #nn.init.xavier_uniform_(m.weight)
        torch.nn.init.normal_(m.weight, mean=0, std=0.02)
        
def sample_data(data, batch_size):
    """
    sample data from numpy array datatype
    """
    return data[np.random.choice(data.shape[0], size=batch_size, replace=False),:]

class GeneratorNet(torch.nn.Module):
    """
    A one hidden-layer generative neural network
    """

    def __init__(self, input_dim, output_dim):
        super(GeneratorNet, self).__init__()
        n_features = input_dim
        n_out = output_dim

        self.hidden0 = nn.Sequential(
            nn.Linear(n_features, 30), nn.ReLU())

        self.out = nn.Sequential(
            nn.Linear(30, n_out))
        
        self.hidden0.apply(weights_init)
        self.out.apply(weights_init)

    def forward(self, x):
        x = self.hidden0(x)
        x = self.out(x)
        return x

class DiscriminatorNet(torch.nn.Module):
    """
    A two hidden-layer discriminative neural network
    """

    def __init__(self, input_dim):
        super(DiscriminatorNet, self).__init__()
        n_features = input_dim
        n_out = 1
	
	# To add dropout in layers, type nn.Dropout(0.3)

        self.hidden0 = nn.Sequential(
            nn.Linear(n_features, 60), 
            nn.LeakyReLU(0.2), nn.Dropout(0.3))

        self.hidden1 = nn.Sequential(
            nn.Linear(60, 30),
            nn.LeakyReLU(0.2), nn.Dropout(0.3))

        self.out = nn.Sequential(
            torch.nn.Linear(30, n_out),
            torch.nn.Sigmoid())

        # Initialise all weights
        
        self.hidden0.apply(weights_init)
        self.hidden1.apply(weights_init)
        self.out.apply(weights_init)
        
    def forward(self, x):
        x = self.hidden0(x)
        x = self.hidden1(x)
        x = self.out(x)
        return x

# Define functions to generate arrays of 0s or 1s
def ones_target(size):
    '''
    Tensor containing ones, with shape = size
    '''
    data = Variable(torch.ones(size, 1))
    return data

def zeros_target(size):
    '''
    Tensor containing zeros, with shape = size
    '''
    data = Variable(torch.zeros(size, 1))
    return data

# Define function to train the discriminator
def train_discriminator(optimizer, discriminator, real_data, fake_data):
    # Define loss function
    loss = nn.BCELoss()

    N = real_data.size(0)
    # Reset gradients
    optimizer.zero_grad()

    # 1.1 Train on Real Data
    prediction_real = discriminator(real_data)
    # Calculate error and backpropagate
    error_real = loss(prediction_real, ones_target(N).to(device))
    error_real.backward()

    # 1.2 Train on Fake Data
    prediction_fake = discriminator(fake_data)
    # Calculate error and backpropagate
    error_fake = loss(prediction_fake, zeros_target(N).to(device))
    error_fake.backward()

    # 1.3 Update weights with gradients
    optimizer.step()

    # Return error and predictions for real and fake inputs
    return error_real + error_fake

# Define function to train the generator
def train_generator(optimizer, discriminator, fake_data):
    # Define loss function
    loss = nn.BCELoss()
    N = fake_data.size(0)
    # Reset gradients
    optimizer.zero_grad()
    # Sample noise and generate fake data
    prediction = discriminator(fake_data)
    #print(prediction.is_cuda)
    # Calculate error and backpropagate
    error = loss(prediction, ones_target(N).to(device))
    error.backward()
    # Update weights with gradients
    optimizer.step()
    # Return error
    return error

def train(generator, discriminator, batch_size, source_tech, target_tech, num_epochs, g_learning_rate,
          d_learning_rate, checkpoint_epoch, techs, path_prefix, path_suffix):
    # Define optimisers for Discriminator and Generator
    d_optimizer = optim.Adam(discriminator.parameters(), lr=d_learning_rate, betas=(0.5, 0.999))
    g_optimizer = optim.Adam(generator.parameters(), lr=g_learning_rate, betas=(0.5, 0.999))
    D_Losses = []
    G_Losses = []
    torch.set_num_threads(2)
    for epoch in range(num_epochs):
        d_loss = 0
        g_loss = 0
        print("Current Epoch: ", epoch)
        
        for _ in range(1,  source_tech.shape[0] // batch_size + 1):
            source_sample = sample_data(source_tech, batch_size)
            target_sample = sample_data(target_tech, batch_size)
            # Convert source and target from eager tensor to native tensor
            source = torch.tensor(source_sample).float().to(device)

            target = torch.tensor(target_sample).float().to(device)

            # source = torch.tensor(source_sample).float()
            #
            # target = torch.tensor(target_sample).float()

            # 1. Train Discriminator
            real_data = target
            # Generate fake data and detach
            # (so gradients are not calculated for generator)
            fake_data = generator(source).detach()
            # Train D
            d_error = train_discriminator(d_optimizer, discriminator, real_data, fake_data)
            d_loss += d_error.item()
            # 2. Train Generator
            source_sample = sample_data(source_tech, batch_size)
            # Convert source and target from eager tensor to native tensor
            source = torch.tensor(source_sample).float().to(device)
            # Generate fake data
            fake_data = generator(source)

            # Train G
            g_error = train_generator(g_optimizer, discriminator, fake_data)
            g_loss += g_error.item()

        d_loss = d_loss * batch_size / source_tech.shape[0]
        g_loss = g_loss * batch_size / source_tech.shape[0]

        D_Losses.append(d_loss)
        G_Losses.append(g_loss)
        
        if ((epoch % checkpoint_epoch) == 0):
            intermediate_path = "{}/Models/{}_to_{}_Generator_{}_{}.pt".format(
                path_prefix, techs[0], techs[1], epoch, path_suffix)
            torch.save(generator, intermediate_path)

    # x = range(epoch + 1)
    # plt.plot(x, D_Losses, label="Discriminator loss")
    # plt.plot(x, G_Losses, label="Generator loss")
    # plt.xlabel("Epoch")
    # plt.ylabel("Loss")
    # plt.legend()
    # plt.show()

    return generator
