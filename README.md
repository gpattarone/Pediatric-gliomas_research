# Code
Useful Code under construction for Image Classification

### Packages
%matplotlib inline
from __future__ import print_function, division

import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim import lr_scheduler
import numpy as np
import torchvision
from torchvision import datasets, models, transforms
import matplotlib.pyplot as plt
import time
import os
import random
import torchvision.transforms.functional as TF
import copy
from resnet import resnet18

plt.ion()   # interactive mode

### Hyperparameters

num_epochs = 200
batch_size = 16
learning_rate = 1e-3

### Neural Network

model = resnet18()
model = model.to(device)

num_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
print('Number of parameters: %d' % num_params)
