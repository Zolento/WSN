import torch
import torch.nn as nn

class Mymodel(nn.Module):
    def __init__(self):
        super(Mymodel, self).__init__()
        emb_dim = 32

    def forward(self, x):  # x= hp of ms; y = hp of pan
        B,L,C = x.size()

        return output
