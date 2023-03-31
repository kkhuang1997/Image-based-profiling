import os
from collections import OrderedDict

import torchvision
from torch.autograd import Variable
from tqdm import tqdm
import numpy as np

from options.test_options import TestOptions
from torch.utils.data import Dataset, DataLoader
from data.image_folder import make_dataset
from data.base_dataset import get_params, get_transform, normalize
import util.util as util
import torch

import pandas as pd
from PIL import Image

opt = TestOptions().parse(save=False)
opt.nThreads = 1   # test code only supports nThreads = 1
opt.batchSize = 1  # test code only supports batchSize = 1
opt.serial_batches = True  # no shuffle
opt.no_flip = True  # no flip

opt.dataroot = "./results/bf2merge_none/test_200/images/"
opt.label_nc = 0
opt.no_instance = True
opt.name = "bf2merge"
opt.resize_or_crop = 'none'
#opt.fineSize = 1024


DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"

## model
class Viability_predictor(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.model = torchvision.models.resnext101_32x8d(pretrained=False)
        # self.model.conv1 = nn.Conv2d(2 * 3, 64, kernel_size=(7, 7), stride=(2, 2), padding=(3, 3), bias=False)
        self.model.fc = torch.nn.Linear(512 * 4, 1)

    def forward(self, x):
        # x = torch.cat([x, y], dim=1)
        x_via = self.model(x)
        return x_via

## dataset
class MapDataset(Dataset):
    def __init__(self, opt):
        self.opt = opt
        ## image
        self.root_dir = opt.dataroot
        self.paths = sorted(make_dataset(self.root_dir))
        ## via
        df_test = pd.read_excel('./datasets/BF_merge_dataset/test_intensity.xls')
        self.rawDataDIR_test = dict(zip(df_test["Label"], df_test["Norm_mean"]))

    def __len__(self):
        return len(self.rawDataDIR_test)

    def __getitem__(self, index):
        ## test image
        path = self.paths[index]
        A = Image.open(path).convert('RGB')
        params = get_params(self.opt, A.size)
        transform = get_transform(self.opt, params, normalize=False)
        A_tensor = transform(A)
        ## via
        image_name = os.path.split(path)[1]
        img_viability = self.rawDataDIR_test[image_name]

        return A_tensor, img_viability

## test
via = Viability_predictor()
via.load_state_dict(torch.load("E:\\pix2pixHD20230105\\datasets\\BF_merge_dataset\\utils_code\\number_regression_demo\\checkpoints\\bf2merge_normed\\90_net_V.pth"))
print("weights have been loaded!")
via.to(DEVICE)

val_dataset = MapDataset(opt = opt)
val_loader = DataLoader(val_dataset, batch_size=opt.batchSize, shuffle=False, num_workers=int(opt.nThreads))

pred_via_list = []
label_via_list = []
if __name__ == "__main__":
    for idx, (x, y) in enumerate(tqdm(val_loader)):
        x = x.to(DEVICE)
        via.eval()
        with torch.no_grad():
            pred_viability = via(x)

        pred_via, label_via = pred_viability.item(), y.item()
        pred_via_list.append(pred_via)
        label_via_list.append(label_via)

    pred_via_array = np.array(pred_via_list)
    label_via_array = np.array(label_via_list)
    np.save("pred_via_array.npy", pred_via_array)
    np.save("label_via_array.npy", label_via_array)




