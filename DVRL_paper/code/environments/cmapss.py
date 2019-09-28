from gym import spaces
from gym.utils import seeding
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

import gym
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import pandas as pd
from hmmlearn import hmm
from sklearn import preprocessing
from scipy.stats import norm
import pomegranate


SMALL_SIZE = 7
MEDIUM_SIZE = 9
BIGGER_SIZE = 11

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


class CMAPSS(gym.Env):

    def __init__(self,
                 data, cycle_num=0, timestep=0):
        df = pd.read_csv('~/Documents/hitachi/CMAPSS/train_FD001.txt', sep=" ", header=None)
        unique_unit_values = df[0].unique() #Number of units
        data_cycles = []
        for unit_num in unique_unit_values:
            data_cycles.append(df[df[0] == unit_num])
        self.data = data_cycles
        self.cycle_num = cycle_num
        self.timestep = timestep

    def resize(self, pos):
        return pos * self.box_scale
        # return (pos - self.box_mean) * 2 / self.box_size

    def observation(self):
        cycle_num = self.cycle_num
        timestep = self.timestep
        return self.data[cycle_num][timestep]

    def get_reward(self, action, state):
        rewards = np.array([[100, 50, 0, -50],[-50, 0, 50, 100]])
        return rewards[action][state]

    def get_state(obs):
        state = 0
        diff = 16
        for i in range(len(statemean)):
            stateDiff = obs - statemean[i]
            stateDiffVal = np.sqrt(np.mean(stateDiff**2))
            if stateDiffVal < diff:
                diff = stateDiffVal
                state = i
        return state

    def _step(self,action):
        unitData = dataT_cycles[self.cycle_num]
        d = False
        if action == 1:
            self.timestep = 0
        else:
            self.timestep = self.timestep+1
        obsNext = unitData.values[self.timestep]
        if self.timestep >= len(unitData) - 1:
            d = True
        s1 = get_state(obsNext)
        r1 = get_reward(action, s1)
        return s1,r1,d,{}

    def _reset(self):
        self.timestep=0
        return self.data[self.cycle_num][self.timestep]

    def _render(self):
        print(self.data[self.cycle_num][self.timestep])

    def _seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]