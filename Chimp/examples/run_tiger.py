'''
File to initialize training.
Contains settings, network definition for Chainer.
Creates the simulator, replay memory, DQN learner, and passes these to the agent framework for training.
'''

import numpy as np

import chainer
import chainer.functions as F
import chainer.links as L
from chainer import cuda, Function, gradient_check, Variable, optimizers, serializers, utils
from chainer import Link, Chain, ChainList
from memories import ReplayMemoryHDF5

from learners import Learner
from agents import DQNAgent

from simulators.pomdp import POMDPSimulator
from simulators.pomdp import TigerPOMDP

print('Setting training parameters...')
# Set training settings
settings = {
    # agent settings
    'batch_size' : 32,
    'print_every' : 5000,
    'save_dir' : 'results/nets_tiger_observation',
    'iterations' : 500000,
    'eval_iterations' : 5000,
    'eval_every' : 5000,
    'save_every' : 5000,
    'initial_exploration' : 10000,
    'epsilon_decay' : 0.0001, # subtract from epsilon every step
    'eval_epsilon' : 0, # epsilon used in evaluation, 0 means no random actions
    'epsilon' : 1.0,  # Initial exploratoin rate
    'model_dims': (1,1),
    'learn_freq' : 1,

    # simulator settings
    'viz' : False,

    # replay memory settings
    'memory_size' : 100000,  # size of replay memory
    'n_frames' : 5,  # number of frames

    # learner settings
    'learning_rate' : 0.001, 
    'decay_rate' : 0.99, # decay rate for RMSprop, otherwise not used
    'discount' : 0.95, # discount rate for RL
    'clip_err' : False, # value to clip loss gradients to
    'clip_reward' : False, # value to clip reward values to
    'target_net_update' : 1000, # update the update-generating target net every fixed number of iterations
    'double_DQN' : False, # use Double DQN (based on Deep Mind paper)
    'optim_name' : 'RMSprop', # currently supports "RMSprop", "ADADELTA" and "SGD"'
    'gpu' : False,
    'reward_rescale': False,

    # general
    'seed_general' : 1723,
    'seed_simulator' : 5632,
    'seed_agent' : 9826,
    'seed_memory' : 7563

    }

print(settings)

np.random.seed(settings["seed_general"])

print('Setting up simulator...')
pomdp = TigerPOMDP( seed=settings['seed_simulator'] )
simulator = POMDPSimulator(pomdp, robs=True)

settings['model_dims'] = simulator.model_dims

print('Initializing replay memory...')
memory = ReplayMemoryHDF5(settings)

print('Setting up networks...')

class Linear(Chain):

    def __init__(self):
        super(Linear, self).__init__(
            l1=F.Bilinear(settings["n_frames"], settings["n_frames"], 200),
            l2=F.Linear(200, 100, wscale=np.sqrt(2)),
            l3=F.Linear(100, 100, wscale=np.sqrt(2)),
            l4=F.Linear(100, 50, wscale=np.sqrt(2)),
            l5=F.Linear(50, simulator.n_actions, wscale = np.sqrt(2))
        )

    def __call__(self, s, action_history):
        h1 = F.relu(self.l1(s,action_history))
        h2 = F.relu(self.l2(h1))
        h3 = F.relu(self.l3(h2))    
        h4 = F.relu(self.l4(h3))
        output = self.l5(h4)
        return output

net = Linear()

print('Initializing the learner...')
learner = Learner(settings)
learner.load_net(net)

print('Initializing the agent framework...')
agent = DQNAgent(settings)

print('Training...')
agent.train(learner, memory, simulator)

print('Loading the net...')
learner = agent.load(settings['save_dir']+'/learner_final.p')

ind_max = learner.val_rewards.index(max(learner.val_rewards))
ind_net = settings['initial_exploration'] + ind_max * settings['eval_every']
agent.load_net(learner,settings['save_dir']+'/net_%d.p' % int(ind_net))

np.random.seed(settings["seed_general"])

print('Evaluating DQN agent...')
print('(reward, MSE loss, mean Q-value, episodes - NA, time)')
reward, MSE_loss, mean_Q_value, episodes, time, paths, actions, rewards = agent.evaluate(learner, simulator, 50000)
print(reward, MSE_loss, mean_Q_value, episodes, time)
