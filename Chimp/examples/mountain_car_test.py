"""
This is a place holder for real unit testing.
Right now we just overfit a simple control problem:
    - the agent tries to get to the top right corner (1,1) of a 2D map
    - action 0 takes it towards (0,0), action 1 takes it toward (1,1)
    - action 1 is optimal for all states
"""

from chimp.learners.chainer_backend import ChainerBackend
from chimp.learners.dqn_learner import DQNLearner
from chimp.utils.policies import DQNPolicy

from chimp.simulators.mdp.mountain_car import MountainCar
from chimp.simulators.mdp.mdp_simulator import MDPSimulator

import numpy as np
import pickle
import pylab as p

import chainer
import chainer.functions as F
import chainer.links as L
from chainer import Chain

settings = {

    # agent settings
    'batch_size' : 32,
    'print_every' : 1000,
    'save_dir' : 'results',
    'iterations' : 2000000,
    'eval_iterations' : 100,
    'eval_every' : 1000,
    'save_every' : 20000,
    'initial_exploration' : 50000,
    'epsilon_decay' : 0.000001, # subtract from epsilon every step
    'eval_epsilon' : 0, # epsilon used in evaluation, 0 means no random actions
    'epsilon' : 1.0,  # Initial exploratoin rate
    'learn_freq' : 1,
    'history_sizes' : (1, 0, 0), # sizes of histories to use as nn inputs (o, a, r)
    'model_dims' : (1,2),

    # simulator settings
    'viz' : False,

    # replay memory settings
    'memory_size' : 100000,  # size of replay memory
    'n_frames' : 1,  # number of frames

    # learner settings
    'learning_rate' : 0.00001,
    'decay_rate' : 0.99, # decay rate for RMSprop, otherwise not used
    'discount' : 0.95, # discount rate for RL
    'clip_err' : False, # value to clip loss gradients to
    'clip_reward' : False, # value to clip reward values to
    'target_net_update' : 2000, # update the update-generating target net every fixed number of iterations
    'double_DQN' : False, # use Double DQN (based on Deep Mind paper)
    'optim_name' : 'ADAM', # currently supports "RMSprop", "ADADELTA", "ADAM" and "SGD"'
    'gpu' : False,
    'reward_rescale': False,

    # general
    'seed_general' : 1723,
    'seed_simulator' : 5632,
    'seed_agent' : 9826,
    'seed_memory' : 7563

    }

class TestNet(Chain):

    def __init__(self):
        super(TestNet, self).__init__(
            l1=F.Linear(settings['model_dims'][1], 20, bias=0.0),
            l2=F.Linear(20, 10, bias=0.0),
            bn1=L.BatchNormalization(10),
            l3=F.Linear(10, 10),
            l4=F.Linear(10, 10),
            bn2=L.BatchNormalization(10),
            lout=F.Linear(10, simulator.n_actions)
        )
        self.train = True
        # initialize avg_var to prevent divide by zero
        self.bn1.avg_var.fill(0.1),
        self.bn2.avg_var.fill(0.1),

    def __call__(self, ohist, ahist):
        h = F.relu(self.l1(ohist))
        h = F.relu(self.l2(h))
        h = self.bn1(h, test=not self.train)
        h = F.relu(self.l3(h))
        h = F.relu(self.l4(h))
        h = self.bn2(h, test=not self.train)
        output = self.lout(h)
        return output


def car_sim(nsteps, simulator, policy, verbose=False):
    mdp = simulator.model

    # re-initialize the model
    simulator.reset_episode()

    rtot = 0.0
    xpos = np.zeros(nsteps)
    vel = np.zeros(nsteps)
    # run the simulation
    input_state = np.zeros((1,2), dtype=np.float32)
    for i in xrange(nsteps):
        state = simulator.get_screenshot()
        input_state[0] = state
        a = policy.action((input_state,None))
        simulator.act(a)
        r = simulator.reward()
        rtot += r
        xpos[i], vel[i] = state
        if simulator.episode_over():
            break
    return rtot, xpos, vel


mdp = MountainCar()
simulator = MDPSimulator(mdp)

net = pickle.load(open("../chimp/pre_trained_nets/mountain_car.net", "rb")) 
backend = ChainerBackend(settings)
backend.set_net(net)
learner = DQNLearner(settings, backend)

policy = DQNPolicy(learner)

r, xtrace, vtrace = car_sim(300, simulator, policy, verbose=True)

p.plot(xtrace); p.plot(10.0*vtrace)
p.show()
