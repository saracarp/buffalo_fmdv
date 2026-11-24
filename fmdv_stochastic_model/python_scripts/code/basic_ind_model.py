# load necessary packages
import abc
import collections
from dataclasses import dataclass

import numpy
import pandas
import scipy
import joblib

import matplotlib.pyplot as plt
import seaborn

import argparse
import json


# Set up classes for data storage

class Buffalo:
    state: str
    age: float

    def age_update(self, time_step):
        self.age += time_step

    def __post_init__(self):
        self.age = float(self.age)

class History:
    '''Record simulation history.'''

    def __init__(self):
        self.time = []
        self.state = []

    def update(self, time, herd):
        '''Update for the current `time` and `state`.'''
        self.time.append(time)

        new_states = collections.Counter(
            'carrier' if b.state in {'carrier_1','carrier_2','carrier_3'} else b.state
            for b in herd)

        self.state.append(new_states)

        #return(new_states['exposed'] + new_states['infectious'] + new_states['carrier'])

    def to_dataframe(self):
        '''Convert to something more useful.'''
        index = pandas.Index(self.time, name='time')
        df = pandas.DataFrame(self.state, index=index)
        report_states = ['maternal_immunity','susceptible','exposed','infectious','carrier','recovered']
        return df.reindex(columns=report_states).fillna(0).astype(int)


# Set up classes for events

class Event(metaclass=abc.ABCMeta):
    '''Base for a kind of event.'''

    @abc.abstractmethod
    def get_waiting_time(self, rng, herd, buff):
        '''Get a random sample of the waiting time to the next event.'''
        raise NotImplementedError

    @abc.abstractmethod
    def get_state_change(self, rng, herd, index):
        '''Get the change to `state` from this event occurring.'''
        raise NotImplementedError

# Waning event
class Waning(Event):
    '''Waning of 'maternal_immunity' to 'susceptible'.'''

    def __init__(self, parameters):
        self.mean_waiting_time = parameters['waning_mean_waiting_time']

    def get_waiting_time(self, rng, herd, buff):
        if buff.state == 'maternal_immunity':
            return rng.exponential(self.mean_waiting_time)
        else:
            return numpy.inf

    def get_state_change(self, rng, herd, index):
        herd[index].state = 'susceptible'
        return herd

# Birth event
class Birth(Event):
    '''Waning of 'maternal_immunity' to 'susceptible'.'''

    def __init__(self, parameters):
        self.mean_waiting_time = parameters['birth_mean_waiting_time']

    def get_waiting_time(self, rng, herd, buff):
        if buff.age > 0:
            return rng.exponential(self.mean_waiting_time)
        else:
            return numpy.inf

    def get_state_change(self, rng, herd, index):
        if herd[index].state in ['maternal_immunity', 'susceptible', 'exposed', 'infectious']:
            herd.append(Buffalo('susceptible', 0))
        else:
            herd.append(Buffalo('maternal_immunity', 0))
        return herd

# Death event
class Death(Event):
    '''Waning of 'maternal_immunity' to 'susceptible'.'''

    def __init__(self, parameters):
        self.mean_waiting_time = parameters['death_mean_waiting_time']

    def get_waiting_time(self, rng, herd, buff):
        return rng.exponential(self.mean_waiting_time)

    def get_state_change(self, rng, herd, index):
        del herd[index]
        return herd

# Infection event
class Infection(Event):
    '''Infection of an individual from susceptible class.'''

    def __init__(self, parameters, herd):
        self.beta_acute = parameters['acute_transmission_rate'] / collections.Counter(i.state for i in herd).total()
        self.beta_carrier = parameters['carrier_transmission_rate'] / collections.Counter(i.state for i in herd).total()

    def get_waiting_time(self, rng, herd, buff):

        if buff.state == 'susceptible':
            infectious = collections.Counter(i.state for i in herd)['infectious']
            carriers = collections.Counter(i.state for i in herd)['carrier_1'] + collections.Counter(i.state for i in herd)['carrier_2'] + collections.Counter(i.state for i in herd)['carrier_3'] 
            foi = self.beta_acute * infectious + self.beta_carrier * carriers
            try:
                self.mean_waiting_time = 1 / foi
            except ZeroDivisionError:
                return numpy.inf
            return rng.exponential(self.mean_waiting_time)
        else:
            return numpy.inf

    def get_state_change(self, rng, herd, index):
        herd[index].state = 'exposed'
        return herd

# Progression event
class Progression(Event):
    '''Progression of infection from latent to infectious period'.'''

    def __init__(self, parameters):
        self.mean_waiting_time = parameters['progression_mean_waiting_time']

    def get_waiting_time(self, rng, herd, buff):
        if buff.state == 'exposed':
            return rng.exponential(self.mean_waiting_time)
        else:
            return numpy.inf

    def get_state_change(self, rng, herd, index):
        herd[index].state = 'infectious'
        return herd

# Recovery event
class Recovery(Event):
    '''Progression of infection from infectious to carrier or recovered class'''

    def __init__(self, parameters):
        self.mean_waiting_time = parameters['recovery_mean_waiting_time']
        self.carrier_probability = parameters['carrier_probability']

    def get_waiting_time(self, rng, herd, buff):
        if buff.state == 'infectious':
            return rng.exponential(self.mean_waiting_time)
        else:
            return numpy.inf

    def get_state_change(self, rng, herd, index):
        herd[index].state = rng.choice(['carrier_1', 'recovered'],
                           p=[self.carrier_probability, 1-self.carrier_probability])
        return herd



# Carrier recover event(s)
class Carrier_Recovery_1(Event):
    '''Progression of infection from infectious to carrier or recovered class'''

    def __init__(self, parameters):
        self.mean_waiting_time = parameters['carrier_recovery_mean_waiting_time']/3

    def get_waiting_time(self, rng, herd, buff):
        if buff.state == 'carrier_1':
            return rng.exponential(self.mean_waiting_time)
        else:
            return numpy.inf

    def get_state_change(self, rng, herd, index):
        herd[index].state = 'carrier_2'
        return herd

class Carrier_Recovery_2(Event):
    '''Progression of infection from infectious to carrier or recovered class'''

    def __init__(self, parameters):
        self.mean_waiting_time = parameters['carrier_recovery_mean_waiting_time']/3

    def get_waiting_time(self, rng, herd, buff):
        if buff.state == 'carrier_2':
            return rng.exponential(self.mean_waiting_time)
        else:
            return numpy.inf

    def get_state_change(self, rng, herd, index):
        herd[index].state = 'carrier_3'
        return herd

class Carrier_Recovery_3(Event):
    '''Progression of infection from infectious to carrier or recovered class'''

    def __init__(self, parameters):
        self.mean_waiting_time = parameters['carrier_recovery_mean_waiting_time']/3

    def get_waiting_time(self, rng, herd, buff):
        if buff.state == 'carrier_3':
            return rng.exponential(self.mean_waiting_time)
        else:
            return numpy.inf

    def get_state_change(self, rng, herd, index):
        herd[index].state = 'recovered'
        return herd



# Creat simulate function

def key_tune(obj):
    return(obj.time)

def simulate(events, time_span, herd_start, seed=None):
    '''Simulate with Gillespie's next-reaction method.
    `events` is a list of all the events that can occur.
    `time_span` is `(time_start, time_end)`.
    `state_start` is the state at `time_start`.
    `seed` is for the random-number generator.'''
    # Simplify by ensuring events is not empty
    assert len(events) > 0
    (time, time_end) = time_span
    # Simplify by ensuring `time_end` is finite
    assert numpy.isfinite(time_end)
    herd = herd_start.copy()
    # Build the random-number generator
    rng = numpy.random.default_rng(seed)
    # Initialize storage of the simulation history
    history = History()
    while True:

        n_inf = history.update(time, herd)

        if time >= time_end:
            break  # Exit the `while` loop
        #if n_inf == 0:
        #   break

        herd_event_next = []
        buff_event_next = collections.namedtuple("buff_event_next", ['buff_index','time','event_type'])

        for i in range(len(herd)):
            waiting_times = [
                event.get_waiting_time(rng, herd, herd[i])
                for event in events
                ]
            t_next = min(waiting_times)
            i_next = numpy.argmin(waiting_times)
            next_e = buff_event_next(i, t_next, events[i_next])
            herd_event_next.append(next_e)

        next_event = min(herd_event_next, key = key_tune)   
        time_next = time + next_event.time
        if time_next <= time_end:
            time = time_next
            # Find the position of the minimum in `waiting_times`
            # to find the corresponding event in `events`
            event_next = next_event.event_type
            herd = event_next.get_state_change(rng, herd, next_event.buff_index)

            for i in herd:
                i.age_update(next_event.time)
        else:
            # The next event is after `time_end`
            time = time_end

    return history.to_dataframe()

def run_one(seed = None):
    rng = numpy.random.default_rng(seed)


    herd_start = []


    for i in range(n_buff-n_inf):
        buff_state = str(rng.choice(states, p = p_start))
        buff_age = rng.uniform(0,14)
        herd_start.append(Buffalo(buff_state, buff_age))

    for i in range(n_inf):
        herd_start.append(Buffalo('infectious', rng.uniform(0,14)))


  # `events` is a list of (instances) of kinds of `Event()`
    events = (
        Waning(parameters),
        Birth(parameters),
        Death(parameters),
        Infection(parameters, herd_start),
        Progression(parameters),
        Recovery(parameters),
        Carrier_Recovery_1(parameters),
        Carrier_Recovery_2(parameters),
        Carrier_Recovery_3(parameters),
        # etc
    )

    result = simulate(events, (0, 365*n_yrs), herd_start, rng)
    #result = result[['maternal_immunity','susceptible','infectious','carrier','recovered']]
    return result

# Read in parameters and set them up correctly

def load_params(path):
    with open(path) as f:
        params = json.load(f)
    # Fail early on bad probabilities
    total = sum(params["p_start"])
    if not numpy.isclose(total, 1.0):
        raise ValueError(f"state_probs sum to {total}, not 1")
    return params


def parse_args():
    parser = argparse.ArgumentParser(description="Run the buffalo disease model.")
    parser.add_argument("params", help="path to JSON parameter file")
    parser.add_argument("--seed", type=int, default=None,
                        help="override the seed in the parameter file")
    return parser.parse_args()

args = parse_args()
parameterss = load_params(args.params)

print(parameters)

