import random
from collections import Counter # Only used for testing

class InvalidProbabilitiesError(ValueError):
    pass

class DiscreteDist(object):
    def __init__(self, values, probs, are_cumulative=False):
        if not are_cumulative and abs(1-sum(probs)) > 2e-3:
            raise InvalidProbabilitiesError("probs do not sum to 1")
        self.values = values
        self.probs = probs
        self.are_cumulative=are_cumulative

    def sample1(self):
        # Sample a single value
        r = random.random()
        cumulative_prob = 0
        for value, prob in zip(self.values, self.probs):
            if self.are_cumulative:
              cumulative_prob = prob
            else:
              cumulative_prob += prob
            if r < cumulative_prob:
                return value
        
        # Just in case the probabilities don't add up
        # to 1 due to floating point arithmetic:
        return self.values[-1]
            
    def sample(self, n=1):
        # Sample n values
        return (self.sample1() for i in xrange(n))

if __name__=="__main__":
    # Check that it works...
    dist = DiscreteDist(values=[1.1, 2.2, 3.3], probs=[0.2, 0.5, 0.3])
    print [x for x in dist.sample(10)]
    print dist.sample1()
    print Counter(dist.sample(10000))
