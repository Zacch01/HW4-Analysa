# Author : Nitzan Tomer & Zaccharie Attias & Sharon Angdo

import math

import sympy as sp
from sympy.utilities.lambdify import lambdify


# Global Variable To Set The Wanted Accuracy Solution
MAX_ERROR = 0.00001


def rootFinder(f, startAt, endAt, selectedMethod):
    """
    Method for getting the functions Roots

    :param f: Our function
    :param startAt: The leftDomain domain of the function
    :param endAt: The rightDomain domain of the function
    :param selectedMethod: The wanted method for finding the root
    """
    # Variables to store our derivative function
    g = f.diff(x)
    h = g.diff(x)

    # Making our function be able to get an X
    f = lambdify(x, f)
    g = lambdify(x, g)
    h = lambdify(x, h)

    # Max iteration for Bisection Method
    maxIteration = int(-(math.log(MAX_ERROR / (endAt - startAt)) / math.log(2))) + 1

    # Divide our function domain range into multiply domains with 0.1 range, then search for each one of them for a root
    while startAt < endAt:

        # In case the function change its sign (Means there's at least one root)
        if f(startAt) * f(startAt + 0.1) < 0:
            if selectedMethod == '1':
                root, iteration = bisectionMethod(f, startAt, startAt + 0.1, maxIteration)

            elif selectedMethod == '2':
                root, iteration = newtonRaphsonMethod(f, g, startAt + 0.05)

            elif selectedMethod == '3':
                root, iteration = secantMethod(f, startAt, startAt + 0.1)

            print('The root --> ' + str(root) + '\tIteration --> ' + str(iteration))

        # In case the derivative function change its sign (Mean there's a possibility for a root)
        if g(startAt) * g(startAt + 0.1) < 0:

            # Getting a possible root (The return of the derivative function can be Root or an Extreme point)
            if selectedMethod == '1':
                possibleRoot, iteration = bisectionMethod(g, startAt, startAt + 0.1, maxIteration)

            elif selectedMethod == '2':
                possibleRoot, iteration = newtonRaphsonMethod(g, h, startAt + 0.05)

            elif selectedMethod == '3':
                possibleRoot, iteration = secantMethod(g, startAt, startAt + 0.1)

            # Checking the possible root is indeed a root
            if abs(f(possibleRoot)) < MAX_ERROR:
                print('The root --> ' + str(possibleRoot) + '\tIteration --> ' + str(iteration))

        # Update our domain for this iteration
        startAt = startAt + 0.1