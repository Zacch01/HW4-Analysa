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

        # Check if the domain borders are roots
        if f(startAt) == 0:
            print('The root --> ' + str(startAt) + '\tIteration --> 0')
            startAt = startAt + 0.1
            continue

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

    if f(endAt) == 0:
        print('The root --> ' + str(endAt) + '\tIteration --> 0')


def bisectionMethod(f, left, right, iterationAllowed):
    """
    Finding the function f root between the domain [leftDomain To rightDomain]

    :param f: Our function
    :param left: The leftDomain domain of the function
    :param right: The rightDomain domain of the function
    :param iterationAllowed: The maximum iteration allowed
    :return: The root of the function if existed, else according failed message
    """
    # Search the root within the maximum allowed iteration
    for i in range(iterationAllowed):

        # Variable to store the middle of the function segment
        middle = left + (right - left) / 2

        # In case we found our root, Return the root and the iteration number
        if abs(f(middle)) < MAX_ERROR:
            return int(middle * 10 ** 5) / 10 ** 5, i + 1

        # In case the root is between the leftDomain To the middle, Update rightDomain to be the middle
        elif f(left) * f(middle) < 0:
            right = middle

        # In case the root is between the middle To the right, Update leftDomain segment to be the middle
        elif f(middle) * f(right) < 0:
            left = middle

    # In case we didn't find the root within the allowed amount iteration, Send fail message and shut down the program
    print('Failed To Find The Root, The Bisection Method Is Not Suitable For This Function')
    exit()


def newtonRaphsonMethod(f, g, currentX):
    """
    Finding the function f root

    :param f: Our function
    :param g: The derivative function of f
    :param currentX: The value of the middle domain range of the function
    :return: The root of the function if existed, else according failed message
    """
    # Search the root within the maximum allowed iteration
    for i in range(100):

        # Variable to store the next X
        nextX = currentX - f(currentX) / g(currentX)

        # In case we found our root, Return the root and the iteration number
        if abs(nextX - currentX) < MAX_ERROR:
            return int(nextX * 10 ** 5) / 10 ** 5, i + 1

        # Update the currentX to be the new one
        currentX = nextX

    # In case we didn't find the root within the allowed amount iteration, Send fail message and shut down the program
    print('Failed To Find The Root, The Newton Raphson Method Is Not Suitable For This Function')
    exit()


def secantMethod(f, previewX, currentX):
    """
    Finding the function f root between the domain range [previewX To currentX]

    :param f: Our function
    :param previewX: The leftDomain domain of the function
    :param currentX: The rightDomain domain of the function
    :return: The root of the function in the domain [previewX To currentX] if existed, else according failed message
    """

    # Search the root within the maximum allowed iteration
    for i in range(100):

        # Variable to store the next X
        nextX = (previewX * f(currentX) - currentX * f(previewX)) / (f(currentX) - f(previewX))

        # In case we found our root, Return the root and the iteration number
        if abs(nextX - currentX) < MAX_ERROR:
            return int(nextX * 10 ** 5) / 10 ** 5, i + 1

        # Update the previewX to be the currentX
        previewX = currentX

        # Update the currentX to be new one
        currentX = nextX

    # In case we didn't find the root within the allowed amount iteration, Send fail message and shut down the program
    print('Failed To Find The Root, The Secant Method Is Not Suitable For This Function')
    exit()


# Our Program Driver
if __name__ == "__main__":
    x = sp.symbols('x')
    function = x ** 4 + x ** 3 - 3 * x ** 2
    domainStart = -3
    domainEnd = 2

    while True:
        userChoice = input('Input the wanted Method\n1 --> Bisection\n2 --> Newton-Raphson\n3 --> Secant\nInput --> ')
        if userChoice != '1' and userChoice != '2' and userChoice != '3':
            print('Invalid Input, Try [One To Three]')

        else:
            print()
            break
    rootFinder(function, domainStart, domainEnd, userChoice)
