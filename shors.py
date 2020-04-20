import numpy as np
from Qubit import *
from math import ceil, log, sqrt
from collections import Counter
from statistics import mode
from fractions import Fraction

"""Simulation of Shor's Algorithm based off of the "Register Class" written by Chandler Jones.
Some code below is not used in the actual "implementation" but still written to remain faithful
to the actual steps of the algorithm

Author: Harrison Costantino"""



#Main loop; N is number to factor
def main(N, attempts=None):

	#quick check of small primes
	#turned off to demonstrate quantum pieces of alg
	if False:
		if N % 2 == 0:
			return 2
		elif N % 3 == 0:
			return 3
		elif N % 5 == 0:
			return 5

	if not attempts:
		attempts = 2*int(log(N, 2))


	guesses = []
	tries = 0
	while len(guesses) == 0 and tries < 2:
		for i in range(attempts):
			print("\n•••••••••••••••••\nIteration {} of {}\n•••••••••••••••••\n".format(i+1, attempts))
			guess = shors_alg(N)
			if guess:
				guesses.append(guess)
		if len(guesses) == 0 and tries == 0:
			print("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nOnly Encountered Degenerate Cases--Trying Again\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
		tries += 1


	try:
		factor = mode(guesses)
	except:
		if len(guesses) > 0:
			print("\nNo Potential Factor Most Likely")
			print("Potential factors are: {}".format(guesses))
			factor = None

	if tries == 2:
		print("\n\nNot able to find any factors after 2 attempts\n\nI'm so sorry I failed you")
		return None
	return factor


#Credit to https://www.geeksforgeeks.org/python-program-for-basic-and-extended-euclidean-algorithms-2/
def euclid_alg(a, b, x=1, y=1):
	if a == 0 :
		x = 0
		y = 1
		return b
	x1 = 1
	y1 = 1
	gcd = euclid_alg(b%a, a, x1, y1)
	x = y1 - (b/a) * x1
	y = x1
	return gcd


#Computes one iteration of Shor's Algorithm
def shors_alg(N):

	a = np.random.randint(1, N)
	print("The random number is {}".format(a))

	#Again turned off to focus on quantum aspects of alg; checks to see if factor randomly chosen
	if False:
		if euclid_alg(a, N) != 1:
			return euclid_alg(a, N)

	K, n = ceil(log(N**2, 2)), ceil(log(N, 2))
	Q = 2**K
	source = Register(K)
	source = source.walsh() #can also use qft; walsh slightly faster & doesnt introduce complex
	target = Register(K)

	#Quantum oracle U_f where f(x) = a^x mod N
	vals = []
	for q in range(Q):
		b = (a**q) % N
		vals.append(b)

	#Storing result of oracle in second register
	tally = Counter(vals)
	amps = []
	for i in range(Q):
		amps.append(sqrt(tally.get(i, 0)/Q))
	target = Register(amplitudes=amps)

	#Choosing an order r and setting second register to align with measurement
	r = target.measure()
	amps = 2**n * [0]
	amps[r] = 1
	target = Register(amplitudes=amps)
	#Second register no longer relevant

	#Collapse first register to be consistent with measurement
	total = 0
	amps = Q*[0]
	for q in range(Q):
		b = (a**q) % N
		if b == r:
			total += 1
			amps[q] = 1
	amps = [sqrt(1/total) * x for x in amps]
	source = Register(amplitudes=amps)

	#Apply QFT to first register to bring out the period
	source = source.QFT()

	#Measurement used in finding the period
	C = source.measure()
	print("Measured {} from the first register and {} from the second".format(r, C))

	#This is "cheating" but cont_fraction_expansion fails for C = 0
	if C == 0:
		print("Retrying this iteration (need non-zero value)\n")
		return shors_alg(N)

	#Determining the period; based off of Box 8.1 in the text; classical
	r = cont_fraction_expansion(C, Q, N)
	print("Found the period to be {}".format(r))

	#check to make sure we're in the good case: r even and a**(r/2) = -1 mod N
	if not ((r % 2 != 0) or ((a**(r//2) % N) == N-1)):
		print("\nDegenerate Case\n")
		return None

	#Find a potential factor
	y = a**(r//2)
	guess1 = euclid_alg(N, y-1)
	guess2 = euclid_alg(N, y+1)

	if guess1 != 1 and guess1 != N:
		print("\nFound potential factor {}\n".format(guess1))
		return guess1
	elif guess2 != 1 and guess2 != N:
		print("\nFound potential factor {}\n".format(guess2))
		return guess2
	print("\nDid not yield non-trivial factor")
	return None

#Computes the period; see Box 8.1/pg 167 in the text
def cont_fraction_expansion(C, Q, N):

	#initalizing the first values
	a_0 = int(C/Q)
	eps_0 = C/Q - a_0
	p_0 = a_0
	q_0 = 1

	a = int(1/(eps_0))
	eps = 1/(eps_0) - a

	#I think if this happens we're in the easy case where we just need to reduce C/Q
	if eps == 0:
		r = Fraction(C, Q).denominator
		return r

	p_1 = a*a_0 + 1
	q_1 = a

	if q_0 < N and N <= q_1:
		return q_0

	a = int(1/(eps))
	eps = 1/eps-a
	p = a*p_1 + p_0
	q = a*q_1 + q_0
	prev_q = q_1
	prev_p = p_1

	#iteration starts with current q_i = q_2; after first loop have q_3
	while not (prev_q < N and N <= q):
		a = int(1/(eps))
		eps = 1/eps-a
		temp_p = p
		temp_q = q
		p = a*p + prev_p
		q = a*q + prev_q
		prev_p = temp_p
		prev_q = temp_q
	return prev_q



#Takes in user input and runs it
N = input("Which number would you like to factor?\nHeavily recommend using numbers less than 33, less than 23 is best.\n")
N = int(N)
attempts = input("\n\nHow many iterations of Shor's Algorithm would you like to run?\nMore gives higher probability of success; type 'd' for the default.\n")
if attempts == "d":
	attempts = None
else:
	attempts = int(attempts)

factor = main(N, attempts)
if factor:
	print("\nWe found the factor {} with corresponding factor {}. There may be others".format(factor, N//factor))
