import os
import hashlib
import gmpy2
from gmpy2 import mpz, random_state
from elliptic import is_point_on_curve,point_addition, scalar_multiplication_double_and_add, n, p, a, b, d, Gx, Gy, Q


def get_random_k():

    # Create a random state object
    seed_value = 2023
    rand_state = random_state(seed_value)

    # Generate a random mpz number of 256 bits
    random_mpz = gmpy2.mpz_rrandomb(rand_state, 256)

    k = random_mpz % (n-2) + 1 # k in [1, n-1]

    return k


# ECDSA signing algorithm using gmpy2
def ecdsa_sign(message):
    # Hash the message using SHA-256
    message_bytes = message.encode('utf-8')
    message_hash = hashlib.sha256(message_bytes).digest()
    z = message_hash.hex().upper()
    z = gmpy2.mpz('0x' + z)

    

    # Generate a random number k from [1, p-1]
    k = get_random_k()
    

    # Calculate the curve point (x1, y1) = k * G
    (x1, y1) = scalar_multiplication_double_and_add(Gx, Gy, k)

    # Calculate r = x1 mod n
    r = x1 % n
    if r == 0:
        return ecdsa_sign(message)  # Start over if r is 0

    # Calculate s = k^-1(z + r*private_key_d) mod n
    s = (gmpy2.invert(k, n) * (z + r * d)) % n
    if s == 0:
        return ecdsa_sign(message) # Start over if s is 0


    # Return the signature pair (r, s)
    return (r, s)


# ECDSA Verification Algorithm
def ecdsa_verify(message, signature):
    r, s = signature

    # Step 1: Check that Q_A is not equal to the identity element O, and its coordinates are otherwise valid.
    if Q == (0, 0) or not is_point_on_curve(*Q):
        print('not on curve')
        return False

    # p * Q should be equal to O
    if scalar_multiplication_double_and_add(*Q, n) != (0, 0):
        print('n * Q != O')
        return False

    # Step 1: Verify that r and s are integers in [1, n - 1].
    if not (1 <= r < n) or not (1 <= s < n):
        print('r or s not in [1, n-1]')
        return False

    # Step 2: Calculate z
    message_bytes = message.encode('utf-8')
    message_hash = hashlib.sha256(message_bytes).digest()
    z = message_hash.hex().upper()
    z = gmpy2.mpz('0x' + z)

    # Step 3: Calculate u1 = z * s^-1 mod n and u2 = r * s^-1 mod n.
    s_inv = gmpy2.invert(s, n)
    u1 = (z * s_inv) % n
    u2 = (r * s_inv) % n

    # Step 5: Calculate the curve point (x1, y1) = u1 * G + u2 * Q_A.
    point_u1G = scalar_multiplication_double_and_add(Gx, Gy, u1)
    point_u2Q = scalar_multiplication_double_and_add(*Q, u2)
    point_x1y1 = point_addition(point_u1G[0], point_u1G[1], point_u2Q[0], point_u2Q[1])

    # If the resulting point is at infinity, the signature is invalid.
    if point_x1y1 == (0, 0):
        print('point at infinity')
        return False


    # Step 6: The signature is valid if r ≡ x1 (mod n), invalid otherwise.
    return r%n == point_x1y1[0] % n


message = "Hello world"

# Sign the message
r, s = ecdsa_sign(message)
print('r: {}'.format(r))
print('s: {}'.format(s))
v = ecdsa_verify(message, (r, s))
print('verified: {}'.format(v))