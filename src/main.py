import gmpy2
from gmpy2 import mpz

# Elliptic curve parameters for secp256k1
p = mpz('0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F')
a = mpz(0)
b = mpz(7)
d = mpz(4064)
Gx = mpz('0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798')
Gy = mpz('0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8')

# Elliptic curve point addition
def point_addition(x1, y1, x2, y2, p, a):
    if x1 == x2 and y1 == y2:
        # Point doubling
        # cath zero division error
        try:
            m = (3 * x1 * x1 + a) * gmpy2.invert(2 * y1, p)
        except ZeroDivisionError:
            print("ZeroDivisionError")
            print("x1: {}".format(x1))
            print("y1: {}".format(y1))
            input()
    else:
        # Point addition
        m = (y2 - y1) * gmpy2.invert(x2 - x1, p)
    
    m = m % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return x3, y3


def scalar_multiplication_double_and_add(x, y, scalar, p, a):
    # Ensure scalar is of type mpz
    scalar = mpz(scalar)

    # Convert the scalar to binary form
    scalar_bin = scalar.digits(2)

    # Initialize the result point as the identity element (0, 0)
    result_x, result_y = mpz(0), mpz(0)

    # Iterate through the bits of the scalar
    for i, bit in enumerate(scalar_bin):
        if i == 0 and bit == '1':
            # If the first bit is 1, set the result point to the original point
            result_x, result_y = x, y
            continue
        elif i == 0 and bit == '0':
            # this should never happen, cast an error
            raise ValueError("Invalid scalar provided: first bit cannot be 0")
        # Double the current point
        result_x, result_y = point_addition(result_x, result_y, result_x, result_y, p, a)

        # If the bit is 1, add the original point
        if bit == '1':
            result_x, result_y = point_addition(result_x, result_y, x, y, p, a)
        # If the bit is 0, continue
        else:
            continue

    return result_x, result_y


# Test the modified function with the same curve and base point
d = 5
result_x_da, result_y_da  = scalar_multiplication_double_and_add(Gx, Gy, d, p, a)

# The restuls of (x, y)
print("Double and add results:")
print("x: {}".format(result_x_da))
print("y: {}".format(result_y_da))