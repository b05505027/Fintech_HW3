import gmpy2
from gmpy2 import mpz

# Elliptic curve parameters for y2=x3+2x+3
# p = mpz(7)
# n = mpz(6)
# a = mpz(2)
# b = mpz(3)
# d = mpz(3)
# Gx = mpz(2)
# Gy = mpz(1)


#Elliptic curve parameters for secp256k1
p = mpz('0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F')
n = mpz('0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141')
a = mpz(0)
b = mpz(7)
d = mpz(4064)
Gx = mpz('0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798')
Gy = mpz('0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8')



def is_point_on_curve(x, y):
    # Check if the point is on the curve
    return (y * y - x * x * x - a * x - b) % p == 0
# Elliptic curve point addition
def point_addition(x1, y1, x2, y2):

    if (x1, y1) == (mpz(0), mpz(0)):
        return (x2, y2)

    elif (x2, y2) == (mpz(0), mpz(0)):
        return (x1, y1)
    
    if x1 == x2 and y1 == y2:
        if y1 == 0:
            return (mpz(0), mpz(0))
        try:
            m = (3 * x1 * x1 + a) * gmpy2.invert(2 * y1, p)
        except ZeroDivisionError:
            print("ZeroDivisionError")
            print("x1: {}".format(x1))
            print("y1: {}".format(y1))
            print("x2: {}".format(x2))
            print("y2: {}".format(y2))
            input()
    elif x1 == x2:
        return (mpz(0), mpz(0))
    else:
        # point addition
        m = (y2 - y1) * gmpy2.invert(x2 - x1, p)
    
    m = m % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return x3, y3


def scalar_multiplication_double_and_add(x, y, scalar):
    # Ensure scalar is of type mpz
    scalar = mpz(scalar)
    if scalar == 0:
        return (mpz(0), mpz(0))
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
        # Double the current point
        result_x, result_y = point_addition(result_x, result_y, result_x, result_y)
        # If the bit is 1, add the original point
        if bit == '1':
            result_x, result_y = point_addition(result_x, result_y, x, y)
        else:
            continue

    return result_x, result_y

Q = scalar_multiplication_double_and_add(Gx, Gy, d)

if __name__ == '__main__':

    # Test the modified function with the same curve and base point
    
    result_x_da, result_y_da  = scalar_multiplication_double_and_add(Gx, Gy, d)
    print("Double and add results:")
    print("x: {}".format(result_x_da))
    print("y: {}".format(result_y_da))
