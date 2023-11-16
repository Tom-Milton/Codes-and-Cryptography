# Run from command line as exampled below. They produce files BasicRhoOutput.txt, FullRhoOutput.txt, and plaintext.txt respectively
# python PollardRho.py exampleInputRho.txt basic
# python PollardRho.py exampleInputRho.txt full
# python PollardRho.py ECDH.txt ecdh

import sys
import random
import math
import des


def extract_parameters(file):
    with open(file, 'r') as f:

        # Extract all lines that aren't commented or new lines, and then remove all whitespace
        lines = [line.strip().replace(' ','') for line in f if (line[0] != '#' and line != '\n')]
        if 'Input:' in lines: lines.remove('Input:')

        # Convert string variables in text file to global variables
        # Don't care about the problems associated with global variables as this code is not 'complex' and these variables never change
        for line in lines:
            name, value = line.split('=')
            try: globals()[name] = eval(value)
            except: globals()[name] = str(value)


def double(P):
    if P is None:  # 2∞=∞
        return None
    elif P == (P[0], -P[1]%p):  # P(-P)=∞
        return None

    x1, y1 = P
    x3 = ((3*(x1**2)+a) * pow(2*y1, -1, p))**2-(2*x1)
    y3 = ((3*(x1**2)+a) * pow(2*y1, -1, p))*(x1-x3)-y1

    return (x3%p, y3%p)


def add(P, Q):
    if P is None or Q is None:  # P+∞ = ∞+P = P
        return P or Q
    elif P == Q:  # P+P=2P
        return double(P)
    elif P == (Q[0], -Q[1]%p):  # P+(-P)=∞
        return None
    
    x1, y1 = P
    x2, y2 = Q
    x3 = ((y2-y1) * pow(x2-x1, -1, p))**2-x1-x2
    y3 = ((y2-y1) * pow(x2-x1, -1, p))*(x1-x3)-y1

    return (x3%p, y3%p)


def double_and_add(P, scalar):
    result = None
    addend = P

    # Iterate through bits of binary scalar from LSB to MSB
    for bit in '{0:b}'.format(scalar)[::-1]:

        # Adds powers of P if present in scalar
        if int(bit) == 1:
            result = add(result, addend)
        
        # Calculates powers of P
        addend = double(addend)

    return result  


def f(P, Q, x, c, d):
    # Partition into 3 sets using x coordinate
    # Functon f from page 3 of http://koclab.cs.ucsb.edu/teaching/ecc/project/2015Projects/Blumenfeld-Presentation.pdf
    if x is None:             # Same as S1
        return add(x, P), (c+1)%n, (d+0)%n
    elif (0 <= x[0] < p//3):          # S1
        return add(x, P), (c+1)%n, (d+0)%n
    elif p//3 <= x[0] < 2*(p//3):     # S2
        return double(x), (2*c)%n, (2*d)%n
    elif 2*(p//3) <= x[0] < p:        # S3
        return add(x, Q), (c+0)%n, (d+1)%n


def pollard_rho(P, Q):
    # Initialise random starting point x and calculate intitial x'
    # x initialised as a scalar multiple of P as described in page 3 of http://koclab.cs.ucsb.edu/teaching/ecc/project/2015Projects/Blumenfeld-Presentation.pdf
    c, d = random.randrange(0, n), 0
    x = double_and_add(P, c)  # x = cP + dQ for random c and d=0
    x_prime, c_prime, d_prime = f(P, Q, x, c, d)  # x' = f(x)

    # Looks for collision with x' moving twice as fast as x
    while (x != x_prime):
        x, c, d = f(P, Q, x, c, d)  # x ← f(x)
        x_temp, c_temp, d_temp = f(P, Q, x_prime, c_prime, d_prime)
        x_prime, c_prime, d_prime = f(P, Q, x_temp, c_temp, d_temp)  # x' ← f(f(x'))
    
    # Search for another collision if c=c' and d=d'
    if (c==c_prime) and (d==d_prime):
        c, c_prime, d, d_prime = pollard_rho(P, Q)

    return c, d, c_prime, d_prime


def full_pollard_rho(P, Q, c, d, c_prime, d_prime):
    gcd = math.gcd(d_prime-d, n)
    # // used to return int instead of float. All numerators always a multiple of gcd so // never rounds
    initial_l = (((c-c_prime)//gcd) * pow((d_prime-d)//gcd, -1, n//gcd)) % (n//gcd)
    # Calculate all possible values of l (there are gcd many possible solutions for l)
    possible_l = [(initial_l + i*(n//gcd))%n for i in range(gcd)]

    # for all possible l, if Q = lP we have found our discrete log
    for i in possible_l: 
        if Q == double_and_add(P, i): l=i

    return l


def decryption(xK, ciphertext):
    # Convert key to binary and fill leading zeros to get 56 bits
    xK = "{0:b}".format(xK).zfill(56)
    # Add an (unused) 1 bit after every 7 bits to get 64 bits
    xK = ''.join(v + '1'*(c % 7 == 6) for c, v in enumerate(xK))
    
    # Convert key and cipertext to bytes
    xK_bytes = int(xK, 2).to_bytes(len(xK) // 8, byteorder='big')
    ciphertext_bytes = bytes.fromhex(ciphertext)

    # Decryption of DES encryption using ECB   
    key = des.DesKey(xK_bytes)
    plaintext = key.decrypt(ciphertext_bytes, padding=True).decode('utf-8')

    # Write plaintext to text file
    with open('plaintext.txt', 'w', encoding='utf-8') as f:
        f.write(plaintext)    


def export_results(c, d, c_prime, d_prime, l=None):
    file_name = 'FullRhoOutput.txt' if l else 'BasicRhoOutput.txt'

    vars = [('p', p), ('a', a), ('b', b), ('P', P), ('n', n), ('Q', Q), ('c', c), ('d', d), ("c'", c_prime), ("d'", d_prime)]
    input_text = ["\n{} = {}".format(i[0], i[1]) for i in (vars[:6])]
    output_text = ["\n{} = {}".format(i[0],str(i[1])) for i in vars[6:]]

    with open(file_name, 'w') as f:
        f.write('Input:' + ''.join(input_text))
        f.write('\n\nCollision:' + ''.join(output_text))
        if l: f.write('\n\nDiscrete logarithm: \nl = ' + str(l))


def main(file, mode):
    extract_parameters(file)

    if mode == 'basic':
        c, d, c_prime, d_prime = pollard_rho(P, Q)
        export_results(c, d, c_prime, d_prime)
    
    if mode == 'full':
        c, d, c_prime, d_prime = pollard_rho(P, Q)
        l = full_pollard_rho(P, Q, c, d, c_prime, d_prime)
        export_results(c, d, c_prime, d_prime, l)
    
    if mode == 'ecdh':
        c, d, c_prime, d_prime = pollard_rho(P, QA)
        dA = full_pollard_rho(P, QA, c, d, c_prime, d_prime)
        xK = double_and_add(QB, dA)[0]
        decryption(xK, ciphertext)

        # dA = 1682779984167835
        # dB = 428971283427559
        # xK = (6714934996831608, 12073846457401645)[0]


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])