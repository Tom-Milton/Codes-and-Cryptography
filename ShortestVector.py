# Run from the command line as exampled below. Press q to terminate and produce ShortestVectorOutput.txt
# python ShortestVector.py latticeBasis.txt difference
# python ShortestVector.py latticeBasis.txt average
# python ShortestVector.py latticeBasis.txt both

import re
import sys
import random
import keyboard
import numpy as np


def extract_basis(file):
    with open(file, 'r') as f:
        # Remove everything but the matrix elements
        nums = re.findall(r'\[([^][]+)\]', f.read())
        # Load matrix from list of strings
        B = np.loadtxt(nums, delimiter=',', dtype=int)
        return B


def modified_average(x, y):
    # Creates new y so that new average is always on the lattice
    new_y = [yi if ((xi + yi) % 2 == 0) else yi-1 if (random.random()<0.5) else yi+1 for xi, yi in zip(x, y)]
    return np.add(x, new_y)//2  # Using floor to return an int. It will never round due to our manipulations above


def export_results(B, u, x, norm):
    vars = [('B', B), ('u', np.atleast_2d(u).T), ('norm', norm), ('x', np.atleast_2d(x).T)]
    text = ["{} = \n{}\n\n".format(i[0], i[1]) for i in vars]

    with open('ShortestVectorOutput.txt', 'w') as f:
        f.write(''.join(text))


def sieving(B, mode):
    # Initialise random vectors
    x_array = [np.random.randint(-10, 10, size=len(B)) for i in range(100)]
    # Calculate lattice points from random vectors
    u_array = [np.matmul(B, x) for x in x_array]
    # Calculate norm of lattice points
    norm_array = [np.linalg.norm(u) for u in u_array]

    # Store (u, x, norm) for lattice point u s.t. u=Bx
    lattice_points = list(zip(u_array, x_array, norm_array))
    # Sort lattice points by norm
    lattice_points = sorted(lattice_points, key=lambda x:x[2])
    # Store shortest lattice point from sorted list
    minimum = lattice_points[0]

    while True:
        # Halt and return on keypress
        if keyboard.is_pressed('q'):
            return minimum

        # Choose random pair of lattice points
        u, v = random.sample(lattice_points, 2)

        # Calculate new u, x, and norm
        if mode == 'difference': new_x = np.subtract(u[1], v[1])
        elif mode == 'average': new_x = modified_average(u[1], v[1])
        elif mode == 'both': new_x = np.subtract(u[1], v[1]) if random.random()<0.5 else modified_average(u[1], v[1])
        new_u = np.matmul(B, new_x)
        new_norm = np.linalg.norm(new_u)

        # Discard vectors that are longer than parents
        if (new_norm > u[2]) or (new_norm > v[2]):
            continue
        # Discard collisions
        elif next((True for elem in lattice_points if np.array_equal(elem[0], new_u)), False):
            continue
        # Discard zero vector
        elif not np.any(new_u):
            continue
        
        else:
            # Otherwise add to lattice_points
            lattice_points.append((new_u, new_x, new_norm))
            # Remove longest lattice point to keep the size of the set constant          
            lattice_points = sorted(lattice_points, key=lambda x:x[2])[:-1]

            # Store point if new minimum is found
            if new_norm < minimum[2]:
                minimum = lattice_points[0]
                print(minimum[2], minimum[0], minimum[1])


def main(file, mode):
    B = extract_basis(file)
    u, x, norm = sieving(B, mode)
    export_results(B, u, x, norm)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])


################################################## OLD CODE ##################################################


def old_modified_average(B, u, x, v, y):
    # Creates new y and then v so that new average is always on the lattice
    new_y = [yi if ((xi + yi) % 2 == 0) else yi-1 if (random.random()<0.5) else yi+1 for xi, yi in zip(x, y)]
    new_v = np.matmul(B, new_y)
    return np.add(u, new_v)//2  # Using floor to return an int. It will never round due to our manipulations


def old_sieving(B, mode):
    B_inverse = np.linalg.inv(B)

    # Initialise random vectors
    x_array = [np.random.randint(-10, 10, size=len(B)) for i in range(100)]
    # Calculate lattice points from random vectors
    u_array = [np.matmul(B, x) for x in x_array]
    # Calculate norm of lattice points
    norm_array = [np.linalg.norm(u) for u in u_array]

    # Store (u, x, norm) for lattice point u s.t. u=Bx
    lattice_points = list(zip(u_array, x_array, norm_array))
    # Sort lattice points by norm
    lattice_points = sorted(lattice_points, key=lambda x:x[2])
    # Store shortest lattice point from sorted list
    minimum = lattice_points[0]

    while True:
        new_lattice_points = []
        norm_sum = 0
        counter = 0

        while len(new_lattice_points) != len(lattice_points):
            counter += 1

            # Add random lattice points if timeout reached
            if counter == 100000:
                missing_x = [np.random.randint(-10, 10, size=len(B)) for i in range(len(lattice_points) - len(new_lattice_points))]
                missing_u = [np.matmul(B, x) for x in missing_x]
                missing_norm = [np.linalg.norm(u) for u in missing_u]
                new_lattice_points = new_lattice_points + list(zip(missing_u, missing_x, missing_norm))
                break

            # Halt and return on keypress
            if keyboard.is_pressed('q'):
                return minimum

            # Choose random pair of lattice points
            u, v = random.sample(lattice_points, 2)

            # Calculate new u, x, and norm
            if mode == 'difference': new_u = np.subtract(u[0], v[0])
            elif mode == 'average': new_u = old_modified_average(B, u[0], u[1], v[0], v[1])
            new_x = np.matmul(B_inverse, new_u).round().astype(int)
            new_norm = np.linalg.norm(new_u)

            # Discard vectors that are longer than parents
            if (new_norm > u[2]) or (new_norm > v[2]):
                continue
            # Discard collisions
            elif next((True for elem in new_lattice_points if np.array_equal(elem[0], new_u)), False):
                continue
            # Discard zero vector
            elif not np.any(new_u):
                continue
            
            else:
                # Otherwise add to lattice_points
                new_lattice_points.append((new_u, new_x, new_norm))
                norm_sum += new_norm

                # Store point if new minimum is found
                if new_norm < minimum[2]:
                    minimum = new_lattice_points[-1]
                    print(minimum[2], minimum[0], minimum[1])
        
        # Throw away previous lattice points and perform sieving on new points
        lattice_points = new_lattice_points