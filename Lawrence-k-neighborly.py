from korthogonal4 import *
from Representative_matrices import get_configuration
import numpy as np


#--------------------------------------------------------------
# Select the parameters of the computation

r = 8  # rank of the matroid
n = 11  # number of elements of the matroid

k_neighborly = 3  # desired neighborliness

# The number of (k_neighborly)-neighborly reorientations
# is obtained by counting (k_neighborly + 1)-orthogonal topes.

k_orthogonal = k_neighborly + 1

#--------------------------------------------------------------
# Get the corresponding matrix and indices

A, indices = get_configuration(r, n)

# Select the corresponding input file
input_file = f"data_r{r}_n{n}.txt"

#--------------------------------------------------------------
# Create the ground set and the supports of bases and circuits

elements = list(range(n))

SB = support_basis_lex(elements, r)
SC = support_circuits(elements, r)

#--------------------------------------------------------------
# Compute the number of k_neighborly reorientations
# for each reorientation class

Counter = 0

output_file = f"results_r{r}_n{n}_k{k_neighborly}.txt"

with open(input_file, "r") as archive, open(output_file, "w") as file:

    for line in archive:

        Counter += 1

        Vector = convert_chirotope(line)

        if len(Vector) != len(indices):
            raise ValueError(
                f"Line {Counter}: expected {len(indices)} signs, "
                f"but received {len(Vector)}."
            )

        np.put(A, indices, Vector)

        Chirotope = []

        for basis in SB:
            product = []

            for i in basis:
                product.append(A[basis.index(i)][i])

            signature = np.prod(product)
            Chirotope.append(signature)

        Circuits = create_circuits(
            Chirotope, elements, r, SB, SC
        )

        X = Orthogonalidad(r, n, Circuits, k_orthogonal)

        print("Class", Counter, X)
        file.write("{}\n".format(X))

#--------------------------------------------------------------

print("\nComputation completed.")
print(f"Number of reorientation classes: {Counter}")
print(f"Results saved in: {output_file}")