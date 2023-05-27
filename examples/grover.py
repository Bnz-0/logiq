#!/bin/env python3
import math
from logiq import ket, Op, CanonBasis, matrix, qbit

# https://en.wikipedia.org/wiki/Grover%27s_algorithm

def grover(db, index_state):
    vect = lambda x: ket(*(1 if i==x else 0 for i in range(len(db))))
    Uw = Op.build(
        {i: -vect(i) if i == index_state else vect(i) for i in range(len(db))}
    )
    Us = Op(2*db.state * ~db.state - matrix.Id(len(db)))

    G = Uw*Us

    for i in range(int(round(math.sqrt(len(db))))):
        G | db

    db.measure(CanonBasis(len(db)))
    # NB: the result should be the measured value, but to distinguish it the basis must have all eigenvalues different.
    return db

def main():
    max_iter_search = 10
    generate_db = lambda: qbit(0,1,1,0,1,0,0,1, normalize=True)
    db_basis = generate_db().basis

    print(f"The database contains those values: {generate_db()}")

    print("searching for |2>... ", end='')
    if any(db_basis[2] == grover(generate_db(), 2).state for _ in range(max_iter_search)):
        print("found!")
    else:
        print("not found, but it should be there...")

    print("searching for |3>... ", end='')
    if any(db_basis[3] == grover(generate_db(), 3).state for _ in range(max_iter_search)):
        print("found, but how since it doesn't exists?!")
    else:
        print("not found")

if __name__ == "__main__":
    main()
