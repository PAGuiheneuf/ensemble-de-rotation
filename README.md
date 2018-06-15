# Ensembles de rotation

Scripts python et code C++ pour calculer les ensembles de rotation en parallèle.

## Usage

Se compile avec:
`g++ EnsRot.cpp -std=c++11 -fopenmp -Ofast -o rotationEnsemble`

Avant de lancer le code, spécifier le nombre de threads OpenMP à utiliser
`export OMP_NUM_THREADS=4`

Lancer le code avec:
`./rotationEnsemble`
