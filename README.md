# Ensembles de rotation

Scripts python et code C++ pour calculer les ensembles de rotation en parallèle.

## Usage

Se compile avec (au choix):
```
g++ EnsRot.cpp -std=c++11 -fopenmp -Ofast -o rotationEnsemble
g++ EnsRot_paral_angle.cpp -std=c++11 -fopenmp -Ofast -o rotationEnsemble
g++ EnsRot_paral_angle_optim.cpp -std=c++11 -fopenmp -Ofast -o rotationEnsemble
```

Avant de lancer le code, spécifier le nombre de threads OpenMP à utiliser
`export OMP_NUM_THREADS=4`

Lancer le code avec:
`./rotationEnsemble N d Tps M nParams`

## Note

La première version (EnsRot.cpp) correspond à une parallélisation en fonction de la valeur du paramètre, à utiliser pour nParams > 1 donc.

Les deuxième et troisième versions correspondent au code parallélisé en angle, à utiliser avec des valeurs de M assez grandes.
Sur ces dernières, la parallélisation n'est pas optimale, il faudrait peut être tester:
* De changer le mode de parallélisation openmp (static -> guided)
* De changer la structure des étiquettes dans les objets Edge2, l'accès se faisant en parallèle, les caches sont surement mal optimisés
* Au lieu de ruser en donnant aux étiquettes une taille égale au nombre d'angles, mieux vaudrait probablement écrire un constructeur de copie de la classe Graphe, afin de pouvoir le passer à chaque thread.