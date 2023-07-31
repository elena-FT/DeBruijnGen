# Projet DeBruijnGen

Le projet DeBruijnGen est un ensemble de scripts Python conçus pour construire, fusionner et découper le graphe de De Bruijn à partir de données génomiques.

## Description

Le projet [Nom du Projet] propose trois scripts Python pour effectuer différentes étapes de traitement sur les données génomiques :

1. `debruijn_build.py`: Ce script prend en entrée un fichier au format FASTA contenant des séquences génomiques et construit le graphe de De Bruijn correspondant. Le graphe de De Bruijn est ensuite enregistré dans un fichier de sortie au format graph.

2. `debruijn_merge.py`: Ce script prend en entrée un fichier graphe de De Bruijn (au format graph) et fusionne les nœuds qui représentent les complémentaires inverses. Le graphe de De Bruijn résultant est enregistré dans un nouveau fichier.

3. `debruijn_cut.py`: Ce script prend en entrée un fichier graphe de De Bruijn (au format graph) et élimine les nœuds "tips" du graphe. Le graphe résultant, sans les nœuds tips, est enregistré dans un nouveau fichier.

## Utilisation

### Compilation (si applicable)

Aucune compilation n'est nécessaire pour les scripts Python. Vous pouvez directement exécuter les commandes suivantes.

### Dépendances

Les scripts nécessitent Python 3 pour s'exécuter. Aucune autre dépendance externe n'est requise.

## Installation

Aucune installation spécifique n'est requise. Les scripts Python peuvent être exécutés directement depuis le répertoire où ils sont situés.

## Instructions d'utilisation

Pour utiliser les scripts, ouvrez un terminal et exécutez les commandes suivantes :

1. Construire le graphe de De Bruijn à partir du fichier FASTA :

```python3 debruijn_build.py reads_gfp_l62_seed5_N50000_norc.fasta 25 reads_gfp_l62_seed5_N50000_norc_raw.graph```


2. Fusionner les nœuds représentant les complémentaires inverses du graphe :

```python3 debruijn_merge.py reads_gfp_l62_seed5_N50000_norc_raw.graph reads_gfp_l62_seed5_N50000_norc_merged.graph```


3. Éliminer les nœuds tips du graphe :

```python3 debruijn_cut.py reads_gfp_l62_seed5_N50000_norc_merged.graph reads_gfp_l62_seed5_N50000_norc_cut.graph```


Assurez-vous d'avoir les droits d'exécution sur les fichiers Python pour pouvoir les lancer.

## Structure des fichiers

Le projet contient les trois fichiers suivants :

- `debruijn_build.py`: Script pour construire le graphe de De Bruijn à partir d'un fichier FASTA.
- `debruijn_merge.py`: Script pour fusionner les nœuds représentant les complémentaires inverses du graphe.
- `debruijn_cut.py`: Script pour éliminer les nœuds tips du graphe.

## Contact

Si vous avez des questions ou des commentaires sur le projet, vous pouvez contacter les auteurs du projet :

- Elena Fouillet (contact@elenafouillet.com)
- Louise Hartmann (louise.hartmann@example.com)



