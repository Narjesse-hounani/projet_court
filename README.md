# Projet : Détection des segments transmembranaires d'une protéine

Ce projet utilise **Biopython** et l’outil externe **mkdssp** pour calculer la surface accessible (ASA) des résidus d’une protéine à partir de son fichier PDB.
---
##  Installation
```bash
# Cloner/copier le projet et entrer dans le dossier
cd projetcourt

# Créer un environnement virtuel
python3 -m venv rendu_env

# Activer l’environnement
source rendu_env/bin/activate    # sous Linux/macOS
rendu_env\Scripts\activate       # sous Windows PowerShell

# Installer les dépendances
pip install -r requirements.txt
'''
---
Le programme mkdssp est requis par Biopython pour calculer la surface accessible.
## Installation
sudo apt update
sudo apt install dssp

#Utilisation du script
Une fois l’environnement activé et mkdssp installé, lancez le script. Il sera demandé le chemin d'un fichier pdb
