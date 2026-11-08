# Pêche BZH

Appli web de pêche pour la côte sud de la Bretagne, de Groix à La Baule : marées, météo croisée entre plusieurs modèles, indice de pêche, carnet de sorties et statistiques.

## Marées
Les prédictions viennent d'une analyse harmonique (UTide) faite sur 1 à 2 ans d'observations des marégraphes REFMAR/SHOM (licence ouverte Etalab) : Lorient, Port-Tudy, Le Crouesty, Larmor-Baden, Vannes, Le Croisic. Le coefficient est calculé à partir de Brest. Les données sont dans `tides.js` et couvrent jusqu'au 31/12/2027. Elles ne servent pas à la navigation : la référence reste le SHOM.

Pour régénérer les données (avant 2028) :
```
pip install utide pandas numpy
python outils/dl.py obs          # télécharge les observations REFMAR (adapter les dates dans le script)
python outils/analyse.py obs res # analyse + validation, écrit res/rapport.json
python outils/build.py res brest_apim.json tides.js  # brest_apim.json : coefficients de Brest de référence (api-maree.fr, clé démo)
```

## Données en ligne
Météo, houle et température de l'eau : Open-Meteo (CC BY 4.0), modèles AROME/ARPEGE, ECMWF, ICON et GFS. Fond de carte : © OpenStreetMap.

Le carnet reste stocké dans le navigateur : exporte-le régulièrement (JSON/CSV).
