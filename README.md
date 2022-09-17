# Wind-engineering

L'objectif de ce cours était d'étudier l'influence du vent sur différentes structures : ponts, lignes haute tension, véhicules (trains, voitures) et bâtiments de grande hauteur. Le cours a été dispensé à l'Université Politechnique de Milan (Politecnico di Milano). Les commentaires sont ainsi en italien et le rapport de projet en anglais. 

Le projet choisi consite à étudier l'influence du vent sur le tablier ("deck") d'un pont. Deux parties composent le projet. Dans un premier temps, l'étude est concentrée sur une seule section du tablier. Des analyses statique puis dynamique permettent de déterminer les déplacements de la section ainsi que ses instabilités (buffeting & flutter instabilities). La seconde partie concerne l'ensemble du tablier divisé en sections. A l'aide d'une analyse modale, les déplacements du tablier entier sont déterminés ainsi que les modes présentant des instabilités. 
L'étude est centrée sur un nombre réduit de modes (14) puisque pour une structure de cette taille, les modes les plus problématiques concernent les basses fréquences. 

Le code contenant l'analyse de la section s'intitule : [Project_section](https://github.com/PalomaPerrin/Wind-engineering/blob/main/Project_section.m)

Le code contenant l'analyse de la du tablier entier s'intitule :  [Project_fullbridge](https://github.com/PalomaPerrin/Wind-engineering/blob/main/Project_fullbridge.m)

Les données et fonctions permettant de faire fonctionner le code sont disponibles au téléchargement dans le même répertoire.
