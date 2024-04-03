
I Choisir la Norme ISO C++20 (/std:c++20):

--> dans les propriétés du projet

II Installation de CPLEX:

1. Installer IBM ILOG CPLEX
   
2. Modifier les propriétés du projet
   
   --> C/C++, Général, Autres répertoires Include : C:\Program Files\IBM\ILOG\CPLEX_Studio221\concert\include;C:\Program Files\IBM\ILOG\CPLEX_Studio221\cplex\include
   
   Attention : le nombre 221 est à adapter selon le nom de votre dossier
   
   --> C/C++, Préprocesseur, Définitions de préprocesseur : WIN32;_CONSOLE;IL_STD;_CRT_SECURE_NO_WARNINGS
   
   Si cela ne fonctionne pas, tester de mettre cela à la place : NDEBUG;_CONSOLE;IL_STD
   
   --> C/C++, Génération de code, Bibliothèque Runtime : DLL multithread (/MD)
   
   --> Editeur de liens, entrée, Dépendznces supplémentaires : C:\Program Files\IBM\ILOG\CPLEX_Studio221\concert\lib\x64_windows_msvc14\stat_mda\concert.lib;C:\Program Files\IBM\ILOG\CPLEX_Studio221\cplex\lib\x64_windows_msvc14\stat_mda\cplex2210.lib;C:\Program Files\IBM\ILOG\CPLEX_Studio221\cplex\lib\x64_windows_msvc14\stat_mda\ilocplex.lib
   
   Attention : le nombre 221 est à adapter selon le nom de votre dossier

III Ajouter un dossier "data" avec les fichiers sources:

--> c'est dans celui-ci que le code cherchera les données.

Présentation global du code:
...

Utilisation de code:
...
