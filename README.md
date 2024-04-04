
Informations générales ____________________________________________________________________________________________________

Présentation global du code:
   
   Ce code a été réalisé dans le but de tester le modèle "Global Spectrum Interpretation" (GSI) sur des jeux de données simulés ou réels. Nous n'avons pas pour objectif d'en faire un outil utilisable pour le moment. Il s'agit d'un travail en cours et par conséquent, le code risque fortement d'évoluer, ainsi que le modèle GSI lui-même. Il n'y a pour l'instant aucune gestion d'erreurs mise en place.

Utilisation de code:

   Cette partie ne fait que résumer les méthodes à utiliser pour tester le modèle mais ne constitue en rien une documentation. La documentation se trouve dans le code lui-même.
   
   Tout le code tourne autour du type "Model". Il faut commencer par créer dans la fonction "Main" une variable de type "Model". A partir d'une instance de cette classe, il est possible d'effectuer les tâches suivantes :
   
   1- Charger un fichier de protéines : méthodes "Load_Proteins"
   
   2- Générer les peptides théoriques : méthodes "In_Silico_Digestion"
   
   3- Construire les spectres théoriques corresondants aux peptides théoriques : méthode "Build_Theoretical_Spectra"
   
   4- Définir les probabilités sur les arêtes protéines-peptides : méthodes "Define_Probabilities"
   
   5- Charger un fichier de spectres expérimentaux : méthodes "Load_Spectra"
   
   6- Générer des spectres simulés à partir des protéines chargées : méthodes "Simulated_Sample"
   
   7- Générer des scores sur les arêtes spectres-peptides : méthodes "Compute_Score" et "Compute_Score_SpecOMS"
   
   8- Générer des scores sur les arêtes spectres-peptides à partir d'un fichier : méthode "Load_Scores"
   
   9- Calculer une solution : méthode "Solve"
   
   10- Afficher la solution : méthode "Print_Solution"
   
   
   Attention : certaines méthodes ne peuvent être appelées que si d'autres l'ont déjà été. De plus, seule l'une des étapes 5 et 6 peut être appliquée (même chose pour 7 et 8).
   
   Il n'existe pour l'instant pas de méthode permettant d'enregistrer la solution dans un fichier mais cela devra très certainement être fait. De même qu'une méthode permettant de définir les probabilités à partir d'un fichier.

Mise en place (! sous visual studio 2022 !) ____________________________________________________________________________________________________

   I Choisir la Norme ISO C++20 (/std:c++20):
   
   --> dans les propriétés du projet
   
   II Installation de CPLEX :
   
   1. Installer IBM ILOG CPLEX
      
   2. Modifier les propriétés du projet
      
      --> C/C++, Général, Autres répertoires Include : C:\Program Files\IBM\ILOG\CPLEX_Studio221\concert\include;C:\Program Files\IBM\ILOG\CPLEX_Studio221\cplex\include
      
      Attention : le nombre 221 est à adapter selon le nom du dossier
      
      --> C/C++, Préprocesseur, Définitions de préprocesseur : WIN32;_CONSOLE;IL_STD;_CRT_SECURE_NO_WARNINGS
      
      Si cela ne fonctionne pas, tester de mettre cela à la place : NDEBUG;_CONSOLE;IL_STD
      
      --> C/C++, Génération de code, Bibliothèque Runtime : DLL multithread (/MD)
      
      --> Editeur de liens, entrée, Dépendances supplémentaires : C:\Program Files\IBM\ILOG\CPLEX_Studio221\concert\lib\x64_windows_msvc14\stat_mda\concert.lib;C:\Program Files\IBM\ILOG\CPLEX_Studio221\cplex\lib\x64_windows_msvc14\stat_mda\cplex2210.lib;C:\Program Files\IBM\ILOG\CPLEX_Studio221\cplex\lib\x64_windows_msvc14\stat_mda\ilocplex.lib
      
      Attention : le nombre 221 est à adapter selon le nom du dossier
   
   III Ajouter un dossier "data" à côté des fichiers sources :
   
   --> c'est dans celui-ci que le code devra pouvoir trouver les fichiers à charger.

