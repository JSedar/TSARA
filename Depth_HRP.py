import os
import sys
from distutils.dir_util import copy_tree
Depth_HRP = sys.argv[0]
arguments = sys.argv[1:]
def recupRepertoires(rep):
    listeRepertoires = [] #Création de la liste qui va contenir les répertoires
    for (repertoires, sousRepertoires, fichiers) in os.walk("."): #boucle qui parcourt tous les dossier et sous dossier
        listeRepertoires.append(repertoires) # Ajout de chaque répertoire trouvé dans la liste des repertoires !!!
    return listeRepertoires # Renvoie de tous les repertoire contenue dans le chemin d'access rensigné !!!

def CreationDossierAmp_cov():
    if os.path.exists("./Amp_cov") == False: #Teste si le dossier Amp_cov Existe // si NON il le Crée si non rien !!!
        os.makedirs("./Amp_cov")     #Création du dossier Amp_cov
        copy_tree("/home/samara/NeST/Tsara/depth", "./Amp_cov") 
    
def Principale():
    valeur = '"'
    while '"' in valeur or "'" in valeur:
        valeur = str(sys.argv[1])
#= str(input("Veuillez entrer le chemin d'acces complet du dossier a traiter (sans les \" \" ou les \') : ")) # récupération du chemin d'acces rensigné par l'utilisateur
    os.chdir(valeur) #Changement du répertoire courant d'execution du programme vers le repertoire du chemin d'accés renseigné
    liste = recupRepertoires(valeur) #Recupération de tous les repertoires(execution de la fonction récupRepertoires ci-dessus) du dossier du chemin d'access rensigné par l'utilisateur
    os.system("echo -e 'sample\tHRP2\tHRP3' > ./Detection_HRP2-3.tsv")
    CreationDossierAmp_cov() #Execution de la fonction CreationDossierAmp_cov ci-dessus
    for i in liste: #Boucle qui parcourt la liste des répertoires
        if os.path.exists(i+"/output_FM_SR_DD_RG.bam") == True: #Teste si dans chaque repertoires, si il existe le fichier recherché  
            os.system("samtools depth " +  i+"/output_FM_SR_DD_RG.bam > ./Amp_cov/" + i.split("/")[1] + "depth ") #Execution de la fonction samtool sur les fichiers qui existe dans les répertoires et renvoie du résultat dans un fichier renommée selon le chemin d'access au fichier dans le dossier Amp_cov 
            os.system("Rscript ./Amp_cov/draw_coverage_plot.R ./Amp_cov/" + i.split("/")[1] + "depth ./Amp_cov/Genes_for_samaraPipeline_2strand.bed ./Amp_cov/ampliconBoudaries.tsv ")
            os.system("echo -n " + i.split("/")[1] + "$'\t' >> ./Detection_HRP2-3.tsv")
            os.system("grep '\%HRP2' ./Amp_cov/" + i.split("/")[1] + "depth | wc -l | awk '{print $1==0?0:1}' | tr '\n' '\t' >> ./Detection_HRP2-3.tsv")
            os.system("grep '\%HRPIII' ./Amp_cov/" + i.split("/")[1] + "depth | wc -l | awk '{print $1==0?0:1}' >> ./Detection_HRP2-3.tsv")
            os.system("rm ./Amp_cov/" + i.split("/")[1] + "depth")
Principale() # Execution de la fonction Principale
os.system("rm ./Amp_cov/draw_coverage_plot.R")
os.system("rm ./Amp_cov/ampliconBoudaries.tsv")
os.system("rm ./Amp_cov/Genes_for_samaraPipeline_2strand.bed")
os.system("rm ./Amp_cov/plotGenes_AV_takenFromSushi.R")

