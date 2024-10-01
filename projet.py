import numpy as np
import pandas as pd
import os
from os import sys
from Bio import SeqIO
from Bio.Blast import NCBIXML
import platform #module qui permet d'accéder aux données du système d'exploitation et la version de l'interpréteur.
import subprocess

def get_data_path():
    """
    fonction qui permet d uniformiser les chemin(path) entre les differents ordinateurs
    """
    return os.path.join(os.getcwd(),"data")

def read_table_tsv(assembly_accession):
    """
    Verifie que le dossier genome existe,et que tsv est un ficher et le lire 
    si il est present.
    param: assembly_accession est l identifiant unique d'un genome donne
    :return: un dataframe pandas contenant les informations sur le genome
    """
    try:
        list_tsv=[x for x in os.listdir(os.path.join(get_data_path(),"genomes",assembly_accession))]
        if list_tsv:
            for tsv in list_tsv:
                if os.path.isfile(os.path.join(get_data_path(),"genomes",assembly_accession,f"annotation_{assembly_accession}.tsv")):
                    genome_df=pd.read_csv(os.path.join(get_data_path(),"genomes",
                                                       assembly_accession,f"annotation_{assembly_accession}.tsv" ),sep="\t")
                else:
                    print(f"{assembly_accession} est introuvable dans :{list_directory}")
            return genome_df
    except FileNotFoundError :
        print(f"chemin du fichier {assembly_accession}:est introuvable")
        
        
def get_protein(protein_id,assembly_accession):
    """
    Fonction qui renvoie ["Protein_Id","Gene_Id", "Gene_Name"] a partir d'un seul de ces 3 critères
    :param protein_id:Protein_Id soit Gene_Id, soit Gene_Name de la protein d'interet
    :param assembly_accession:est l identifiant unique d un genome_query
    :return index_protein
    """
    genome_df=read_table_tsv(assembly_accession)
    index_protein=[]
    index =genome_df.loc[(genome_df["Protein_Id"]==protein_id)|
                         (genome_df["Gene_Id"]==protein_id)|
                         (genome_df["Gene_Name"]==protein_id)].index[0]
    for protein_id_name in ['Protein_Id','Gene_Id', "Gene_Name"]:
        index_protein.append(genome_df.loc[genome_df.index == index][protein_id_name].iloc[0])
    return index_protein


def read_seq_protein(protein_id,assembly_accession):
    """
    Fonction qui verifie que protein.faa est un fichier et le parse et renvoi un seqricord(mais n'a pas vraiment d'utilite autre)
    :param protein_id:id de la protein d interet
    :param assembly_accession:est l’identifiant unique d'un genome_query
    :return seq_records
    """
    try: 
        protein_faa_file = os.path.join(get_data_path(), "genomes",assembly_accession,"protein.faa")
        if os.path.isfile(protein_faa_file):
            seq_records = [seq_record for seq_record in SeqIO.parse(protein_faa_file, "fasta") if seq_record.id == protein_id]
            return seq_records
        else:
            return "Le fichier protein.faa n est pas valide pour cette souche"
    except FileNotFoundError:
        return f"Chemin du fichier {protein_faa_file} introuvable"

def run_blastp(protein_id,assembly_accession,subject_assembly_accession,my_os_type):
    """
    Retourne le chemin blastp en fonction du type de systeme d'exploitation.
    ensuite verifie que le fichier protein.faa est un fichier et le parse et 
    enregistre le seqricord sous format fasta dans un nouveau fichier et execute le blast
    :param protein_id:id de la protein d'interet
    :param assembly_accession:est l identifiant unique d un genome_query
    :param subject_assembly_accession: est l’identifiant unique d'un genome_subject
    :param my_os_type:le systeme d exploitation
    :return: tableau de syntenie_df
    """
    if  my_os_type == "linux":
        blastp_path=os.path.join(get_data_path(),"Blast+/linux/bin/blastp")
    elif  my_os_type == "macosx":
        blastp_path=os.path.join(get_data_path(),"Blast+/macosx/bin/blastp")
    elif  my_os_type == "win":
        blastp_path="blastp"
    else:
        return f"os choisi{my_os_type} n'est pas valide"
    protein_faa_file = os.path.join(get_data_path(), "genomes",assembly_accession,"protein.faa")
    if os.path.isfile(protein_faa_file):
            seq_records = [seq_record for seq_record in SeqIO.parse(protein_faa_file, "fasta") if seq_record.id == protein_id]
            if seq_records:
                SeqIO.write(seq_records, os.path.join("data",protein_id+".fasta"), "fasta")
    query_file = os.path.join(get_data_path(),protein_id+".fasta")
    subject_db = os.path.join(get_data_path(), "genomes",subject_assembly_accession,subject_assembly_accession)
    evalue = "0.001"
    out_file ="result_blast.xml"
    outfmt = "5"
    cmd_line = [blastp_path, "-query", query_file, "-db", subject_db, "-evalue", evalue, "-out", out_file, "-outfmt", outfmt]
    subprocess.run(cmd_line, shell=True)
    return f"{out_file}"

def save_to_blast(protein_id,assembly_accession):
    """
    fonction qui enregitre le resultat d'un blast sous format xml
    :param protein_id:id de la proteine et sert aussi a nommer le fichier.xml
    :param assembly_accession:est l identifiant unique du génome et sert aussi a nommer le fichier.xml
    :return le chemin du fichier resultat du blast enregistre
    """
    fh= open("result_blast.xml","r")
    with open(os.path.join(get_data_path(),f"{protein_id}_{assembly_accession}.xml"),"w") as save_to:
        save_to.write(fh.read())
        save_to.close()
    fh.close()
    return os.path.join(get_data_path(),f"{protein_id}_{assembly_accession}.xml")

def info_best_hit(protein_id,assembly_accession):
    """
     Extrait les informations sur le meuilleur hit proteine a partir de l accession du genome 
    et de l id de cette proteine, elle englobe la fonction best_hit qui lit le fichier resultat du blast
    et return le best_hit
    :param protein_id:
    :param accession_genome:
    :return: un dictionnaire avec comme clef (la position,Similarity,best_hit_id),les informations sur le meuilleur
    """
    def best_hit():
        result_blastp=save_to_blast(protein_id,assembly_accession)
        result_handle=open(result_blastp)
        blast_records = NCBIXML.parse(result_handle)
        best_hit = None
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if best_hit is None or hsp.expect < best_hit.hsps[0].expect:
                        best_hit = alignment
                        best_hit.hsps = [hsp]
        result_handle.close()
        return best_hit
    best_hit=best_hit()
    if best_hit:
        strand=best_hit.hsps[0].strand
        begin=best_hit.hsps[0].sbjct_start
        end=best_hit.hsps[0].sbjct_end
        similarity = best_hit.hsps[0].identities / best_hit.hsps[0].align_length * 100
        info_hit=best_hit.title
        best_hit_id=best_hit.hit_id.split("|")
        info_best_hit={"position":(begin,end,strand),"Similarity":similarity,
                       "best_hit_id":best_hit_id[1],"info_hit":info_hit}
        return info_best_hit
    else:
        return "aucun hit trouve"

def get_proteins_upstream_downstream(protein_id,assembly_accession):
    """
    Fonction qui renvoie ["Protein_Id","Gene_Id", "Gene_Name"] a partir d un seul de ces 4 critères
    les index en amont(upstream) et aval(downstream) d une proteine dans un tableau genome_df.
    :param protein_id:Protein_Id soit Gene_Id, soit Gene_Name de la protein d'interet
    :param assembly_accession:est l identifiant unique d'un genome donne(pour le tableau genome_df)
    :return: un dico contenant les valeurs list de index_upstream, et list denindex_downstream.
    """
    genome_df=read_table_tsv(assembly_accession)
    index = genome_df.loc[(genome_df["Protein_Id"] == protein_id) | 
                          (genome_df["Gene_Id"] == protein_id) | (genome_df["Gene_Name"] == protein_id)].index[0]
    index_upstream_pro = [] #valeur du dictionnaire
    index_downstream_pro = [] #valeur du dictionnaire

    for protein_id_name in ["Protein_Id","Gene_Id", "Gene_Name"]:
        index_upstream_pro.append(genome_df.loc[genome_df.index == index - 1][protein_id_name].iloc[0])
        index_downstream_pro.append(genome_df.loc[genome_df.index == index + 1][protein_id_name].iloc[0])  

    if genome_df.iloc[index].Strand == "-": # on inverse upstream, downstream si brin -
        index_upstream_pro, index_downstream_pro = index_downstream_pro, index_upstream_pro
    dico_upstream_downstream={"index_upstream_pro":index_upstream_pro,"index_downstream_pro":index_downstream_pro}
    return dico_upstream_downstream

def info_best_hit_upstream(protein_id,assembly_accession,subject_assembly_accession,my_os_type):
    protein_id_upstream=get_proteins_upstream_downstream(protein_id,assembly_accession)["index_upstream_pro"][0]
    # Lance le blastp pour la proteine en amont(upstream)
    run_blastp(protein_id_upstream,assembly_accession,subject_assembly_accession,my_os_type)
    #sauvegarde le resultat du blastp de la proteine upstream
    save_to_blast(protein_id_upstream,assembly_accession)
    #cherche les informations sur la meuilleur hit dans le resultat du blastp de upstream
    info_best_hit_upstream=info_best_hit(protein_id_upstream,assembly_accession)
    return info_best_hit_upstream
    
def info_best_hit_downstream(protein_id,assembly_accession,subject_assembly_accession,my_os_type):
    protein_id_downstream=get_proteins_upstream_downstream(protein_id,assembly_accession)["index_downstream_pro"][0]
    # Lance le blastp pour la proteine en amont(downstream)
    blastp_downstream=run_blastp(protein_id_downstream,assembly_accession,subject_assembly_accession,my_os_type)
    #sauvegarde le resultat du blastp de la proteine downstream
    save_to_blast(protein_id_downstream,assembly_accession)
    #cherche les informations sur la meuilleur hit dans le resultat du blastp de downstream
    info_best_hit_downstream=info_best_hit(protein_id_downstream,assembly_accession)
    return info_best_hit_downstream
    
    
def get_position(protein_id,assembly_accession):
    df_genome=read_table_tsv(assembly_accession)
    """
    Fonction qui a partir de id d'une protéine renvoie la position du gène qui code pour celle-ci sur le genome
    :param assembly_accession: identifiant unique associe a un genome donne
    :param protein_id:Protein_Id soit Gene_Id, soit Gene_Name
    :return: tuple position = (begin, end, strand)
    """
    line = df_genome.loc[(df_genome["Protein_Id"] == protein_id)|
                         (df_genome["Gene_Id"] == protein_id) | (df_genome["Gene_Name"] == protein_id)]
    if line.shape[0] != 0: # si le tableau n'est pas vide (il contient une ligne qui correspond à notre protéine)
        return line.Begin.iloc[0], line.End.iloc[0], line.Strand.iloc[0]
    else:
        raise Exception(f"La proteine {protein_id} ne se trouve pas dans le genome{assembly_accession}")
        
        
def get_synteny_table(protein_id,assembly_accession,subject_assembly_accession,my_os_type):
    """
    Fonction qui renvoie tableau de synténie sous forme de df et stocke ce df dans le dossier data/synteny_table
    :param protein_id:id de la protein d interet
    :param assembly_accession:est l’identifiant unique d'un génome_query
    :param subject_assembly_accession: est l identifiant unique d un génome_subject
    :param my_os_type:le systeme d exploitation
    :return: tableau de synténie avec pour colonne ["Protein, "begin, "End", "Strand", "Similarity"] 
    et pour lignes ["upstream","blastp protein en amont","protein_of_interest",
    "blastp protein d'interet","downstream","blastp proteine en aval"]
    """
    #Recupere sequence de la proteine d interet
    seq_protein=read_seq_protein(protein_id,assembly_accession)
   # Charge le tableau genome_df du genome d interet
    genome_df=read_table_tsv(assembly_accession)
     # Recupere l'index de la proteine d interet dans ce tableau
    index_protein=get_protein(protein_id,assembly_accession)
    #recupere la position du gene sur le genome dans ce tableau
    position_protein=get_position(protein_id,assembly_accession)
    # Recupere l'index de la proteine en amont dans ce tableau
    index_upstream_pro=get_proteins_upstream_downstream(protein_id,assembly_accession)["index_upstream_pro"]
    # Recupere l'index de la proteine en aval dans ce tableau
    index_downstream_pro=get_proteins_upstream_downstream(protein_id,assembly_accession)["index_downstream_pro"]
    # Lance le blastp de la proteine d interet
    run_blastp(protein_id,assembly_accession,subject_assembly_accession,my_os_type)
    #donne les infos du meuilleur hit de la proteine d interet
    info_hit=info_best_hit(protein_id,assembly_accession)
    #recuperation des infos du meuileur hit de la proteine en amont
    protein_id_upstream=get_proteins_upstream_downstream(protein_id,assembly_accession)["index_upstream_pro"][0]
    position_upstream=get_position(protein_id_upstream,assembly_accession)
    run_blastp(protein_id_upstream,assembly_accession,subject_assembly_accession,my_os_type)
    save_to_blast(protein_id_upstream,assembly_accession)
    info_best_hit_upstream=info_best_hit(protein_id_upstream,assembly_accession)
    #recuperation des infos du meuileur hit de la proteine en aval
    protein_id_downstream=get_proteins_upstream_downstream(protein_id,assembly_accession)["index_downstream_pro"][0]
    position_downstream=get_position(protein_id_downstream,assembly_accession)
    blastp_downstream=run_blastp(protein_id_downstream,assembly_accession,subject_assembly_accession,my_os_type)
    save_to_blast(protein_id_downstream,assembly_accession)
    info_best_hit_downstream=info_best_hit(protein_id_downstream,assembly_accession)
    position_blastp_protein=get_position(protein_id,subject_assembly_accession)
    protein_id_blastp_up=info_best_hit_upstream["best_hit_id"]
    position_blastp_up=get_position(protein_id_blastp_up,subject_assembly_accession)
    protein_id_blastp_down=info_best_hit_downstream["best_hit_id"]
    position_blastp_downstream=get_position(protein_id_blastp_down,subject_assembly_accession)
    #si tout marche bien on commence la creation du tableau
    #clonne du tableau DataFrame
    Relative_position=["upstream","blastp_up","protein_of_interest","blastp_prot","downstream","blastp_down"]
    #chacunes des termes representent un lignes du DataFrame
    Protein=pd.Series([index_upstream_pro,info_best_hit_upstream["best_hit_id"],
                       index_protein,info_hit["best_hit_id"],index_downstream_pro,
                       info_best_hit_downstream["best_hit_id"]],index=Relative_position)
    
    Begin=pd.Series(
        [position_upstream[0],position_blastp_up[0],
         position_protein[0],position_blastp_protein[0],position_downstream[0],
         position_blastp_downstream[0]],index=Relative_position)
    
    End=pd.Series([position_upstream[1],position_blastp_up[1],
                   position_protein[1],position_blastp_protein[1],position_downstream[1],
                   position_blastp_downstream[1]],index=Relative_position)
    
    Strand=pd.Series([position_upstream[2],position_blastp_up[2],
                      position_protein[2],position_blastp_protein[2],position_downstream[2],
                      position_blastp_downstream[2]],index=Relative_position)
    
    Similarity=pd.Series(["NaN",info_best_hit_upstream["Similarity"],"NaN",
                          info_hit["Similarity"],"NaN",info_best_hit_downstream
                          ["Similarity"]],index=Relative_position)
    #le tableau DataFrame nomme synteny_df
    synteny_df= pd.DataFrame({"Protein":Protein, "Begin": Begin, "End": End,"Strand":Strand,"Similarity":Similarity})
    #verifie que le dossier synteny_table existe
    folder =os.path.join(get_data_path(), "synteny_table")
    #si le chemin du dossier n existe pas
    if os.path.exists(folder) == False:
        #il cree le chemin
        os.mkdir(folder)
        #et stock le resultat du tableau syntenie_df dans le dossier syntenie_table sous format xlsx (excel)
    synteny_df.to_excel(os.path.join(folder,f"{protein_id}_{assembly_accession}_{subject_assembly_accession},.xlsx"))
    return synteny_df