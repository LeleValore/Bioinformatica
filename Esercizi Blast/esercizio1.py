""" Eseguire l'algoritmo BLAST PROTEIN PROTEIN avendo la sequenza di query in formato FASTA"""
import ssl # per disabilitare la verifica SSL
ssl._create_default_https_context = ssl._create_unverified_context



from Bio.Blast import NCBIWWW
from Bio import SeqIO   

# Ho il formato FASTA in una variabile stringa
fasta_string = """ MAFPRPKRPAPAQEAATEGPSAASGVPQTGPGREVAATRPKTTKSGKALAKTRWVEPQNVVAAAAAKAKM
ATSIPEPEGAAAATAQHSAEPWARMGGKRTKKSKHLDDEYESSEEERETPAVPPTWRASQPSLTVRAQLA
PRPPMAPRSQIPSRHVLCLPPRNVTLLQERANKLVKYLMIKDYKKIPIKRADMLKDVIREYDEHFPEIIE
RATYTLEKKFGIHLKEIDKEEHLYILVCTRDSSARLLGKTKDTPRLSLLLVILGVIFMNGNRASEAVLWE
ALRKMGLRPGVRHPFLGDLRKLITDDFVKQKYLEYKKIPNSNPPEYEFLWGLRARHETSKMRVLRFIAQN
QNRDPREWKAHFLEAVDDAFKTMDVDMAEEHARAQMRARMNIGDEALIGRWSWDDIQVELLTWDEDGDFG
DAWARIPFAFWARYHQYILNSNRANRRATWRAGVSSGTNGGASTSVLDGPSTSSTIRTRNAARAGASFFS
WIQHR"""

# Scrivo la stringa in un file temporaneo
with open("temp_query.fasta", "w") as fasta_file:
    fasta_file.write(">query_sequence\n")
    fasta_file.write(fasta_string)  

# Leggo la sequenza dal file FASTA
record = SeqIO.read("temp_query.fasta", format="fasta") 
sequence = str(record.seq)
# Eseguo il BLAST
result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
# Salvo i risultati in un file
with open("blast_results.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()   
print("BLAST eseguito e risultati salvati in 'blast_results.xml'")

# Rimuovo il file temporaneo
import os
os.remove("temp_query.fasta")   
