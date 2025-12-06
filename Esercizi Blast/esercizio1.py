""" Eseguire l'algoritmo BLAST basato sul tipo di sequenza specificato
    In questo file .py possiamo direttamente inserire il FASTA della Proteina o Nucleotide
    e verrà eseguito il BLAST appropriato (blastp, blastn, ecc.)
    I risultati verranno salvati in un file XML chiamato 'blast_results_{tipo}.xml'
"""
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from Bio.Blast import NCBIWWW
from Bio import SeqIO
import os

# Dizionario con le sequenze di test per <MAGED4 in Homo sapiens>
sequences_MAGED4 = {
    "protein": """ MAFPRPKRPAPAQEAATEGPSAASGVPQTGPGREVAATRPKTTKSGKALAKTRWVEPQNVVAAAAAKAKM
ATSIPEPEGAAAATAQHSAEPWARMGGKRTKKSKHLDDEYESSEEERETPAVPPTWRASQPSLTVRAQLA
PRPPMAPRSQIPSRHVLCLPPRNVTLLQERANKLVKYLMIKDYKKIPIKRADMLKDVIREYDEHFPEIIE
RATYTLEKKFGIHLKEIDKEEHLYILVCTRDSSARLLGKTKDTPRLSLLLVILGVIFMNGNRASEAVLWE
ALRKMGLRPGVRHPFLGDLRKLITDDFVKQKYLEYKKIPNSNPPEYEFLWGLRARHETSKMRVLRFIAQN
QNRDPREWKAHFLEAVDDAFKTMDVDMAEEHARAQMRARMNIGDEALIGRWSWDDIQVELLTWDEDGDFG
DAWARIPFAFWARYHQYILNSNRANRRATWRAGVSSGTNGGASTSVLDGPSTSSTIRTRNAARAGASFFS
WIQHR""",
    "nucleotide": "ATGAAATTCCCACGACCAAAACGCCCCGCACCCGCACAGGAAGCCGCCACAGAGCCGCCCAGCGCATCC",
}
sequence_P53 = {
    "protein": """MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA
PRVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKT
CPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRN
TFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACAGR
DRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALEL
KDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
"""
}
sequence_P73 = {
    'protein': """MAQSTATSPDGGTTFEHLWSSLEPDSTYFDLPQSSRGNNEVVGGTDSSMDVFHLEGMTTSVMAQFNLLSS
TMDQMSSRAASASPYTPEHAASVPTHSPYAQPSSTFDTMSPAPVIPSNTDYPGPHHFEVTFQQSSTAKSA
TWTYSPLLKKLYCQIAKTCPIQIKVSTPPPPGTAIRAMPVYKKAEHVTDVVKRCPNHELGRDFNEGQSAP
ASHLIRVEGNNLSQYVDDPVTGRQSVVVPYEPPQVGTEFTTILYNFMCNSSCVGGMNRRPILIIITLEMR
DGQVLGRRSFEGRICACPGRDRKADEDHYREQQALNESSAKNGAASKRAFKQSPPAVPALGAGVKKRRHG
DEDTYYLQVRGRENFEILMKLKESLELMELVPQPLVDSYRQQQQLLQRPSHLQPPSYGPVLSPMNKVHGG
MNKLPSVNQLVGQPPPHSSAATPNLGPVGPGMLNNHGHAVPANGEMSSSHSAQSMVSGSHCTPPPPYHAD
PSLVSFLTGLGCPNCIEYFTSQGLQSIYHLQNLTIEDLGALKIPEQYRMTIWRGLQDLKQGHDYSTAQQL
LRSSNAATISIGGSGELQRQRVMEAVHFRVRHTITIPNRGGPGGGPDEWADFGFDLPDCKARKQPIKEEF
TEAEIH
"""
}
sequence_BAX = {
    'protein': """DMFSDGNFNWVRVVALFYFAS"""
}

# Mapping tra tipo sequenza e parametri BLAST
blast_config = {
    "protein": {
        "program": "blastp",
        "database": "nr",
        "description": "Protein-Protein BLAST"
    },
    "nucleotide": {
        "program": "blastn",
        "database": "nt",
        "description": "Nucleotide-Nucleotide BLAST"
    },
}

def run_blast(sequence_type="protein"):
    """
    Esegue BLAST basato sul tipo di sequenza specificato.
    
    Args:
        sequence_type (str): Tipo di sequenza - 'protein' o 'nucleotide'
    """
    #NB: Se cambi tipo di sequenza, cambia anche il nome della variabile "sequence_BAX..ecc.. basandoti 
    # sul dizionario definito sopra)
    print(f"Avvio BLAST per il tipo di sequenza: {sequence_type}")
    
    # Valida il tipo di sequenza
    if sequence_type not in sequence_BAX:
        print(f"Errore: tipo di sequenza '{sequence_type}' non valido.")
        print(f"Tipi disponibili: {list(sequence_BAX.keys())}")
        return
    
    fasta_string = sequence_BAX[sequence_type]
    config = blast_config[sequence_type]
    
    print(f"\n{'='*60}")
    print(f"Esecuzione: {config['description']}")
    print(f"Program: {config['program']} | Database: {config['database']}")
    print(f"{'='*60}\n")
    
    # Scrivo la sequenza in un file temporaneo
    temp_file = f"temp_query_{sequence_type}.fasta"
    with open(temp_file, "w") as fasta_file:
        fasta_file.write(f">{sequence_type}_sequence\n")
        fasta_file.write(fasta_string)
    
    # Leggo la sequenza dal file FASTA
    record = SeqIO.read(temp_file, format="fasta")
    sequence = str(record.seq)
    
    # Eseguo BLAST
    try:
        print("Invio della richiesta BLAST ai server NCBI...")
        result_handle = NCBIWWW.qblast(config["program"], config["database"], sequence)
        
        # Salvo i risultati in un file XML
        output_file = f"blast_results_BAX_{sequence_type}.xml"
        with open(output_file, "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()
        
        print(f"✓ BLAST completato!")
        print(f"✓ Risultati salvati in '{output_file}'")
        
    except Exception as e:
        print(f"✗ Errore durante l'esecuzione di BLAST: {e}")
    
    finally:
        # Rimuovo il file temporaneo
        if os.path.exists(temp_file):
            os.remove(temp_file)

# Esegui il BLAST desiderato
if __name__ == "__main__":
    # Specifica il tipo di sequenza che vuoi analizzare
    # Scegli tra: "protein" o "nucleotide"
    sequence_type = "protein"  
    
    run_blast(sequence_type)