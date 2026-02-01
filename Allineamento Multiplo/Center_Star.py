from Bio import Align
import numpy as np

def center_star_msa(sequences):
    # 1. Configurazione dell'aligner (Parametri simili alle slide)
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2

    n = len(sequences)
    
    # 2. Calcolo della Matrice delle Distanze (Score)
    # Troviamo la sequenza centrale Sc che massimizza la somma dei punteggi
    scores = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            score = aligner.score(sequences[i], sequences[j])
            scores[i, j] = scores[j, i] = score

    # Scelta di Sc
    center_idx = np.argmax(scores.sum(axis=1))
    center_seq = sequences[center_idx]
    
    # 3. Costruzione dell'MSA (Inizia con la sequenza centrale)
    # msa_list conterrà le sequenze come liste di caratteri
    msa_list = [list(center_seq)]
    others = [i for i in range(n) if i != center_idx]
    
    for idx in others:
        query_seq = sequences[idx]
        # Allineamento tra Sc (originale) e la nuova sequenza
        alignment = aligner.align(center_seq, query_seq)[0]
        
        # Estraiamo le stringhe allineate (es. 'CCTGCT-GCAG' e 'GATG-TG-CAG')
        aligned_c, aligned_q = alignment.target, alignment.query
        # Nota: In versioni recenti di Biopython usiamo l'indice degli allineamenti
        c_str, q_str = format(alignment, "fasta").split('\n')[1], format(alignment, "fasta").split('\n')[3]

        # 4. Propagazione dei Gap ("Once a gap, always a gap")
        new_msa = []
        new_q = []
        p_aln = 0 # Puntatore per l'allineamento corrente
        p_msa = 0 # Puntatore per l'MSA esistente

        while p_aln < len(c_str) or p_msa < len(msa_list[0]):
            char_c = c_str[p_aln] if p_aln < len(c_str) else None
            char_m = msa_list[0][p_msa] if p_msa < len(msa_list[0]) else None
            
            if char_c == '-' and char_m == '-':
                # Entrambi hanno un gap: aggiungi normalmente
                for i in range(len(msa_list)):
                    if len(new_msa) <= i: new_msa.append([])
                    new_msa[i].append(msa_list[i][p_msa])
                new_q.append(q_str[p_aln])
                p_aln += 1; p_msa += 1
            elif char_c == '-':
                # Nuovo gap in Sc: aggiungilo a TUTTE le sequenze nell'MSA
                for i in range(len(msa_list)):
                    if len(new_msa) <= i: new_msa.append([])
                    new_msa[i].append('-')
                new_q.append(q_str[p_aln])
                p_aln += 1
            elif char_m == '-':
                # Gap già presente nell'MSA: aggiungilo alla nuova sequenza
                for i in range(len(msa_list)):
                    if len(new_msa) <= i: new_msa.append([])
                    new_msa[i].append(msa_list[i][p_msa])
                new_q.append('-')
                p_msa += 1
            else:
                # Nessun gap: aggiungi i caratteri correnti
                for i in range(len(msa_list)):
                    if len(new_msa) <= i: new_msa.append([])
                    new_msa[i].append(msa_list[i][p_msa])
                new_q.append(q_str[p_aln])
                p_aln += 1; p_msa += 1
        
        msa_list = new_msa + [new_q]

    return ["".join(seq) for seq in msa_list], center_idx

# Esempio con le sequenze delle slide
seqs = ["CCTGCTGCAG", "GATGTGCCG", "GATGTGCAG", "CCGCTAGCAG", "CCTGTAGG"]
risultato, c_idx = center_star_msa(seqs)

print(f"Sequenza centrale scelta: S{c_idx+1}")
for i, s in enumerate(risultato):
    print(f"S_allineata: {s}")