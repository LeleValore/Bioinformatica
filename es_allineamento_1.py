# Esercizio 1:
# Date due sequenze T = AG e S = ACG, con i seguenti punteggi Match = +1, Mismatch = -1, Gap = -2.

# 1. Costruire la matrice di punteggio con stampa della matrice 
# 2. Eseguire traceBack con stampa della matrice di traceBack
# 3. Fornire l'allineamnto ottimale

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def build_matrice(seq1, seq2, match, mismatch, gap):
    """ Costruisce la matrice di punteggio per l'allineamento globale"""

    n = len(seq1) + 1
    m = len(seq2) + 1 

    # Inizializza la matrice

    matrice = [[0 for j in range(m)] for i in range(n)] 
    # Inizializza la prima riga e colonna con gli zeri come da definizione della matrice dei punteggi
    for i in range(n):
        matrice[i][0] = 0
        for j in range(m):
            matrice[0][j] = 0
    # Compila la matrice
    for i in range(1,n):
        for j in range(1,m):
            if seq1[i-1] == seq2[j-1]:
                score_diagonale = matrice[i-1][j-1] + match
            else:
                score_diagonale = matrice[i-1][j-1] + mismatch

            score_sopra = matrice[i-1][j] + gap
            score_sinistra = matrice[i][j-1] + gap

            matrice[i][j] = max(score_diagonale, score_sopra, score_sinistra, 0) # Aggiunto 0 per allineamento locale
    return matrice

def print_matrice(matrice):
    """Stampa la matrice di punteggio"""
    for row in matrice:
        print("\t".join(map(str, row)))

def trace_back(matrice, seq1, seq2, match, mismatch, gap):
    """"Esegue il trace back per trovare l'allineamento ottimale"""
    aligne_sequence1 = ""
    aligne_sequence2 = ""
    i = len(seq1)
    j = len(seq2)

    while i > 0 and j > 0:
        current_score = matrice[i][j]
        diagonal_score = matrice[i-1][j-1]
        score_sopra = matrice[i-1][j]
        score_sinistra = matrice[i][j-1]
        if current_score == diagonal_score + (match if seq1[i-1] == seq2[j-1]else mismatch):
            aligne_sequence1 = seq1[i-1] + aligne_sequence1
            aligne_sequence2 = seq2[j-1] + aligne_sequence2
            i -= 1 
            j -= 1
        elif current_score == score_sopra + gap:
            aligne_sequence1 = seq1[i-1] + aligne_sequence1
            aligne_sequence2 = "-" + aligne_sequence2
            i -= 1
        elif current_score == score_sinistra + gap:
            alignes_sequence1 = "-" + aligne_sequence1
            alignes_sequence2 = seq2[j-1] + aligne_sequence2
            j -= 1
        else:
            break
    return aligne_sequence1, aligne_sequence2


def print_traceback_matrice(matrice, seq1, seq2, match, mismatch, gap):
    """" 
    Stampa la matrice di trace back come un insimee di lettere dove:
    D = Diagonale,
    S = Sopra,
    L = Sinistra
    """

    n = len(seq1) + 1
    m = len(seq2) + 1 

    for i in range(1,n):
        row = []
        for j in range(1,m):
            current_score = matrice[i][j]
            diagonal_score = matrice[i-1][j-1]
            score_sopra = matrice[i-1][j]
            score_sinistra = matrice[i][j-1]

            if current_score == diagonal_score + (match if seq1[i-1] == seq2[j-1] else mismatch):
                row.append("D")
            elif current_score == score_sopra + gap:
                row.append("S")
            elif current_score == score_sinistra + gap:
                row.appendd("L")
            else:
                row.append("0")

        print("\t".join(row))




def main():

    """Definire due sequenze T, S con lettere A, C, G, T"""

    seq1 = "AGTAAT"
    seq2 = "TGAAAT"

    match = 1 
    missmatch = -1 
    gap = -2 

    #1. Costruire la matrice di punteggio con stampa della matrice
    matrice = build_matrice(seq1,seq2, match, missmatch, gap)
    print("Matrice di punteggio:")
    print_matrice(matrice)

    #2. Eseguire traceBack con stampa della matrice di traceBack
    print("\nMatrice di trace back:")
    print_traceback_matrice(matrice,seq1,seq2, match, missmatch,gap)
    aligned_seq1, aligned_seq2 = trace_back(matrice, seq1, seq2, match, missmatch, gap)

    # 3. Allineamento Ottimale
    print("\nAllineamento ottimale:")
    print(aligned_seq1)
    print(aligned_seq2)


if __name__ == "__main__":
    main()  