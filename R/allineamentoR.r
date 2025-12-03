### Allinemaento di sequenze in R: 

# Per fare l'allineamento in R dobbiamo utilizzare il pacchetto pwalign di Bioconductor.
# Se non lo avete ancora installato, potete farlo con i seguenti comandi:


# Assicurati che BiocManager sia installato e caricato
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# Esegui l'installazione di Biostrings

# Crea una directory della libreria personale (se non esiste)
.libPaths("~/R/library")
# Riprova l'installazione
BiocManager::install("Biostrings")
BiocManager::install("pwalign")

# Carichiamo il pacchetto
library(pwalign)
library(Biostrings) 

# Creiamo due sequenze di esempio
seq1 <- DNAString("AGCTTAGCTAGCTAGCTAGC")
seq2 <- DNAString("AGCTAGCTAGCTAGCTAGC")    
# Eseguiamo l'allineamento globale utilizzando l'algoritmo di Needleman-Wunsch
alignment <- pairwiseAlignment(seq1, seq2, type = "global")
# Visualizziamo il risultato dell'allineamento
print(alignment)