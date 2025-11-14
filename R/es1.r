esponente <- function(base, exp = 3 ){
  risultato <- base ^ exp
  return(risultato)
}

somma <- function(a,b){
  return(a+b)
}


# Creazione di un vettore di sequenze di DNA
sequences <- c("ATGCGTA", "CGTACGTAGC", "TTAGGC", "GGCAT")
# Calcolo della lunghezza di ciascuna sequenza
lengths <- nchar(sequences)
# Stampa delle lunghezze
print(lengths)


for (i in 1:length(sequences)) {
  cat("La lunghezza della sequenza", sequences[i], "è", lengths[i], "nucleotidi.\n")
}


while(lengths[1] < 10) {
  sequences[1] <- paste0(sequences[1], "A")
  lengths[1] <- nchar(sequences[1])
  cat("La nuova lunghezza della prima sequenza è", lengths[1], "nucleotidi.\n")
}

print(esponente(5)) # 5 ^ 3 = 125
print(somma(3.5, 6)) # capisci i tipi di dato in R 


###VETTORI 

vet1 <- c(2,3,4,5)

print(vet1) # print del vettore

# posso anche stampare un singolo elemento del vettore 
print(vet1[2]) # stampa il secondo elemento del vettore
print(vet1[1:3]) # stampa i primi tre elementi del vettore
print(length(vet1)) # stampa la lunghezza del vettore


# posso anche sommare due vettori
vet2 <- c(10,20,30,40)
somma_vettori <- vet1 + vet2
print(somma_vettori) # stampa la somma dei due vettori



# posso anche utilizzare indici negatici per escludere elementi
print(vet1[-2]) # stampa il vettore senza il secondo elemento


# posso anche inserire elementi in un vettore
vet1 <- c(vet1, 6) # aggiunge 6 alla fine del vettore
print(vet1) # stampa il vettore aggiornato    

vet1 <- c(1, vet1) # aggiunge 1 all'inizio del vettore
print(vet1) # stampa il vettore aggiornato  

# posso anche rimuovere elementi da un vettore
vet1 <- vet1[-3] # rimuove il terzo elemento del vettore
print(vet1) # stampa il vettore aggiornato

# posso anche modificare elementi di un vettore
vet1[2] <- 10 # modifica il secondo elemento del vettore
print(vet1) # stampa il vettore aggiornato  


# Tutte le operazioni sui vettori sono le seguenti:
# Creazione: c(), seq(), rep()
# Accesso: []
# Modifica: []
# Aggiunta: c()
# Rimozione: []
# Lunghezza: length() 
# Operazioni aritmetiche: +, -, *, /, ^
# Funzioni utili: sum(), mean(), median(), sd(), var()    
# Funzioni di ordinamento: sort(), order()
# Funzioni di ricerca: which(), match(), %in%
# Funzioni di manipolazione: rev(), unique(), intersect(), union(), setdiff()
# Funzioni di aggregazione: tapply(), aggregate()
# Funzioni di trasformazione: as.numeric(), as.character(), as.factor()
# Funzioni di combinazione: cbind(), rbind()
# Funzioni di subsetting: subset()



###MATRICI 
matrice1 <- matrix(1:9, nrow=3, ncol=3, byrow=FALSE)
print(matrice1) # stampa la matrice
# posso anche accedere a singoli elementi della matrice
print(matrice1[2,3]) # stampa l'elemento alla seconda riga e terza colonna
print(matrice1[1,]) # stampa la prima riga
print(matrice1[,2]) # stampa la seconda colonna 
# posso anche modificare elementi della matrice
matrice1[3,1] <- 10 # modifica l'elemento alla terza riga e prima colonna
print(matrice1) # stampa la matrice aggiornata
# posso anche aggiungere righe o colonne alla matrice
matrice1 <- rbind(matrice1, c(11,12,13)) # aggiunge una nuova riga alla matrice
print(matrice1) # stampa la matrice aggiornata
matrice1 <- cbind(matrice1, c(14,15,16,17)) # aggiunge una nuova colonna alla matrice
print(matrice1) # stampa la matrice aggiornata
# posso anche rimuovere righe o colonne dalla matrice
matrice1 <- matrice1[-2, ] # rimuove la seconda riga della matrice
print(matrice1) # stampa la matrice aggiornata
matrice1 <- matrice1[, -3] # rimuove la terza colonna della matrice
print(matrice1) # stampa la matrice aggiornata      
# Tutte le operazioni sulle matrici sono le seguenti:
# Creazione: matrix()
# Accesso: [, ]
# Modifica: [, ]
# Aggiunta: rbind(), cbind()          
# Rimozione: [, -]
# Dimensioni: dim(), nrow(), ncol()
# Operazioni aritmetiche: +, -, *, /, ^
# Funzioni utili: rowSums(), colSums(), rowMeans(), colMeans()    
# Funzioni di ordinamento: apply(), order()
# Funzioni di ricerca: which(), match(), %in%
# Funzioni di manipolazione: t(), diag(), solve()
# Funzioni di aggregazione: aggregate()
# Funzioni di trasformazione: as.data.frame(), as.matrix()
# Funzioni di combinazione: rbind(), cbind()
# Funzioni di subsetting: subset()  

###LISTE
lista1 <- list(nome="Emanuele", età=30, voti=c(28,30,27), matrice=matrix(1:4, nrow=2))
print(lista1) # stampa la lista
# posso anche accedere a singoli elementi della lista
print(lista1$nome) # stampa l'elemento "nome"

print(lista1[[3]]) # stampa il terzo elemento della lista


# posso anche modificare elementi della lista
lista1$età <- 31 # modifica l'elemento "età"
print(lista1) # stampa la lista aggiornata



# posso anche aggiungere elementi alla lista
lista1$città <- "Roma" # aggiunge un nuovo elemento "città"
print(lista1) # stampa la lista aggiornata



# posso anche rimuovere elementi dalla lista
lista1$matrice <- NULL # rimuove l'elemento "matrice"
print(lista1) # stampa la lista aggiornata





# Tutte le operazioni sulle liste sono le seguenti:
# Creazione: list()
# Accesso: $, [[ ]]
# Modifica: $, [[ ]]
# Aggiunta: $
# Rimozione: $ <- NULL
# Lunghezza: length()
# Funzioni utili: names(), unlist(), lapply(), sapply()
# Funzioni di ordinamento: order()
# Funzioni di ricerca: which(), match(), %in%        
# Funzioni di manipolazione: c(), append()
# Funzioni di aggregazione: rbind(), cbind()
# Funzioni di trasformazione: as.data.frame(), as.list()
# Funzioni di combinazione: c()
# Funzioni di subsetting: subset()
l <- list(a=1, b=2, c=3)
print(l[1])     # stampa la prima componente come lista
print(l[[1]])   # stampa il valore della prima componente
