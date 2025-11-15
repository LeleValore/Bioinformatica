# Funzione lapply: Applicazione di una funzione a ogni elemento di una lista
# Creazione di una lista di vettori numerici
lista_numeri <- list(a = 1:5, b = 6:10, c = 11:15)  
# Funzione per calcolare la somma di un vettore
calcola_somma <- function(x) {
  return(sum(x))
}
print(lista_numeri)
# Applicazione della funzione a ogni elemento della lista usando lapply
somme <- lapply(lista_numeri, calcola_somma)
# Stampa dei risultati
print(somme)    

#Funzione sapply(): Simile a lapply ma restituisce un vettore o una matrice
# Creazione di una lista di vettori numerici
lista_numeri2 <- list(a = 1:5, b = 6:10, c = 11:15)  
# Applicazione della funzione a ogni elemento della lista usando sapply
somme2 <- sapply(lista_numeri2, calcola_somma)
# Stampa dei risultati
print(somme2)   





# Data-Frame : Una delle strutture dati più comuni in R per memorizzare dati tabulari è una tabella bidimensionale per rappresentare dati.
# le righe rappresentano le osservazioni le colonne rappresentano le variabili

# si utilizza la funzione data.frame() per creare un data frame

# Creazione di un data frame
dati <- data.frame(
  Nome = c("Alice", "Bob", "Charlie"),
  Età = c(25, 30, 35),
  Altezza = c(165, 180, 175)
)   


# Stampa del data frame
print(dati)



# Accesso a una colonna specifica
print(dati$Nome) # stampa la colonna "Nome"
print(dati[, "Età"]) # stampa la colonna "Età"


# Accesso a una riga specifica
print(dati[2, ]) # stampa la seconda riga


# Accesso a un elemento specifico
print(dati[1, "Altezza"]) # stampa l'altezza di Alice



# Aggiunta di una nuova colonna
dati$Peso <- c(55, 75, 68) # aggiunge la colonna "Peso"
print(dati) # stampa il data frame aggiornato




# Aggiunta di una nuova riga
nuova_riga <- data.frame(Nome = "David", Età = 28, Altezza = 172, Peso = 70)
dati <- rbind(dati, nuova_riga) # aggiunge la nuova riga
print(dati) # stampa il data frame aggiornato   



# Modifica di un elemento specifico
dati[2, "Età"] <- 31 # modifica l'età di Bob
print(dati) # stampa il data frame aggiornato



# Rimozione di una colonna
dati$Peso <- NULL # rimuove la colonna "Peso"
print(dati) # stampa il data frame aggiornato


# Funzioni rbind() e cbind(): Aggiunta di righe e colonne a matrici rbind() aggiunge righe, cbind() aggiunge colonne
# Creazione di una matrice
matrice1 <- matrix(1:12, nrow=4, ncol=3)
print(matrice1) # stampa la matrice     
# Aggiunta di una nuova riga
nuova_riga <- c(13,14,15)
matrice1 <- rbind(matrice1, nuova_riga) # aggiunge la nuova riga
print(matrice1) # stampa la matrice aggiornata
# Aggiunta di una nuova colonna
nuova_colonna <- c(16,17,18,19,20)
matrice1 <- cbind(matrice1, nuova_colonna) # aggiunge la nuova colonna
print(matrice1) # stampa la matrice aggiornata  


####Grafici in R : Creazione di grafici semplici utilizzando le funzioni di base di R

# Creazione di un vettore di dati
x <- c(1, 2, 3, 4, 5)
y <- c(2, 3, 5, 7, 11)


# Creazione di un grafico a dispersione (scatter plot)
plot(x, y, main="Grafico a dispersione", xlab="Asse X", ylab="Asse Y", col="blue", pch=19)



# Aggiunta di una linea di regressione
model <- lm(y ~ x)
abline(model, col="red")    




# Creazione di un istogramma
dati_istogramma <- c(1,2,2,3,3,3,4,4,4,4,5,5,5,5,5) 
hist(dati_istogramma, main="Istogramma", xlab="Valori", ylab="Frequenza", col="lightgreen", border="black")



# Creazione di un grafico a barre
categorie <- c("A", "B", "C", "D")
valori <- c(10, 15, 7, 20)
barplot(valori, names.arg=categorie, main="Grafico a barre", xlab="Categorie", ylab="Valori", col="orange", border="blue")


# Creazione di un grafico a torta
fette <- c(10, 20, 30, 40)
etichette <- c("Fetta A", "Fetta B", "Fetta C", "Fetta D")
pie(fette, labels=etichette, main="Grafico a torta", col=rainbow(length(fette)))        



# Funzioni lines() e points(): Aggiunta di linee e punti a grafici esistenti

# Creazione di un grafico vuoto
plot(1:10, 1:10, type="n", main="Grafico con lines() e points()", xlab="Asse X", ylab="Asse Y")             
# Aggiunta di una linea
lines(1:10, (1:10)^2, col="blue", lwd=2) # linea blu spessa
# Aggiunta di punti
points(1:10, (1:10)^2, col="red", pch=19) # punti rossi solidi       

#Conteggio dati nel dataset con barplot:
dati <- c("A", "B", "A", "C", "B", "A", "D", "C", "B", "A")
conteggio <- table(dati) # conta le occorrenze di ogni categoria
barplot(conteggio, main="Conteggio delle categorie", xlab="Categorie", ylab="Frequenza", col="lightblue")           


