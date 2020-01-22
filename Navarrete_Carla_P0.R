#1. A partir de las siguientes secuencias concatenadas de RNA que te
#envío como archivo adjunto en formato fasta encuentra (usando librerías
#especializadas) la secuencia correspondiente de aminoácidos.

library(Biostrings)
seq_RNA<-readRNAStringSet("CARLA/first (1).fasta")
seq_RNA
seq_AA<-translate(seq_RNA)
seq_AA #HAPPY NEW YEAR


#2. Escoge dos problemas de la plataforma Rosalind de entre los siguientes: 

#PROBLEMA 1: Counting DNA Nucleotides <http://rosalind.info/problems/dna/>

#Solucion 1 (con libreria: usando alphabetFrequency)
secuencia<-DNAString("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC")
secuencia
alphabetFrequency(secuencia, baseOnly=TRUE) #Determinarla cantidad de A, T, C y G de la secuencia
#baseOnly=TRUE hace que solo te arroje el no. que hay de cada nucleótido y no de todas las letras de aa  

#Solución 2 (con libreria: usando matchPattern)
m1<-length(matchPattern("A", secuencia)) #Determinar con "lenght" la cantidad de letras en la secuencia que coincidan (matchPattern) con A
m1
m2<-length(matchPattern("C", secuencia)) #Determinar la cantidad de letras en la secuencia que coincidan con C
m2
m3<-length(matchPattern("G", secuencia)) #Determinar la cantidad de letras en la secuencia que coincidan con G
m3
m4<-length(matchPattern("T", secuencia)) #Determinar la cantidad de letras en la secuencia que coincidan con T
m4

tabla<-c(A=m1, C=m2, G=m3, T=m4) #ordenar los resultados en una tabla
tabla

#Solucion 3 (sin libreria especializada)
sec<-c("A","G","C","T","T","T","T","C","A","T","T","C","T","G","A","C","T","G","C","A","A","C","G","G","G","C","A","A","T","A","T","G","T","C","T","C","T","G","T","G","T","G","G","A","T","T","A","A","A","A","A","A","A","G","A","G","T","G","T","C","T","G","A","T","A","G","C","A","G","C")
sec

A<-"A"
A<-countMatches(A, sec)
T<-"T"
T<-countMatches(T, sec)
C<-"C"
C<-countMatches(C, sec)
G<-"G"
G<-countMatches(G, sec)

todos<-c("A"=A, "T"=T, "G"=G, "C"=C)
todos

#PROBLEMA 2: Complementing a Strand of DNA <http://rosalind.info/problems/revc/>

#Solucion 1 (con libreria: usando reverseComplement)
seq_1<-DNAString("AAAACCCGGT")
seq_1
reverseComplement(seq_1) #obtener la secuencia reverso complemento de seq_1

#Solución 2 (sin libreria especializada)
#Cambiar las T por A de la secuencia original seq_1
seq_2<-gsub("A", "T", seq_1)
seq_2

#Cambiar las C por G de la secuencia nueva seq_2
seq_3<-gsub("C", "G", seq_2)
seq_3<-DNAString(seq_3)
length(seq_3)

#Cambiar las G por C de la secuencia nueva seq_3
seq_3[8:9]="C"
seq_3

#Cambiar las T por A de la secuencia nueva seq_3
seq_3[10]="A"
seq_3 #esta es la secuencia complemento

#Voltear la secuencia
rev_comp1<-reverse(seq_3) #puedes usar la funcion reverse (no es de Biostrings)
rev_comp1 #esta es la secuencia reverso complemento

rev_comp2<-seq_3[10:1] # o tambien puedes listar las letras de la secuencia de atras hacia adelante
rev_comp2 #esta es la secuencia reverso complemento


#PROBLEMA 3: Computing GC Content <http://rosalind.info/problems/gc/>

#Solucion 1 (con libreria: usando letterFrequency)
Rosalind1<-DNAStringSet("CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG")
Rosalind2<-DNAStringSet("CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC")
Rosalind3<-DNAStringSet("CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")

Todos_juntos<-c(Rosalind1, Rosalind2, Rosalind3)
Todos_juntos
nombres<-c("Rosalind_6404", "Rosalind_5959", "Rosalind_0808")
names(Todos_juntos)<-nombres
Todos_juntos #las 3 secuencias y su nombre en un solo objeto

for(i in 1:length(Todos_juntos)){print(letterFrequency(Todos_juntos[i], letters = "GC")/width(Todos_juntos[i])*100)
  } #ciclo for que permite calcular el GC% de cada secuencia del objeto "Todos_juntos"

contenido_GC<-c(53.75, 53.57, 60.91) #resultados del GC% de cada secuencia

if(contenido_GC[1]>contenido_GC[3]&contenido_GC[2]){
  print(paste("La secuencia", nombres[1], "posee el mayor contenido GC, que tiene un valor de", contenido_GC[1],"%"))
} else if(contenido_GC[2]>contenido_GC[3]&contenido_GC[1]){
    print(paste("La secuencia", nombres[2], "posee el mayor contenido GC, que tiene un valor de", contenido_GC[2],"%"))
} else if(contenido_GC[3]>contenido_GC[2]&contenido_GC[1]){
    print(paste("La secuencia", nombres[3], "posee el mayor contenido GC, que tiene un valor de", contenido_GC[3],"%"))
  } #te arroja el nombre de la secuencia con mayor GC%

#Solucion 2 (con libreria: usando vcountPattern)
GC_R1<-((vcountPattern("C", Rosalind1)+vcountPattern("G", Rosalind1))/width(Rosalind1))*100
GC_R1
GC_R2<-((vcountPattern("C", Rosalind2)+vcountPattern("G", Rosalind2))/width(Rosalind2))*100
GC_R2
GC_R3<-((vcountPattern("C", Rosalind3)+vcountPattern("G", Rosalind3))/width(Rosalind3))*100
GC_R3

if(GC_R1>GC_R3&GC_R2){print(paste("La secuencia Rosalind_6404 posee el mayor contenido GC, que tiene un valor del 53.75%"))
  } else if (GC_R2>GC_R3&GC_R1){print(paste("La secuencia Rosalind_5959 posee el mayor contenido GC, que tiene un valor del 53.57%"))
  } else if (GC_R3>GC_R2&GC_R1){print(paste("La secuencia Rosalind_0808 posee el mayor contenido GC, que tiene un valor del 60.91%"))
  }

#Solucion 3 (sin libreria especializada)
Rosalind1<-c("C","C","T","G","C","G","G","A","A","G","A","T","C","G","G","C","A","C","T","A","G","A","A","T","A","G","C","C","A","G","A","A","C","C","G","T","T","T","C","T","C","T","G","A","G","G","C","T","T","C","C","G","G","C","C","T","T","C","C","C","T","C","C","C","A","C","T","A","A","T","A","A","T","T","C","T","G","A","G","G")
x<-"C"
x<-countMatches(x, Rosalind1)
y<-"G"
y<-countMatches(y, Rosalind1)
cont_GC<-(x+y)/length(Rosalind1)*100

Rosalind2<-c("CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC")
#haciendo lo mismo el resultado sería 53.57

Rosalind3<-c("CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")
#haciendo lo mismo el resultado sería 60.91

nombres<-c("Rosalind_6404", "Rosalind_5959", "Rosalind_0808")
contenido_GC<-c(53.75, 53.57, 60.91)

if(contenido_GC[1]>contenido_GC[3]&contenido_GC[2]){
  print(paste("La secuencia", nombres[1], "posee el mayor contenido GC, que tiene un valor de", contenido_GC[1],"%"))
} else if(contenido_GC[2]>contenido_GC[3]&contenido_GC[1]){
  print(paste("La secuencia", nombres[2], "posee el mayor contenido GC, que tiene un valor de", contenido_GC[2],"%"))
} else if(contenido_GC[3]>contenido_GC[2]&contenido_GC[1]){
  print(paste("La secuencia", nombres[3], "posee el mayor contenido GC, que tiene un valor de", contenido_GC[3],"%"))
}
