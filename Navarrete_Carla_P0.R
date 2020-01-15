#1. A partir de las siguientes secuencias concatenadas de RNA que te
#envío como archivo adjunto en formato fasta encuentra (usando librerías
#especializadas) la secuencia correspondiente de aminoácidos.

library(Biostrings)
seq_RNA<-readRNAStringSet("CARLA/first (1).fasta")
seq_RNA
seq_AA<-translate(seq_RNA)
seq_AA

#2. Escoge dos problemas de la plataforma Rosalind de entre los siguientes: 
#Counting DNA Nucleotides <http://rosalind.info/problems/dna/>
secuencia<-DNAStringSet("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC")
secuencia
alphabetFrequency(secuencia, baseOnly=TRUE)


#Complementing a Strand of DNA <http://rosalind.info/problems/revc/>,
#Computing GC Content <http://rosalind.info/problems/gc/>

Rosalind1<-DNAStringSet("CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG")
Rosalind2<-DNAStringSet("CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC")
Rosalind3<-DNAStringSet("CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")

Todos_juntos<-c(Rosalind1, Rosalind2, Rosalind3)
Todos_juntos
nombres<-c("Rosalind_6404", "Rosalind_5959", "Rosalind_0808")
names(Todos_juntos)<-nombres
Todos_juntos
contenido_GC<-letterFrequency(Todos_juntos, letters="CG")
contenido_GC

for(i in 1:length(Todos_juntos)){
  print(paste("el nombre de la secuencia es", nombres, "y el GC% es", contenido_GC, "%"))
