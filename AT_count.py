#we turn on python
python

#we repeat lines 2, 5 and 8 from GC_count.py to import SeqIO, GC and pandas
from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd

#we import GC123 from Bio.SeqUtils - that counts GC not only in total, but also from the basis on the 1st, 2nd and 3rd position
from Bio.SeqUtils import GC123

#we load the first set of sequences - the whole mitochondrial genomes in the file genomes-for-AT_content.fa
filepath='/home/tom/Documents/Amblycera-mt-genomes/AT_content/genomes-for-AT_content.fa'

#we split it into separate sequences
seq_objects=SeqIO.parse(filepath,'fasta')

#from those sequences, we filter only the sequences themselves (without names of samples)
sequences=[seq for seq in seq_objects]

#we control the number of sequences
number_of_sequences=len(sequences)
print(number_of_sequences)

#now we make loop computing the AT contents, first part of the loop repeats from GC_count.py
for seq in sequences:
	#be VERY careful and responsible with the indention, it's an extremely IMPORTANT thing in python
	seq_id=seq.id
	sequence=seq.seq
	gc_content=GC(sequence)
	#we add computation of AT content as 100 - gc_content
	at_content=100 - gc_content
	#now we round the AT content for 2 decimal places
	at_content=round(at_content,2)
	#we print the rounded AT content
	print(seq_id,at_content)
	#and hit ENTER to see what will happen
	
	#WONDERFUL!!! It seems it works.
	
#now, let's use the same loop for computing AT contents also for 1st, 2nd and 3rd position (function GC123)
for seq in sequences:
	seq_id=seq.id
	sequence=seq.seq
	gc123_content=GC123(sequence)
	#first value of the chain will be GC content for the whole genomes, we will use it to count AT contents for the whole genomes
	at_content=100 - gc123_content[0]
	#second value will be the GC content of the first position
	at_content_1=100 - gc123_content[1]
	#and so on for the second and third position
	at_content_2=100 - gc123_content[2]
	at_content_3=100 - gc123_content[3]
	#now we can print it altogether
	print(seq_id,at_content,at_content_1,at_content_2,at_content_3)
	
#GREAT!!! It seems it works again, apart from that we forgot to round it. We can round it now.
at_content=round(at_content,2)
at_content_1=round(at_content_1,2)
at_content_2=round(at_content_2,2)
at_content_3=round(at_content_3,2)
#and try to print it again, it should look better now
print(seq_id,at_content,at_content_1,at_content_2,at_content_3)
#OOOOPS!!! It works only for the last sequence, to make it everywhere, we would need to repeat all the loop. Let's do it, but instead of printing, export the results directly to data rows:

#first, we must create the empty data rows
seq_ids=[]
at_contents=[]
at1_contents=[]
at2_contents=[]
at3_contents=[]

#the we do the computing loop
for seq in sequences:
	seq_id=seq.id
	sequence=seq.seq
	gc123_content=GC123(sequence)
	at_content=100 - gc123_content[0]
	at_content=round(at_content,2)
	at_content_1=100 - gc123_content[1]
	at_content_1=round(at_content_1,2)
	at_content_2=100 - gc123_content[2]
	at_content_2=round(at_content_2,2)
	at_content_3=100 - gc123_content[3]
	at_content_3=round(at_content_3,2)
	#and we append all computed AT contents to the data rows
	seq_ids.append(seq_id)
	at_contents.append(at_content)
	at1_contents.append(at_content_1)
	at2_contents.append(at_content_2)
	at3_contents.append(at_content_3)
	#now we print it to visualize the course of the process
	print(seq_id,at_content,at_content_1,at_content_2,at_content_3)
	#we hit ENTER and EVERYTHING WORKS!!!
	
### now we need to add an information about AT content of coding regions. The coding regions are in a separate file, so we must load that file first.
filepath='/home/tom/Documents/Amblycera-mt-genomes/AT_content/coding_regions.fa'
seq_objects=SeqIO.parse(filepath,'fasta')
sequences=[seq for seq in seq_objects]

#we create empty data row for the coding regions AT contents
at_coding_regions=[]

#and repeat the same computing loop as before
for seq in sequences:
	seq_id=seq.id
	sequence=seq.seq
	gc_content=GC(sequence)
	at_content=100 - gc_content
	at_content=round(at_content,2)
	at_coding_regions.append(at_content)
	print(seq_id,at_content)
	
	###Okey, after several trials it works, it is important that 2 different variables must always have 2 different names, and a data row cannot have the same name as a variable
	
#now, let's put all that to a nice dataframe. Firstly create an empty dataframe.
dataframe=pd.DataFrame()

#then load there all data rows we have
dataframe['Sequence_ID']=seq_ids

### OOOPPPSSS!!!! Something doesn't work!!!! It seems that python forgot the datarows from the first loop (with whole genomes). So we will need to repeat it again. Also, we discovered that the SeQIDs in the coding regions are slightly different, so we will let them for later.
filepath='/home/tom/Documents/Amblycera-mt-genomes/AT_content/genomes-for-AT_content.fa'
seq_objects=SeqIO.parse(filepath,'fasta')
sequences=[seq for seq in seq_objects]

seq_ids=[]
at_contents=[]
at1_contents=[]
at2_contents=[]
at3_contents=[]

for seq in sequences:
	seq_id=seq.id
	sequence=seq.seq
	gc123_content=GC123(sequence)
	at_content=100 - gc123_content[0]
	at_content=round(at_content,2)
	at_content_1=100 - gc123_content[1]
	at_content_1=round(at_content_1,2)
	at_content_2=100 - gc123_content[2]
	at_content_2=round(at_content_2,2)
	at_content_3=100 - gc123_content[3]
	at_content_3=round(at_content_3,2)
	seq_ids.append(seq_id)
	at_contents.append(at_content)
	at1_contents.append(at_content_1)
	at2_contents.append(at_content_2)
	at3_contents.append(at_content_3)
	print(seq_id,at_content,at_content_1,at_content_2,at_content_3)
	
#now we create empty dataframe and put there the values from this loop
dataframe=pd.DataFrame()
dataframe['Sequence_ID']=seq_ids
dataframe['AT_contents']=at_contents
dataframe['AT_contents_1']=at1_contents
dataframe['AT_contents_2']=at2_contents
dataframe['AT_contents_3']=at3_contents

#now we can try to print the dataframe
print(dataframe)

#or its shape
print(dataframe.shape)

#and export it to a csv file, we don't need to keep the orders (so index=False)
outputfile='/home/tom/Documents/Amblycera-mt-genomes/AT_content/at_content.csv'
dataframe.to_csv(outputfile,index=False)


#that's it. Now we have already fixed the coding regions file, so we can try to connect data from two files into one dataframe. First all the AT contents that we computed above:
filepath='/home/tom/Documents/Amblycera-mt-genomes/AT_content/genomes-for-AT_content.fa'
seq_objects=SeqIO.parse(filepath,'fasta')
sequences=[seq for seq in seq_objects]

seq_ids=[]
at_contents=[]
at1_contents=[]
at2_contents=[]
at3_contents=[]

for seq in sequences:
	seq_id=seq.id
	sequence=seq.seq
	gc123_content=GC123(sequence)
	at_content=100 - gc123_content[0]
	at_content=round(at_content,2)
	at_content_1=100 - gc123_content[1]
	at_content_1=round(at_content_1,2)
	at_content_2=100 - gc123_content[2]
	at_content_2=round(at_content_2,2)
	at_content_3=100 - gc123_content[3]
	at_content_3=round(at_content_3,2)
	seq_ids.append(seq_id)
	at_contents.append(at_content)
	at1_contents.append(at_content_1)
	at2_contents.append(at_content_2)
	at3_contents.append(at_content_3)
	print(seq_id,at_content,at_content_1,at_content_2,at_content_3)
	
#now the coding regions
filepath='/home/tom/Documents/Amblycera-mt-genomes/AT_content/coding_regions.fa'
seq_objects=SeqIO.parse(filepath,'fasta')
sequences=[seq for seq in seq_objects]

at_coding_regions=[]

for seq in sequences:
	seq_id=seq.id
	sequence=seq.seq
	gc_coding_rgs=GC(sequence)
	at_codings_rgs=100 - gc_coding_rgs
	at_codings_rgs=round(at_codings_rgs,2)
	at_coding_regions.append(at_codings_rgs)
	print(seq_id,at_codings_rgs)
	
#we create and empty dataframe
dataframe=pd.DataFrame()

#and load there all datarows in the order we want
dataframe['Sequence_ID']=seq_ids
dataframe['AT_contents']=at_contents
dataframe['AT_coding_regions']=at_coding_regions
dataframe['AT_contents_1']=at1_contents
dataframe['AT_contents_2']=at2_contents
dataframe['AT_contents_3']=at3_contents

#we try to print it
print(dataframe)

#HOOOORRRRAAAAAYYYYY!!!!!! It works! So we export it to a csv file.
outputfile='/home/tom/Documents/Amblycera-mt-genomes/AT_content/at_content.csv'
dataframe.to_csv(outputfile,index=False)
