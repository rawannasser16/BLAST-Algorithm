sequence = open('sequence.fasta', 'r')
# step2
unneeded_var = "RLLA"
def getseq():
    for line in sequence:
     if line[0] != '>':
        new_sequence = line.replace("RLLA", "XXXX")
        print(new_sequence)
    return new_sequence    
      
# step3
arr = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z',
       'X']

f = open('blosum.txt', 'r')
x = f.readlines()
f.close()
l = []
h = [""]

for i in x:
    singleRow = i.split(" ")
    for j in singleRow:
        if j != "" and j != "\n":
            h.append(j)
    l.append(h)
    h = []


def Ditionary():
    list = []
    Dic = {}
    for a in word:  # 1st word
        for b in a:  # 1st letter in 1st word
            for c in arr:  # 1st aminoacid
                new_a = a.replace(b, c)  # PQG replace p with each amino acid c in arr (keeping rest of word as it is)
                # calculate score of new_a compared to a (original one) from matrix

                Dic[new_a] = Local_Aligment_Protein(new_a, a)
    print(Dic)
    return Dic
    


def Thershold(ther, Dic):
    list=[]
    for key in Dic:
        if Dic[key] > ther:
            list.append(key)
            
    print(list)
    return list



def getString_Database(seed,sequence):
    list = []
    for line in open('sequence.fasta', 'r'):
        seq = line.strip()
        codon=[seq[i:i+3] for i in range(0,len(seq),3)]
    print(codon)
    for i in range(len(seed)):
        word=seed[i] in codon 
        if(word == True):
            list.append(seed[i])
    return(list)
    
    
   
def extend(trueSeed,HSP):
    for i in range(len(HSP)):
        if(HSP[i] == trueSeed):
            first = list(HSP[i-1])
            second = list(HSP[i+1])
            print(first[2]+trueSeed+second[0])
    
    


def Local_Aligment_Protein(first_seq_pro, second_seq_pro):
    len_first_seq = len(first_seq_pro)
    len_second_seq = len(second_seq_pro)
    max_value = 0
    maxcell = (0, 0)
    Aligment = [[0 for first_seq_pro in range(len_first_seq + 1)] for second_seq_pro in range(len_second_seq + 1)]
    for i in range(len_second_seq + 1):
        Aligment[i][0] = 0
    for j in range(len_first_seq + 1):
        Aligment[0][j] = 0
    for i in range(1, len_second_seq + 1, 1):
        for j in range(1, len_first_seq + 1, 1):

            First_aligment_index = l[0].index(first_seq_pro[j - 1])
            Second_aligment_index = l[0].index(second_seq_pro[i - 1])
            score = l[First_aligment_index][Second_aligment_index]
            Aligment[i][j] = max(Aligment[i - 1][j] - 10, Aligment[i][j - 1] - 10, Aligment[i - 1][j - 1] + int(score),
                                 0)
            if max_value <= Aligment[i][j]:
                max_value = Aligment[i][j]
                maxcell = (i, j)

    Gap_First_seq = ""
    Gap_Second_seq = ""

    i, j = maxcell
    while (Aligment[i][j] > 0):
        up = Aligment[i - 1][j] - 10
        left = Aligment[i][j - 1] - 10

        First_aligment_index = l[0].index(first_seq_pro[j - 1])
        Second_aligment_index = l[0].index(second_seq_pro[i - 1])
        score = l[First_aligment_index][Second_aligment_index]

        diagonal = Aligment[i - 1][j - 1] + int(score)
        if Aligment[i][j] == diagonal:
            Gap_First_seq += first_seq_pro[j - 1]
            Gap_Second_seq += second_seq_pro[i - 1]

            if first_seq_pro[j - 1] == second_seq_pro[i - 1]:
                i -= 1
                j -= 1
        elif Aligment[i][j] == up:

            Gap_Second_seq += second_seq_pro[i - 1]
            i -= 1
        else:
            Gap_First_seq += first_seq_pro[j - 1]

            j -= 1
        return max_value


choice = ""
Query = ""
word=[]
seq=[]
wordThershold = 0
wordLenght = 0
HSPThershold = 0
sequence=getseq()
choice = input("Enter to Blast")
if (choice == "Blast"):
    Query = input("Enter Query")
    wordLenght = input("wordLenght")
    wordLenght = int(wordLenght)
    wordThershold = input("wordThershold")
    wordThershold = int(wordThershold)
    I = 0
    I = int(I)
   
#    LADESRRLPTTTYLSEYYVAGEFPTADGPLPSRSPVLAANDGCDAEAGSGPPFFHVKK
    for I in range (len(Query)):
        if I+3	<= len(Query):
           word.append(Query[I:I+3])
    print(word)
   
Dic = Ditionary()
HSP=Thershold(wordThershold, Dic)
trueSeed = getString_Database(HSP,sequence)
extend(trueSeed,HSP)
