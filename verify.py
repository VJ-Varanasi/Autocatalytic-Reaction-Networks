from utils import *

success = glob.glob("Data/Dict-Matrix-*[0-9]-RAF.txt")

print("SUCCESS")
for file in success:
    
    n = int(file.split("-")[2])
    f = float(file.split("-")[3])

    file1 = open(file, 'r')
    Lines = file1.readlines()

    Cs = [ast.literal_eval(line) for line in Lines]
    X,F,R = create_XFR(n)
    
    incorrect = 0
    
    for C in Cs:
        if RAF(X,F,R,C) !=1:
            incorrect +=1
        
    print("({},{}): {} across {}".format(n,f, incorrect/len(Cs),len(Cs) ))


print("FAILURE")
failure = glob.glob("Data/Dict-Matrix-*-Non-RAF.txt")

for file in failure:
    
    n = int(file.split("-")[2])
    f = float(file.split("-")[3])

    file1 = open(file, 'r')
    Lines = file1.readlines()

    Cs = [ast.literal_eval(line) for line in Lines]
    X,F,R = create_XFR(n)

    print(n,f)
    if len(Cs)!=0:
        incorrect = 0
        
        for C in Cs:
            if RAF(X,F,R,C) !=0:
                incorrect +=1
            
        print("({},{}): {} across {}".format(n,f, incorrect/len(Cs),len(Cs) ))



