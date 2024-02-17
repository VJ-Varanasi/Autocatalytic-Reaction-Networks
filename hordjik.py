
from utils import *


def RAF_sim(n, fs, reps):
    if n == 2:
        t = 1
    else:
        t = 2

    
    output = []
    
    for f in fs:
        X,F,R = create_XFR(n, 0, 0)
        R_size = len(R)
        p= f/R_size
        #ps = np.linspace(0,p_max, reps)
        Xs = []
        Ys= []
        raf_count = 0
        for j in range(reps):
            X,F,R= create_XFR(n)
            p = f/len(R)
            C = C = create_catalysts(X, len(R),p)

            raf_count += RAF(X,F,R,C)
        #output.append([n, f, raf_count/reps])
        output.append(raf_count/reps)

    return output
            # Xs.append(p*R_size)
            # Ys.append(raf_count/Ns[i])
    
    #plt.plot(Xs,Ys, "-o", label= "n = {}".format(ns[i]))

    


reps = ast.literal_eval(sys.argv[1])
ns = ast.literal_eval(sys.argv[2])
f_range = ast.literal_eval(sys.argv[3])

fs = np.arange(f_range[0], f_range[1], f_range[2])
#ns = np.arange(n_range[0], n_range[1]+1, 1)





pool = multiprocess.Pool(processes=int(os.getenv('SLURM_CPUS_ON_NODE')))

if __name__ == '__main__':
    with pool as p:
        vals = []
        
        for i in range(len(ns)):
            vals = vals + [(ns[i], fs, reps[i])]
            
        results = p.starmap(RAF_sim, vals)

df = pd.DataFrame(results, index = ns, columns = fs)
#print(df)

plt.figure(figsize=(8,8))

for i in df.index:
    #print(df.loc[i])
    plt.plot(fs, df.loc[i], label = i)

plt.legend()
plt.xlabel("Expected Number of Catalysts")
plt.ylabel("Probability of RAF Set")
plt.title("RAF Probability Simulation")
plt.savefig("Images/PhaseTransition")
