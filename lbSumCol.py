import sys, time, subprocess
import os, random, math
from subprocess import STDOUT, check_output

PATH_EXE_MOMC = "./MoMC"
TMP_DIR = "tmp"
RESULT_DIR = "results"
TIME_LIMIT = 600
MAX_SIZE_GRAPH = 3000

def get_info_from_MoMCoutput(filename):
    iss = []
    maxCLQ, nbMaxCLQ, time = None, None, -1.
    with open(filename, 'r') as f:
        for line in f:
            words = line.strip().split()
            if words[0] == 'M':
                iss.append([int(x) for x in words[1:]])
            if words[0] == 's':
                maxCLQ = int(words[4])
                if len(words) == 15:
                    nbMaxCLQ = int(words[6])
                    time = float(words[10])
                else:
                    time = float(words[8])
    return maxCLQ, nbMaxCLQ, time, iss

def build_graph_IS_max(iss):
    n = len(iss)
    #print(n)
    graph = [[] for _ in range(n)]
    for i,k in enumerate(iss):
        d = {x:True for x in k}
        for j in range(i+1,n):
            l = iss[j]
            edge = False
            for x in l:
                if x in d:
                    edge = True
                    break
            if edge:
                graph[i].append(j)
                graph[j].append(i)
    return graph

def build_complement_graph_IS_max(iss):
    n = len(iss)
    #print(n)
    graph = [[] for _ in range(n)]
    for i,k in enumerate(iss):
        d = {x:True for x in k}
        for j in range(i+1,n):
            l = iss[j]
            edge = False
            for x in l:
                if x in d:
                    edge = True
                    break
            if not edge:
                graph[i].append(j)
                graph[j].append(i)
    return graph

def readDIMACS(filename):
    g = []
    with open(filename, 'r') as f:
        for line in f:
            words = line.strip().split()
            if len(words) > 0 :
                if words[0] == 'p':
                    n = int(words[2])
                    m = int(words[3])
                    g = [[] for _ in range(n)]
                if words[0] == 'e':
                    u = int(words[1])-1
                    v = int(words[2])-1
                    if v not in g[u]:
                        g[u].append(v)
                    if u not in g[v]:
                        g[v].append(u)
    return g

def writeDIMACS(g, filename, textComent=''):
    with open(filename, 'w') as f:
        n = len(g)
        m = 0
        text = ""
        for v, lv in enumerate(g):
            for u in lv:
                if v < u:
                    m += 1
                    text += "e {} {}\n".format(v+1, u+1)
        for tc in textComent : f.write("c "+tc+"\n")
        if n != 0 :
            f.write("c density {:.3f}\n".format(m/n/(n+1)*2))
        f.write("p edge {} {}\n".format(n, m))
        f.write(text)
    print("Graph file {} is created.".format(filename))

def complement(graph):
    n = len(graph)
    return [[i for i in range(n) if i != v and i not in l] for v,l in enumerate(graph)]

def read_write_complement_graph(filename_in):
    graphname = os.path.basename(filename_in)
    #path = os.path.dirname(filename_in)
    if not os.path.exists(TMP_DIR):
        print("Diretory {} is created.".format(TMP_DIR))
        os.makedirs(TMP_DIR)
    filename_out = "{}/{}_c.col".format(TMP_DIR, graphname[:-4])
    g = readDIMACS(filename_in)
    gc = complement(g)
    writeDIMACS(gc, filename_out, textComent=['Complement graph of '+ graphname])
    return filename_out

def get_outputname(nameinstance, only_one=True):
    graphname = os.path.basename(nameinstance)
    if not os.path.exists(RESULT_DIR):
        os.makedirs(RESULT_DIR)
        print("Diretory {} is created.".format(RESULT_DIR))
    if only_one:
        outputfile = "{}/{}_one.txt".format(RESULT_DIR, graphname[:-4])
    else:
        outputfile = "{}/{}_all.txt".format(RESULT_DIR, graphname[:-4])
    return outputfile

def run_MoMC(nameinstance, time_limit, only_one=True):
    fileoutput = get_outputname(nameinstance, only_one)
    if only_one:
        exe = [PATH_EXE_MOMC, nameinstance, ">", fileoutput]
    else:
        #nb_max_clq = 10000
        exe = [PATH_EXE_MOMC, nameinstance, "-a",  "{}".format(MAX_SIZE_GRAPH), ">", fileoutput]
    cmd = ' '.join(exe)
    #print(cmd)
    #alpha, dt = None, -1
    maxCLQ, nbMaxCLQ, dtime, iss = None, None, -1., []
    t0 = time.time()
    #print('t0=',t0)
    proc = subprocess.Popen(
        cmd,
        stderr=subprocess.STDOUT,  # Merge stdout and stderr
        stdout=subprocess.PIPE,
        shell=True)
    try:
        print("Result file {} .....".format(fileoutput), end="", flush=True)
        outs, errs = proc.communicate(timeout=time_limit)
        print(" is created.")
        #print("outs:",outs)
        #print("errs:",errs)
        dt = time.time() - t0
        maxCLQ, nbMaxCLQ, dtime, iss = get_info_from_MoMCoutput(fileoutput)
        #print(dtime, dt)
        #print(dt, alpha)
        #key = outs.strip()
    except subprocess.TimeoutExpired:
        dtime = time_limit
        print('Time limit exceeded.')
    proc.kill()
    #print('Time limit exceeded.')
    return fileoutput, maxCLQ, nbMaxCLQ, dtime, iss

def ub_alpha(g):
    n = len(g)
    d = [n-1-len(v) for v in g]
    d.sort()
    d.reverse()
    alpha = n
    for k,dk in enumerate(d):
        if k+1 <= dk:
            alpha = k+1
            #print(k,dk,alpha)
        else:
            break
    return alpha

def sigmaM00(n, alpha):
    return sigmaM0(n, alpha, int(n/alpha))

def sigmaM0(n, alpha, m):
    return sigmaM(n, alpha, 0, m)

def sigmaM(n, alpha, s, m):
    if m*alpha > n:
        m = int(n/alpha)
    q = (n-m*alpha) // (alpha-1)
    r = (n-m*alpha)-q*(alpha-1)
    nbcol =  m + q + 1 if r > 0 else m + q
    if nbcol >= s:
        return int(m*(m+1)/2*alpha + q*(2*m+q+1)/2*(alpha-1) + (m+q+1)*r)
    else:
        return int(s*(s+1)/2) + sigmaM0(n-s, alpha-1, m)

def ub_chi(n, alpha, m):
    if m*alpha > n:
        m = int(n/alpha)
    q = (n-m*alpha) // (alpha-1)
    r = (n-m*alpha)-q*(alpha-1)
    return m + q + 1 if r > 0 else m + q

def isInt(x):
    x_int, x_is_int = None, False
    try:
        x_int = int(x)
        x_is_int = True
    except ValueError:
        pass
    return x_int, x_is_int

def printtime(dt1):
    return "{:.0f}".format(dt1) if dt1 >= 1 else ("{:.1f}".format(dt1) if dt1 > 0.1 else "$<0.1$")

def printbest(x, best):
    return "\\bf{}".format(x) if x == best else "{}".format(x)

def outlatex(graphname, n, density, alpha1, dt1, nb, dt2, alpha_G, dt3,
             nbChro, ubChi, lbme, em0, em, chro_sum_int, chro_sum_is_int,
             best_ocs_int, best_ocs_is_int):
    #dd1 = "{:.0f}".format(dt1) if dt1 >= 1 else ("{:.1f}".format(dt1) if dt1 > 0.1 else "<0.1")
    best = max([lbme, em0, em, chro_sum_int if chro_sum_is_int else 0, best_ocs_int if best_ocs_is_int else 0])
    text = "{} & {} & {:.2f} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {}\\\ ".format(graphname[:-4], n, density, alpha1, printtime(dt1), nb, printtime(dt2), alpha_G, 0 if alpha_G == 1 else printtime(dt3), nbChro, ubChi, printbest(chro_sum_int,best) if chro_sum_is_int else '?', printbest(best_ocs_int,best) if best_ocs_is_int else '?', printbest(lbme,best), printbest(em0,best), printbest(em,best))
    text = text.replace('_', '\_')
    print(text)

def main():
    if not (len(sys.argv) == 4 or len(sys.argv) == 6):
        print("python lbSumCol.py <graph instance> <max seconds for MoMC> <lower bound of chromatic number>\nor")
        print("python lbSumCol.py <graph instance> <max seconds for MoMC> <lower bound of chromatic number> <chromatic sum if known> <best old chromatic sum>")
        print("exemples : python lbSumCol.py dimacs/DSJC1000.9.col 100 216")
        print("           python lbSumCol.py dimacs/DSJC1000.9.col 100 0 ? ?")
        exit(0)
    nameinstance = sys.argv[1]
    maxseconds = int(sys.argv[2])
    chro_nb = int(sys.argv[3])
    chro_sum, best_old_chro_sum = "?", "?"
    if len(sys.argv) == 6:
        chro_sum = sys.argv[4]
        best_old_chro_sum = sys.argv[5]
    chro_sum_int, chro_sum_is_int = isInt(chro_sum)
    best_ocs_int, best_ocs_is_int = isInt(best_old_chro_sum)
    graphname = os.path.basename(nameinstance)
    print("graph name             = {}".format(graphname))
    #print("graphpath   = {}".format(nameinstance))
    print("max seconds for MoMC   = {}".format(maxseconds))
    print("LB of chromatic number = {}".format(chro_nb))
    #print("filename_out  = {}".format(filename_out))
    print("chromatic sum          = {}".format(chro_sum_int if chro_sum_is_int else "unknown"))
    print("best old chromatic sum = {}".format(best_ocs_int if  best_ocs_is_int else "unknown"))
    #t0 = time.time()
    print("Reading graph..... ", end='', flush=True)
    g = readDIMACS(nameinstance)
    print("OK", flush=True)
    n = len(g)
    m = sum([len(x) for x in g]) / 2
    density = m / (n*(n-1)/2)
    #print(graphname, n, m, density)
    apha_ub = None
    maxIS, nbMaxIS, dtime2, nbMaxIS2, dtime3 = None, None, -1., None, -1.
    print("Creating complement graph..... ", end='', flush=True)
    nameinstance_c = read_write_complement_graph(nameinstance)
    _, maxIS, _, dtime, _ = run_MoMC(nameinstance_c, maxseconds, only_one=True)
    #print(maxIS, dtime)
    print("IS max     = {} in {}s".format(maxIS, dtime), flush=True)
    if dtime < maxseconds:
         outputname, _, _, _, _ = run_MoMC(nameinstance_c, maxseconds, only_one=False)
         _, nbMaxIS, dtime2, iss = get_info_from_MoMCoutput(outputname)
         print("nb IS max  = {} in {}s".format(nbMaxIS, dtime2))
         if nbMaxIS is not None and 1 < nbMaxIS < MAX_SIZE_GRAPH:
             filename_graph_iss_c = "{}/{}_is_c.col".format(TMP_DIR, graphname[:-4])
             if not os.path.exists(filename_graph_iss_c ):
                 print("Creating complement graph of max IS ..... ", end="", file=sys.stdout, flush=True)
                 dt_gc0 = time.time()
                 graph_iss_c = build_complement_graph_IS_max(iss)
                 dt_gc = time.time() - dt_gc0
                 print("is done ({:.2f}s)".format(dt_gc), file=sys.stdout, flush=True)
                 writeDIMACS(graph_iss_c, filename_graph_iss_c, textComent=['Complement graph of maximum IS of grpah '+ graphname])
             outputname2, nbMaxIS2, _, dtime3, _ = run_MoMC(filename_graph_iss_c, maxseconds, only_one=True)
             if nbMaxIS2 is None:
                 nbMaxIS2 = nbMaxIS
         else:
             if nbMaxIS is None:
                 nbMaxIS = int(n/maxIS)
             nbMaxIS2 = nbMaxIS
         print("nb IS max* = {} in {}s".format(nbMaxIS2, dtime3))
         ubChi = ub_chi(n, maxIS, nbMaxIS2)
         lbme0= sigmaM00(n, maxIS)
         s = math.ceil(n/maxIS)
         lbme0p= sigmaM(n, maxIS, s, 100000000)
         lbme = sigmaM(n, maxIS, chro_nb, 100000000)
         em0  = sigmaM0(n, maxIS, nbMaxIS2)
         em0p = sigmaM(n, maxIS, s, nbMaxIS2)
         em   = sigmaM(n, maxIS, chro_nb, nbMaxIS2)

         print("UB chi =", ubChi)
         #print("LBME0  =", lbme0)
         #print("LBME0' =", lbme0p)
         print("LBME   =", lbme)
         print("EM0    =", em0)
         print("EM0'   =", em0p)
         print("EM     =", em)
         #output for our article
         if len(sys.argv) == 6:
             outlatex(graphname, n, density, maxIS, dtime, nbMaxIS, dtime2,
                      nbMaxIS2, dtime3, chro_nb, ubChi, lbme, em0, em,
                      chro_sum_int, chro_sum_is_int, best_ocs_int,
                      best_ocs_is_int)

if __name__ == "__main__":
    main()
