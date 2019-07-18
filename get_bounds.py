import sys, math
import lbSumCol as lb

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("python get_bounds.py <nb vertex> <maximum IS size> <nb of maximum IS> <LB chromatic nbr>")
        print("exemple : python get_bounds.py 1000 15 6 ?")
        exit(0)

    n = int(sys.argv[1])
    maxIS = int(sys.argv[2])
    nbMaxIS = int(sys.argv[3])
    lb_chro_nb,  lb_chro_nb_is_int = lb.isInt(sys.argv[4])
    if not lb_chro_nb_is_int:
        lb_chro_nb = 0

    lb_chi0 = math.ceil(n/maxIS)
    ub   = lb.ub_chi(n, maxIS, nbMaxIS)
    em0  = lb.sigmaM0(n, maxIS, nbMaxIS)
    em0p = lb.sigmaM(n, maxIS, lb_chi0, nbMaxIS)
    em   = lb.sigmaM(n, maxIS, lb_chro_nb, nbMaxIS)
    lbme = lb.sigmaM(n, maxIS, lb_chro_nb, 100000000)

    print("ub chi =", ub)
    print("LBME   =", lbme)
    print("EM0    =", em0)
    print("EM0'   =", em0p)
    print("EM     =", em)
