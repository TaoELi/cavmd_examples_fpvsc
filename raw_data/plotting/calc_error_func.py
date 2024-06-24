import numpy as np

cavity_1 = {"ngrid" : 36, "lx_max" : 36, "centered" : False}
cavity_2 = {"ngrid" : 36, "lx_max" : 108, "centered" : False}
cavity_3 = {"ngrid" : 6, "lx_max" : 36, "centered" : False}
cavity_4 = {"ngrid" : 36, "lx_max" : 36, "centered" : True}


def calc_error_func(cavity):
    ngrid = cavity["ngrid"]
    xgrid = np.linspace(1, ngrid, ngrid) / (ngrid + 1)
    if cavity["centered"]:
        xgrid = np.array([0.5] * ngrid)
    print("Ngrid = ", ngrid)
    print("molecular distribution is ", xgrid)
    lx_max = cavity["lx_max"]
    print("lx_max = ", lx_max)

    # Calculate the error function
    # Error(kx) = 1/ngrid sum_lx sum_kx' sin(kx'*xi) * sin(kx*xi)
    err_lst = []
    for i in range(1, lx_max+1):
        err = 0.0
        for j in range(1, lx_max+1):
            if j != i:
                # calculate the error function
                kx = np.pi * i
                kx_prime = np.pi * j
                err += np.sum(np.sin(xgrid * kx_prime) * np.sin(xgrid * kx))
        err_lst.append(err / ngrid)
    err_lst = np.array(err_lst)

    print("error functions are ", err_lst)
    print("avg of error functions is ", np.mean(err_lst))
    print("absolute avg of error functions is ", np.mean(np.abs(err_lst)))
    print("non-zero error functions are ", err_lst[np.abs(err_lst) > 1e-5])
    print("non-zero error function percentage ", np.size(err_lst[np.abs(err_lst) > 1e-5]) / lx_max * 100)
    print("\n\n\n")


calc_error_func(cavity_1)
calc_error_func(cavity_2)
calc_error_func(cavity_3)
calc_error_func(cavity_4)