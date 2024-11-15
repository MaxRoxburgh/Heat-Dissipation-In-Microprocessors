from Calculators import gradient_nat_convection, new_imaginary_points, grad_forced_convection, mirror_array, average_temp
from Initializing import initialise_boundary, np
from MatrixFunctions import update_T_matrix_method
from Plotting import trim


def update_with_source(T_k, T_0, q, h, T_a=293):
    # This new update function will use the fact that there will be a source term at each of the points
    # it will also update the new T_0 with the BC
   
    T_kP1 = T_0.copy()
   
    K_si = 150
   
    length = len(T_k)
    width = len(T_k[0].A1)
   
    for i in range(1, length - 1):
        for j in range(1, width - 1):
            ## i and j now index the entire inside of the array and does not touch any of the edges
           
            # the k_si and T_a which are global variables come from dimensionless and h = 1/stepsize(m)
            T_kP1[i,j] += 1/4*(T_k[i,j+1] + T_k[i,j-1] + T_k[i+1,j] + T_k[i-1,j]) + 1/4 * q *h **2 / (K_si*T_a)
               
    return T_kP1

def update_without_source(T_k, T_0):
    T_kP1 = T_0.copy()
   
    length = len(T_k)
    width = len(T_k[0].A1)
   

    for i in range(1, length - 1):
        for j in range(1, width - 1):
            ## i and j now index the entire inside of the array and does not touch any of the edges
           
            # no source term
            T_kP1[i,j] += 1/4*(T_k[i,j+1] + T_k[i,j-1] + T_k[i+1,j] + T_k[i-1,j])

    return T_kP1

def update_boundary_convection(T_k, K, T_a = 293, h=(1/6)*1e-3 ,nat=True):
    '''
    This func can be used for all of the different arrays of materials
    it calculates the update on the boundary due to convection
    the update for the boudary with different materials will be dealt with later and overwrite the relevant bdd values
    '''
    top = T_k[1].A1[1:-1]
    second_top = T_k[2].A1[1:-1]
    if nat:
        grad_top = np.array([gradient_nat_convection(i, T_a, K, h) for i in top])
    else:
        grad_top = np.array([grad_forced_convection(i, T_a, K) for i in top])
    new_top = new_imaginary_points(grad_top, second_top, h)
   
    col_list = T_k.transpose()
    left = col_list[1].A1[1:-1]
    second_left = col_list[2].A1[1:-1]
    if nat:
        grad_left = np.array([gradient_nat_convection(i, T_a, K, h) for i in left])
    else:
        grad_left = np.array([grad_forced_convection(i, T_a, K) for i in left])
    new_left = new_imaginary_points(grad_left, second_left, h)
   
    bottom = T_k[-2].A1[1:-1]
    second_bottom = T_k[-3].A1[1:-1]
    if nat:
        grad_bottom = np.array([gradient_nat_convection(i, T_a, K, h) for i in bottom])
    else:
        grad_bottom = np.array([grad_forced_convection(i, T_a, K) for i in bottom])
    new_bottom = new_imaginary_points(grad_bottom, second_bottom, h)
   
    right = col_list[-2].A1[1:-1]
    second_right = col_list[-3].A1[1:-1]
    if nat:
        grad_right = np.array([gradient_nat_convection(i, T_a, K, h) for i in right])
    else:
        grad_right = np.array([grad_forced_convection(i, T_a, K) for i in right])
    new_right = new_imaginary_points(grad_right, second_right,  h)
   
    T_0_new = initialise_boundary(new_top, new_left, new_right, new_bottom, T_a, [0,0,0,0], Get_T_K=False)
   
    return T_0_new

def calc_update_shared_bdd(T_below, T_above, k_below, k_above, T_0_below, T_0_above, T_a = 293, h=(1/6)*1e-3):
    # Pass in the second to last row from each and relative k values, this is the T_above/below variables
    # Finds the array of shared points
    # updates the two T_0 matricies by replacing the relevant values
   
    row_below, col_below = len(T_0_below), len(T_0_below[0].A1)
    row_above, col_above = len(T_0_above), len(T_0_above[0].A1)
    diff = col_above - col_below # it is assumed that the above is always bigger than the below
   
    above_middle_points = int(diff/2+1)
   
    # Takes relevant points and calculates the new imaginary point assuming gradients in and out must be equal
    new_points = (k_below*T_below[1:-1] + k_above*T_above[above_middle_points:-above_middle_points])/(k_below+k_above)
   
    # The new points are then used to rewrite the relevant parts of the T_0 matricies for the two materials
    T_below_boundary = np.matrix(np.append(np.append(T_0_below[0,0],new_points),T_0_below.A1[col_below-1:])).reshape(row_below,col_below)
   
    T_above_boundary = np.append(T_0_above.A1[:-int(col_above-1-diff/2)], np.append(new_points,T_0_above.A1[int((row_above*col_above)-1-diff/2):]))
    T_above_boundary = np.matrix(T_above_boundary).reshape(row_above,col_above)
   
    return T_below_boundary, T_above_boundary # These will be the new T_0's for the above and below case

def heat_sinc_boundary(T_below, T_above_list, T_0_below, T_0_above_list, spacing=20, fin_size=10):
   
    number_of_fins = len(T_above_list)
    T_above = [0]
    for i, item in enumerate(T_above_list): # makes the T_above with all of the boundary points
        if i+1 != number_of_fins:
            T_above = np.append(np.append(T_above, item[1:-1]),[0 for i in range(spacing)])
        else:
            T_above = np.append(np.append(T_above, item[1:-1]),[0])
       
    new_points = T_above
   
    row_below, col_below = len(T_0_below), len(T_0_below[0].A1)
   
    # sorting out the below boundarys
    T_below_boundary_top = T_0_below[0].A1
    for i in range(number_of_fins-1):
        T_below_boundary_top[1+i*(fin_size+spacing):(1+fin_size)+i*(fin_size+spacing)] = new_points[1+i*(fin_size+spacing):(1+fin_size)+i*(fin_size+spacing)]        
    T_below_boundary_top[col_below-(1+fin_size):col_below-1] = new_points[col_below-(1+fin_size):col_below-1]
    T_below_boundary = np.matrix(np.append(T_below_boundary_top,T_0_below.A1[col_below:])).reshape(row_below,col_below)
   
    new_points = T_below
   
    for index, T_0_above in enumerate(T_0_above_list):
        T_0_above[-1][1:-1] = np.array(new_points[1+index*(fin_size+spacing):(1+fin_size)+index*(fin_size+spacing)])
   
    return T_below_boundary, T_0_above_list

# def heat_sinc_boundary(T_below, T_above_list, T_0_below, T_0_above_list, spacing=20, fin_size=10):
   
#     number_of_fins = len(T_above_list)
#     T_above = [0]
#     for i, item in enumerate(T_above_list): # makes the T_above with all of the boundary points
#         if i+1 != number_of_fins:
#             T_above = np.append(np.append(T_above, item[1:-1]),[0 for i in range(spacing)])
#         else:
#             T_above = np.append(np.append(T_above, item[1:-1]),[0])
       
#     new_points = (T_below + T_above)/(2)
   
#     row_below, col_below = len(T_0_below), len(T_0_below[0].A1)
   
#     # sorting out the below boundarys
#     T_below_boundary_top = T_0_below[0].A1
#     for i in range(number_of_fins-1): # this is not general for different dimensions, but is for more fins
#         T_below_boundary_top[1+i*(fin_size+spacing):(1+fin_size)+i*(fin_size+spacing)] = new_points[1+i*(fin_size+spacing):(1+fin_size)+i*(fin_size+spacing)]        
#     T_below_boundary_top[col_below-(1+fin_size):col_below-1] = new_points[col_below-(1+fin_size):col_below-1]
#     T_below_boundary = np.matrix(np.append(T_below_boundary_top,T_0_below.A1[col_below:])).reshape(row_below,col_below)
   
   
#     for index, T_0_above in enumerate(T_0_above_list):
#         T_0_above[-1][1:-1] = np.array(new_points[1+index*(fin_size+spacing):(1+fin_size)+index*(fin_size+spacing)])
   
#     return T_below_boundary, T_0_above_list

def full_matrix_iteration_method(T_0_pr_update, T_0_cs_update, T_0_hs_update, T_0_fn_update, T_k_pr_update
                                 , T_k_cs_update, T_k_hs_update, T_k_fn_update, h, M_inv_pr, M_inv_cs
                                 , M_inv_hs, M_inv_fn, vol_pr, vol_cs, vol_sinc, vol_fins, path, title, iterations=150000, scaling=3, nat=True):
   
    from MatrixFunctions import update_T_matrix_method    
    from Calculators import average_temp
    import pandas as pd
    import Plotting as plot
   
    q = 0.5*1e9 # W/m^3
    K_si = 150 # W/mK
    K_ceramic = 230 # W/mK
    K_alu = 250 # W/mK
    n_fins = len(T_0_fn_update)
   
   
    avT_k_1 = max(T_0_pr_update[0].A1)
    avT_k_1_2 = max(T_0_pr_update[0].A1)
    diff_array_big = []
    av_array_big = []
    av_pr_array = []
    av_cs_array = []
    av_hs_array = []
    every = iterations//15
    T_sys = []
   
    for i in range(iterations):
       
        # update boundaries convection
        T_0_pr_update, T_0_cs_update = update_boundary_convection(T_k_pr_update, K_si, nat=nat), update_boundary_convection(T_k_cs_update, K_ceramic, nat=nat)
        T_0_hs_update = update_boundary_convection(T_k_hs_update, K_alu, nat=nat)
        T_0_fn_update = np.array([update_boundary_convection(np.matrix(fin_k), K_alu, nat=nat) for fin_k in T_k_fn_update])
   
       
        # update boudaries shared
        T_0_pr_update, T_0_cs_update = calc_update_shared_bdd(T_k_pr_update[2].A1, T_k_cs_update[-3].A1, K_si, K_ceramic, T_0_pr_update, T_0_cs_update)
        T_0_cs_update, T_0_hs_update = calc_update_shared_bdd(T_k_cs_update[2].A1, T_k_hs_update[-3].A1, K_ceramic, K_alu, T_0_cs_update, T_0_hs_update)
        T_above_list = np.array([fin_k[-2] for fin_k in mirror_array(T_k_fn_update, True)])
        T_below = T_k_hs_update[1].A1
        T_0_fn_update_2 = mirror_array(T_0_fn_update, True)
        T_0_hs_update, T_0_fn_update = heat_sinc_boundary(T_below, T_above_list, T_0_hs_update, T_0_fn_update_2, spacing=2*scaling, fin_size=1*scaling)
        T_0_fn_update = np.array([T_0_fn_update[i] for i in range(n_fins)])
       
       
        # update inside temps
        T_k_pr_update = update_T_matrix_method(T_0_pr_update, M_inv_pr, q, h)
        T_k_cs_update = update_T_matrix_method(T_0_cs_update, M_inv_cs)
        T_k_hs_update = update_T_matrix_method(T_0_hs_update, M_inv_hs)
        T_k_fn_update = np.array([update_T_matrix_method(np.matrix(fin_0), M_inv_fn) for fin_0 in T_0_fn_update])
       
        if i%every == 0:
            average_T_pr = average_temp(T_k_pr_update)
            average_T_cs = average_temp(T_k_cs_update)
            average_T_fn = sum([average_temp(np.matrix(T_k)) for T_k in T_k_fn_update])/len(T_k_fn_update)
            average_T_hs = average_temp(T_k_cs_update)
           
            average_fin_hs = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs)/(2*vol_fins+vol_sinc)
            print('\n\n\n-------------------------------------------------\n\n\n')
            print('T:     pr     |     cs     |      hs      |\n', average_T_pr, average_T_cs, average_fin_hs)
            ambient = min(T_k_fn_update[0][0])
            plot.plot_whole_system(T_k_cs_update, T_k_pr_update, mirror_array(T_k_fn_update, True), T_k_hs_update, ambient, 6)    
           
            avT_k = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs + vol_cs*average_T_cs + vol_pr* average_T_pr)/(2*vol_fins+vol_sinc+vol_cs+vol_pr)
            diff = avT_k-avT_k_1
            print(f'\nAverage T, difference from {every} runs ago: {i}\n',avT_k, diff)
            avT_k_1 = avT_k
               
           
        if i%int(every//10) == 0 and i>0:
           
           
            average_T_pr = average_temp(T_k_pr_update)
            average_T_cs = average_temp(T_k_cs_update)
            average_T_fn = sum([average_temp(np.matrix(T_k)) for T_k in T_k_fn_update])/len(T_k_fn_update)
            average_T_hs = average_temp(T_k_cs_update)
            average_fin_hs = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs)/(2*vol_fins+vol_sinc)
            avT_k_2 = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs + vol_cs*average_T_cs + vol_pr* average_T_pr)/(2*vol_fins+vol_sinc+vol_cs+vol_pr)
            diff_2 = avT_k_2-avT_k_1_2
            diff_array_big.append(diff_2)
            av_array_big.append(avT_k_2)
            av_pr_array.append(average_T_pr)
            av_cs_array.append(average_T_cs)
            av_hs_array.append(average_fin_hs)
            avT_k_1_2 = avT_k_2
           
        if i%5000 == 0 and i>0:
            average_T_pr = average_temp(T_k_pr_update)
            average_T_cs = average_temp(T_k_cs_update)
            average_T_fn = sum([average_temp(np.matrix(T_k)) for T_k in T_k_fn_update])/len(T_k_fn_update)
            average_T_hs = average_temp(T_k_cs_update)
            T_sys.append((2*vol_fins*average_T_fn + vol_sinc*average_T_hs + vol_cs*average_T_cs + vol_pr* average_T_pr)/(2*vol_fins+vol_sinc+vol_cs+vol_pr))
            # print('\n', T_sys, '\n')
           
            if len(T_sys)>2:
                check = np.array(T_sys)[-2:]
                if abs(check[1]-check[0]) < 1e-3/(n_fins**2):
                   
                    print(f'\n\n     ::::   Converged   ::::\n    After {i} iterations\n')
                    ambient = min(T_k_fn_update[0][0])
                    plot.plot_whole_system(T_k_cs_update, T_k_pr_update, mirror_array(T_k_fn_update, True), T_k_hs_update, ambient, 6)
                    plot.save(T_k_pr_update, T_k_cs_update, T_k_hs_update, T_k_fn_update, title, f"{path}")
                    df = pd.DataFrame(np.matrix([diff_array_big, av_array_big, av_pr_array, av_cs_array, av_hs_array]).transpose())
                    df.columns = ["diff_array_big", "av_array_big", "av_pr_array", "av_cs_array", "av_hs_array"]
                    df.to_csv(f'{path}/data_{title}', index=False)
                    print('\n\n\n-------------------------------------------------\n')
                    print('-----------------------fin-----------------------')
                    print('\n-------------------------------------------------\n\n\n')
                    return
   
   
    plot.save(T_k_pr_update, T_k_cs_update, T_k_hs_update, T_k_fn_update, title, f"{path}")
           
           
           
           
    df = pd.DataFrame(np.matrix([diff_array_big, av_array_big, av_pr_array, av_cs_array, av_hs_array]).transpose())
    df.columns = ["diff_array_big", "av_array_big", "av_pr_array", "av_cs_array", "av_hs_array"]
    df.to_csv(f'{path}/data_{title}', index=False)
    print('\n\n\n-------------------------------------------------\n')
    print('-----------------------fin-----------------------')
    print('\n-------------------------------------------------\n\n\n')
   
   
    plot.save(T_k_pr_update, T_k_cs_update, T_k_hs_update, T_k_fn_update, title, f"{path}")
   
# def update_system_n_times(n, scaling, T_k_pr, T_k_cs, T_k_hs, T_k_fn, T_0_pr, T_0_cs, T_0_hs, T_0_fn, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat=True):
#     q = 0.5*1e9 # W/m^3
#     K_si = 150 # W/mK
#     K_ceramic = 230 # W/mK
#     K_alu = 250 # W/mK
#     n_fins = len(T_0_fn)
#     h = (1/scaling)*1e-3
   
    
#     for i in range(n):
#         # update boundaries convection
#         T_0_pr, T_0_cs = update_boundary_convection(T_k_pr, K_si, nat=nat, h=h), update_boundary_convection(T_k_cs, K_ceramic, nat=nat, h=h)
#         T_0_hs = update_boundary_convection(T_k_hs, K_alu, nat=nat, h=h)
#         T_0_fn = np.array([update_boundary_convection(np.matrix(fin_k), K_alu, nat=nat, h=h) for fin_k in T_k_fn])
   
       
#         # update boudaries shared
#         T_0_pr, T_0_cs = calc_update_shared_bdd(T_k_pr[2].A1, T_k_cs[-3].A1, K_si, K_ceramic, T_0_pr, T_0_cs, h=h)
#         T_0_cs, T_0_hs = calc_update_shared_bdd(T_k_cs[2].A1, T_k_hs[-3].A1, K_ceramic, K_alu, T_0_cs, T_0_hs, h=h)
#         T_above_list = np.array([fin_k[-3] for fin_k in mirror_array(T_k_fn, True)])
#         T_below = T_k_hs[-3].A1
#         T_0_fn_2 = mirror_array(T_0_fn, True)
#         T_0_hs, T_0_fn = heat_sinc_boundary(T_below, T_above_list, T_0_hs, T_0_fn_2, spacing=2*scaling, fin_size=1*scaling)
#         T_0_fn = np.array([T_0_fn[i] for i in range(n_fins)])
       
       
#         # update inside temps
#         T_k_pr = update_T_matrix_method(T_0_pr, M_inv_pr, q, h)
#         T_k_cs = update_T_matrix_method(T_0_cs, M_inv_cs)
#         T_k_hs = update_T_matrix_method(T_0_hs, M_inv_hs)
#         T_k_fn = np.array([update_T_matrix_method(np.matrix(fin_0), M_inv_fn) for fin_0 in T_0_fn])
       
#     average_T_pr, vol_pr = average_temp(T_k_pr, True)
#     average_T_cs, vol_cs = average_temp(T_k_cs, True)
#     average_T_fn = sum([average_temp(np.matrix(T_k)) for T_k in T_k_fn])/len(T_k_fn)
#     average_T_hs, vol_sinc= average_temp(T_k_cs, True)
#     vol_fins = len(np.matrix(T_0_fn[0]).A1)*len(T_0_fn)*2
   
#     av_T = (average_T_pr*vol_pr + average_T_cs*vol_cs + average_T_hs*vol_sinc +average_T_fn*vol_fins)/(vol_pr + vol_cs + vol_sinc + vol_fins)
   
#     return T_k_pr, T_k_cs, T_k_hs, T_k_fn, T_0_pr, T_0_cs, T_0_hs, T_0_fn, av_T

def update_system_n_times(n, scaling, T_k_pr, T_k_cs, T_k_hs, T_k_fn, T_0_pr, T_0_cs, T_0_hs, T_0_fn, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat=True, fin_length = 30,spacing=2):
    q = 0.5*1e9 # W/m^3
    K_si = 150 # W/mK
    K_ceramic = 230 # W/mK
    K_alu = 250 # W/mK
    n_fins = len(T_0_fn)
    h = (1/scaling)*1e-3
   
    
    for i in range(n):
        # update boundaries convection
        T_0_pr, T_0_cs = update_boundary_convection(T_k_pr, K_si, nat=nat, h=h), update_boundary_convection(T_k_cs, K_ceramic, nat=nat, h=h)
        T_0_hs = update_boundary_convection(T_k_hs, K_alu, nat=nat, h=h)
        T_0_fn = np.array([update_boundary_convection(np.matrix(fin_k), K_alu, nat=nat, h=h) for fin_k in T_k_fn])
   
       
        # update boudaries shared
        T_0_pr, T_0_cs = calc_update_shared_bdd(T_k_pr[2].A1, T_k_cs[-3].A1, K_si, K_ceramic, T_0_pr, T_0_cs, h=h)
        T_0_cs, T_0_hs = calc_update_shared_bdd(T_k_cs[2].A1, T_k_hs[-3].A1, K_ceramic, K_alu, T_0_cs, T_0_hs, h=h)
        T_above_list = np.array([fin_k[-2] for fin_k in mirror_array(T_k_fn, True)])
        T_below = T_k_hs[1].A1
        T_0_fn_2 = mirror_array(T_0_fn, True)
        T_0_hs, T_0_fn = heat_sinc_boundary(T_below, T_above_list, T_0_hs, T_0_fn_2, spacing=spacing*scaling, fin_size=1*scaling)
        T_0_fn = np.array([T_0_fn[i] for i in range(n_fins)])
       
       
        # update inside temps
        T_k_pr = update_T_matrix_method(T_0_pr, M_inv_pr, q, h)
        T_k_cs = update_T_matrix_method(T_0_cs, M_inv_cs)
        T_k_hs = update_T_matrix_method(T_0_hs, M_inv_hs)
        T_k_fn = np.array([update_T_matrix_method(np.matrix(fin_0), M_inv_fn) for fin_0 in T_0_fn])
       
    average_T_pr, vol_pr = average_temp(T_k_pr, True)
    average_T_cs, vol_cs = average_temp(T_k_cs, True)
    average_T_fn = sum([average_temp(np.matrix(T_k)) for T_k in T_k_fn])/len(T_k_fn)
    average_T_hs, vol_sinc= average_temp(T_k_hs, True)
    vol_fins = len(trim(np.matrix(T_0_fn[0])).A1)*len(T_0_fn)*2
   
    av_T = (average_T_pr*vol_pr + average_T_cs*vol_cs + average_T_hs*vol_sinc +average_T_fn*vol_fins)/(vol_pr + vol_cs + vol_sinc + vol_fins)
   
    return T_k_pr, T_k_cs, T_k_hs, T_k_fn, T_0_pr, T_0_cs, T_0_hs, T_0_fn, av_T


def update_pr_cs_n_times(n, scaling, T_k_pr, T_k_cs, T_0_pr, T_0_cs, M_inv_pr, M_inv_cs, nat=True):
    q = 0.5*1e9 # W/m^3
    K_si = 150 # W/mK
    K_ceramic = 230 # W/mK
    # K_alu = 250 # W/mK
    h = (1/scaling)*1e-3

   
    for i in range(n):
        # update boundaries convection
        T_0_pr, T_0_cs = update_boundary_convection(T_k_pr, K_si, nat=nat, h=h), update_boundary_convection(T_k_cs, K_ceramic, nat=nat, h=h)  
       
        # update boudaries shared
        T_0_pr, T_0_cs = calc_update_shared_bdd(T_k_pr[2].A1, T_k_cs[-3].A1, K_si, K_ceramic, T_0_pr, T_0_cs, h=h)      
       
        # update inside temps
        T_k_pr = update_T_matrix_method(T_0_pr, M_inv_pr, q, h)
        T_k_cs = update_T_matrix_method(T_0_cs, M_inv_cs)
       
    average_T_pr, vol_pr = average_temp(T_k_pr, True)
    average_T_cs, vol_cs = average_temp(T_k_cs, True)
   
    av_T =(average_T_cs*vol_cs + average_T_pr*vol_pr)/(vol_cs+vol_pr)
    return T_k_pr, T_k_cs, T_0_pr, T_0_cs, av_T
