#Copyright 2017 Hang Yin
#Transformation between stress and strain.
import numpy as np
import math
import matplotlib.pyplot as plt
v_judge=1
while v_judge:
    v = float(input('Please input Poisson ratio:(-1<v<0.5):'))
    if v<=-1 or v>=0.5:
        print('Retry please.')
    else:
        v_judge=0
E=float(input('Young\'s modulus is:' ))
def solve_for_given_stress(dir_vec,a):
    length=math.sqrt(dir_vec[0][0]*dir_vec[0][0]+dir_vec[0][1]*dir_vec[0][1]+dir_vec[0][2]*dir_vec[0][2])
    for i in range(0,3):
        dir_vec[0][i]=dir_vec[0][i]/length
    stress_x=a[0][0]*dir_vec[0][0]+a[0][1]*dir_vec[0][1]+a[0][2]*dir_vec[0][2]
    stress_y = a[1][0] * dir_vec[0][0] + a[1][1] * dir_vec[0][1] + a[1][2] * dir_vec[0][2]
    stress_z = a[2][0] * dir_vec[0][0] + a[2][1] * dir_vec[0][1] + a[2][2] * dir_vec[0][2]
    stress_square=stress_x*stress_x+stress_y*stress_y+stress_z*stress_z
    normal_stress=a[0][0]*dir_vec[0][0]*dir_vec[0][0]+a[1][1]*dir_vec[0][1]*dir_vec[0][1]+a[2][2]*dir_vec[0][2]*dir_vec[0][2]+\
        2*a[0][1]*dir_vec[0][0]*dir_vec[0][1]+2*a[0][2]*dir_vec[0][0]*dir_vec[0][2]+2*a[2][1]*dir_vec[0][2]*dir_vec[0][1]
    shear_stress=math.sqrt(stress_square-math.pow(normal_stress,2))
    print('The normal stress is {} Pa.'.format(normal_stress))
    print('The shear stress is {} Pa.'.format(shear_stress))
    return 0
def solve_for_given_strain(dir_vec,a):
    length=math.sqrt(dir_vec[0][0]*dir_vec[0][0]+dir_vec[0][1]*dir_vec[0][1]+dir_vec[0][2]*dir_vec[0][2])
    for i in range(0,3):
        dir_vec[0][i]=dir_vec[0][i]/length
    strain_x=a[0][0]*dir_vec[0][0]+a[0][1]*dir_vec[0][1]+a[0][2]*dir_vec[0][2]
    strain_y = a[1][0] * dir_vec[0][0] + a[1][1] * dir_vec[0][1] + a[1][2] * dir_vec[0][2]
    strain_z = a[2][0] * dir_vec[0][0] + a[2][1] * dir_vec[0][1] + a[2][2] * dir_vec[0][2]
    strain_square=strain_x*strain_x+strain_y*strain_y+strain_z*strain_z
    normal_strain=a[0][0]*dir_vec[0][0]*dir_vec[0][0]+a[1][1]*dir_vec[0][1]*dir_vec[0][1]+a[2][2]*dir_vec[0][2]*dir_vec[0][2]+\
        2*a[0][1]*dir_vec[0][0]*dir_vec[0][1]+2*a[0][2]*dir_vec[0][0]*dir_vec[0][2]+2*a[2][1]*dir_vec[0][2]*dir_vec[0][1]
    shear__strain=math.sqrt(strain_square-math.pow(normal_strain,2))
    shear_strain=2*shear__strain
    print('The normal strain is {}.'.format(normal_strain))
    print('The shear strain is {}.'.format(shear_strain))
    return 0
def principal_stress(a):
    b, c = np.linalg.eig(a)
    d = [[0 for i in range(4)] for j in range(3)]
    for i in range(0, 3):
        d[i][0] = b[i]
        d[i][1] = c[i][0]
        d[i][2] = c[i][1]
        d[i][3] = c[i][2]
    for i in range(0, 3):
        c[i] = c[i] / math.sqrt(c[i][0] * c[i][0] + c[i][1] * c[i][1] + c[i][2] * c[i][2])
    d.sort()
    temp=[0 for i in range(4)]
    temp=d[2]
    d[2]=d[0]
    d[0]=temp
    print('The principal stresses and their principal direction vectors are:',end=' ')
    for i in range(0, 3):
        print('\nPrincipal stress[{0}] ={1} Pa.'.format(i+1, d[i][0]))
        print('Its principal direction vector:')
        for j in range(1, 4):
            print('{}'.format(d[i][j]), end=' ')
    r1=(d[0][0]-d[2][0])/2
    r2=(d[0][0]-d[1][0])/2
    r3=(d[1][0]-d[2][0])/2
    theta=np.linspace(0,2*np.pi,3600)
    x1=(d[0][0]+d[2][0])/2+r1*np.cos(theta)
    y1=r1*np.sin(theta)
    x2=(d[0][0]+d[1][0])/2+r2*np.cos(theta)
    y2=r2*np.sin(theta)
    x3=(d[1][0]+d[2][0])/2+r3*np.cos(theta)
    y3=r3*np.sin(theta)
    plt.plot(x1,y1)
    plt.plot(x2,y2)
    plt.plot(x3,y3)
    plt.show()
    return 0
def principal_strain(a):
    b, c = np.linalg.eig(a)
    d = [[0 for i in range(4)] for j in range(3)]
    for i in range(0, 3):
        d[i][0] = b[i]
        d[i][1] = c[i][0]
        d[i][2] = c[i][1]
        d[i][3] = c[i][2]
    for i in range(0, 3):
        c[i] = c[i] / math.sqrt(c[i][0] * c[i][0] + c[i][1] * c[i][1] + c[i][2] * c[i][2])
    d.sort()
    temp = [0 for i in range(4)]
    temp = d[2]
    d[2] = d[0]
    d[0] = temp
    print('The principal strains and their principal direction vectors are:',end=' ')
    for i in range(0, 3):
        print('\nPrincipal strain[{0}] ={1}'.format(i+1, d[i][0]))
        print('Its principal direction vector:')
        for j in range(1, 4):
            print('{}'.format(d[i][j]), end=' ')
    r1 = (d[0][0] - d[2][0]) / 2
    r2 = (d[0][0] - d[1][0]) / 2
    r3 = (d[1][0] - d[2][0]) / 2
    theta = np.linspace(0, 2 * np.pi, 3600)
    x1 = (d[0][0] + d[2][0]) / 2 + r1 * np.cos(theta)
    y1 = r1 * np.sin(theta)
    x2 = (d[0][0] + d[1][0]) / 2 + r2 * np.cos(theta)
    y2 = r2 * np.sin(theta)
    x3 = (d[1][0] + d[2][0]) / 2 + r3 * np.cos(theta)
    y3 = r3 * np.sin(theta)
    plt.plot(x1, y1)
    plt.plot(x2, y2)
    plt.plot(x3, y3)
    plt.show()
    return 0
def tranformation_from_stress_to_strain(a,E,v):
     G=E/(2+2*v)
     eps_x=(a[0][0]-v*a[1][1]-v*a[2][2])/E
     eps_y=(a[1][1]-v*a[2][2]-v*a[0][0])/E
     eps_z=(a[2][2]-v*a[0][0]-v*a[1][1])/E
     gam__xy=a[0][1]/(2*G)
     gam__xz=a[0][2]/(2*G)
     gam__yz=a[1][2]/(2*G)
     b=np.zeros((3,3))
     b[0][0]=eps_x
     b[1][1]=eps_y
     b[2][2]=eps_z
     b[0][1]=gam__xy
     b[1][0]=gam__xy
     b[0][2]=gam__xz
     b[2][0]=gam__xz
     b[1][2]=gam__yz
     b[2][1]=gam__yz
     return b
def tranformation_from_strain_to_stress(a, E, v):
    G = E / (2 + 2 * v)
    sig_x=((1-v*v)*a[0][0]+(v+v*v)*(a[1][1]+a[2][2]))*E/(1-2*v*v*v-3*v*v)
    sig_y=((1-v*v)*a[1][1]+(v+v*v)*(a[0][0]+a[2][2]))*E/(1-2*v*v*v-3*v*v)
    sig_z=((1-v*v)*a[2][2]+(v+v*v)*(a[1][1]+a[0][0]))*E/(1-2*v*v*v-3*v*v)
    tau_xy=G*2*a[0][1]
    tau_xz=G*2*a[0][2]
    tau_yz=G*2*a[1][2]
    b = np.zeros((3, 3))
    b[0][0] = sig_x
    b[1][1] = sig_y
    b[2][2] = sig_z
    b[0][1] = tau_xy
    b[1][0] = tau_xy
    b[0][2] = tau_xz
    b[2][0] = tau_xz
    b[1][2] = tau_yz
    b[2][1] = tau_yz
    return b
type_judge=1
while type_judge:
    judge = input('What do you know? Stress(ess) or strain(ain)?')
    if judge!='ess'and judge!='ain':
        print('Retry please.')
    else:
        type_judge=0
if judge=='ess':
    print('Please complete the stress matrix :')
    stress = np.zeros((3, 3))
    sig_x = float(input('sigma_x='))
    stress[0][0] = sig_x
    sig_y = float(input('sigma_y='))
    stress[1][1] = sig_y
    sig_z = float(input('sigma_z='))
    stress[2][2] = sig_z
    tau_xy=float(input('tau_xy='))
    stress[0][1]=tau_xy
    stress[1][0]=tau_xy
    tau_xz=float(input('tau_xz='))
    stress[0][2]=tau_xz
    stress[2][0]=tau_xz
    tau_yz=float(input('tau_yz='))
    stress[1][2]=tau_yz
    stress[2][1]=tau_yz
    print('The stress matrix is :')
    print(stress)
    print('The direction vector for the stress you want to solve is :')
    dir_vec=np.zeros((1,3))
    for i in range(0,3):
        dir_vec[0][i]=float(input())
    solve_for_given_stress(dir_vec,stress)
    principal_stress(stress)
    print('\nThe strain metrix is :\n',tranformation_from_stress_to_strain(stress,E,v))
    print('The direction vector for the strain you want to solve is :')
    dir_vec_1 = np.zeros((1, 3))
    for i in range(0, 3):
        dir_vec_1[0][i] = float(input())
    solve_for_given_strain(dir_vec_1,tranformation_from_stress_to_strain(stress,E,v))
    principal_strain(tranformation_from_stress_to_strain(stress,E,v))
if judge=='ain':
    print('Please complete the strain matrix :')
    strain = np.zeros((3, 3))
    eps_x = float(input('epsilon_x='))
    strain[0][0] = eps_x
    eps_y = float(input('epsilon_y='))
    strain[1][1] = eps_y
    eps_z = float(input('epsilon_z='))
    strain[2][2] = eps_z
    gam__xy = float(input('gamma__xy='))
    strain[0][1] = gam__xy
    strain[1][0] = gam__xy
    gam__xz = float(input('gamma__xz='))
    strain[0][2] = gam__xz
    strain[2][0] = gam__xz
    gam__yz = float(input('gamma__yz='))
    strain[1][2] = gam__yz
    strain[2][1] = gam__yz
    print('The strain matrix is :')
    print(strain)
    print('The direction vector for the strain you want to solve is :')
    dir_vec = np.zeros((1, 3))
    for i in range(0,3):
        dir_vec[0][i]=float(input())
    solve_for_given_strain(dir_vec, strain)
    principal_strain(strain)
    print('\nThe stress metrix is :\n',tranformation_from_strain_to_stress(strain,E,v))
    print('The direction vector for the stress you want to solve is :')
    dir_vec_1 = np.zeros((1, 3))
    for i in range(0, 3):
        dir_vec_1[0][i] = float(input())
    solve_for_given_stress(dir_vec_1, tranformation_from_strain_to_stress(strain,E,v))
    principal_stress(tranformation_from_strain_to_stress(strain,E,v))
