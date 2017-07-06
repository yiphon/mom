#(c)  5.14.2017  Hang Yin
#Transformation between stress and strain.

from tkinter import *
import numpy as np
import math
import matplotlib.pyplot as plt
#一个字符高约19像素，宽10像素

root=Tk()
root.title('Transformation between strain and stress')
t1=Text(root,height=6,width=60,font=('Times New Roman',12))   #此处高和宽是字符高宽的倍数，不是像素
t1.place(x=0,y=0)
t1.insert(1.0,'Poisson ratio is (-1 < v < 0.5) :\n\n')
t1.insert(END,'Young\'s modulus is (Pa) :\n\n')
t1.insert(END,'What do you know? Stress(ess) or strain(ain)?\n')

en1=Entry(root,width=8,bd=5)    #字符高宽。
en1.place(x=200,y=0)
en2=Entry(root,width=15,bd=5)
en2.place(x=200,y=35)
en3=Entry(root,width=5,bd=5)
en3.place(x=300,y=75)

def setParameter():  #设置参数。
    t2.delete('0.0',END)
    v=en1.get()
    E=en2.get()
    MType=en3.get()
    v=float(v)
    E=float(E)

    def completeStressMetrix():  #完成应力矩阵。
        stress = np.zeros((3, 3))
        sigma_x=enEss1.get()
        stress[0][0]=float(sigma_x)
        sigma_y=enEss2.get()
        stress[1][1]=float(sigma_y)
        sigma_z=enEss3.get()
        stress[2][2]=float(sigma_z)
        tau_xy=enEss4.get()
        stress[0][1]=float(tau_xy)
        stress[1][0]=float(tau_xy)
        tau_xz=enEss5.get()
        stress[0][2]=float(tau_xz)
        stress[2][0]=float(tau_xz)
        tau_yz=enEss6.get()
        stress[2][1]=float(tau_yz)
        stress[1][2]=float(tau_yz)

        t3.delete('0.0',END)
        t3.insert(0.0,'The stress matrix is:\n\n')
        t3.insert(END,'\t%15f\t\t%15f\t\t%15f\n\n'%(stress[0][0],stress[0][1],stress[0][2]))
        t3.insert(END,'\t%15f\t\t%15f\t\t%15f\n\n'%(stress[1][0],stress[1][1],stress[1][2]))
        t3.insert(END,'\t%15f\t\t%15f\t\t%15f\n'%(stress[2][0],stress[2][1],stress[2][2]))
        t3.insert(END,'\n The direction vector for the stress you want to solve is (press Enter,  split with comma) :\n\n\n')
        var=StringVar()
        var.set('0.0')
        endir1=Entry(root,width=30,bd=5,textvariable=var)
        endir1.place(x=110,y=565)

        def solveForGivenStress():  #求解给定平面上的正应力与切应力。
            endir1.delete(0,END)

            def return_endir1(en): #将回车与程序的下一部分建立联系，使程序继续下去。
                t3.delete('0.0',END)
                t3.insert(0.0, 'The stress matrix is:\n\n')
                t3.insert(END, '\t%15f\t\t%15f\t\t%15f\n\n' % (stress[0][0], stress[0][1], stress[0][2]))
                t3.insert(END, '\t%15f\t\t%15f\t\t%15f\n\n' % (stress[1][0], stress[1][1], stress[1][2]))
                t3.insert(END, '\t%15f\t\t%15f\t\t%15f\n' % (stress[2][0], stress[2][1], stress[2][2]))
                t3.insert(END,
                          '\n The direction vector for the stress you want to solve is (press Enter,  split with comma) :\n\n\n')
                dir_vec_ess=np.zeros((1,3))
                dir_vec1=endir1.get()
                dir_vec1=dir_vec1.split(',')
                for i in range(3):
                    dir_vec1[i]=float(dir_vec1[i])
                    dir_vec_ess[0][i]=dir_vec1[i]

                length=math.sqrt(dir_vec_ess[0][0]*dir_vec_ess[0][0]+dir_vec_ess[0][1]*dir_vec_ess[0][1]+dir_vec_ess[0][2]*dir_vec_ess[0][2])
                for i in range(0,3):
                    dir_vec_ess[0][i]=dir_vec_ess[0][i]/length
                stress_x=stress[0][0]*dir_vec_ess[0][0]+stress[0][1]*dir_vec_ess[0][1]+stress[0][2]*dir_vec_ess[0][2]
                stress_y = stress[1][0] * dir_vec_ess[0][0] + stress[1][1] * dir_vec_ess[0][1] + stress[1][2] * dir_vec_ess[0][2]
                stress_z = stress[2][0] * dir_vec_ess[0][0] + stress[2][1] * dir_vec_ess[0][1] + stress[2][2] * dir_vec_ess[0][2]
                stress_square=stress_x*stress_x+stress_y*stress_y+stress_z*stress_z
                normal_stress=stress[0][0]*dir_vec_ess[0][0]*dir_vec_ess[0][0]+stress[1][1]*dir_vec_ess[0][1]*dir_vec_ess[0][1]+stress[2][2]*dir_vec_ess[0][2]*dir_vec_ess[0][2]+\
                    2*stress[0][1]*dir_vec_ess[0][0]*dir_vec_ess[0][1]+2*stress[0][2]*dir_vec_ess[0][0]*dir_vec_ess[0][2]+2*stress[2][1]*dir_vec_ess[0][2]*dir_vec_ess[0][1]
                try:
                    shear_stress=math.sqrt(stress_square-math.pow(normal_stress,2))
                except ValueError:
                    shear_stress=0
                t3.insert(END,'\tThe normal stress is {} Pa.\n'.format(normal_stress))
                t3.insert(END,'\tThe shear stress is {} Pa.\n\n'.format(shear_stress))

                def principal_stress():  #求解主应力与主方向，最大切应力及所在平面方向。
                    t4.delete('0.0',END)
                    b, c = np.linalg.eig(stress)
                    d = [[0 for i in range(4)] for j in range(3)]
                    for i in range(0, 3):
                        d[i][0] = b[i]
                        d[i][1] = c[i][0]
                        d[i][2] = c[i][1]
                        d[i][3] = c[i][2]
                    t4.insert(1.0, '\tThe principal stresses and their principal direction vectors are:\n\n')
                    for i in range(0, 3):
                        c[i] = c[i] / math.sqrt(c[i][0] * c[i][0] + c[i][1] * c[i][1] + c[i][2] * c[i][2])
                    d.sort()
                    temp=[0 for i in range(4)]
                    temp=d[2]
                    d[2]=d[0]
                    d[0]=temp
                    for i in range(0, 3):
                        t4.insert(END,'\tPrincipal stress[{0}] ={1} Pa.\n'.format(i+1, d[i][0]))
                        t4.insert(END,'\tIts principal direction vector:\n')
                        for j in range(1, 4):
                            t4.insert(END,'\t{}\n'.format(d[i][j]))
                        t4.insert(END,'\n')
                    tau=d[0][0]-d[2][0]
                    t4.insert(END,'\tMaximum shear stress is {0} Pa.\n'.format(tau/2))
                    t4.insert(END, '\tIts principal direction vector:\n')
                    for i in range(0,3):
                        t4.insert(END,'\t{}\n'.format(-(d[0][i+1]+d[2][i+1])/math.sqrt(2)))
                def tranformation_from_stress_to_strain():  #实现应力应变矩阵的转换。
                    t5.delete('0.0',END)
                    G = E / (2 + 2 * v)
                    eps_x = (stress[0][0] - v * stress[1][1] - v * stress[2][2]) / E
                    eps_y = (stress[1][1] - v * stress[2][2] - v * stress[0][0]) / E
                    eps_z = (stress[2][2] - v * stress[0][0] - v * stress[1][1]) / E
                    gam__xy = stress[0][1] / (2 * G)
                    gam__xz = stress[0][2] / (2 * G)
                    gam__yz = stress[1][2] / (2 * G)
                    strain = np.zeros((3, 3))
                    strain[0][0] = eps_x
                    strain[1][1] = eps_y
                    strain[2][2] = eps_z
                    strain[0][1] = gam__xy
                    strain[1][0] = gam__xy
                    strain[0][2] = gam__xz
                    strain[2][0] = gam__xz
                    strain[1][2] = gam__yz
                    strain[2][1] = gam__yz
                    t5.insert(0.0, 'The strain metrix is :\n\n')
                    t5.insert(END, '%15f\t\t%15f\t\t%15f\n\n' % (strain[0][0], strain[0][1], strain[0][2]))
                    t5.insert(END, '%15f\t\t%15f\t\t%15f\n\n' % (strain[1][0], strain[1][1], strain[1][2]))
                    t5.insert(END, '%15f\t\t%15f\t\t%15f\n' % (strain[2][0], strain[2][1], strain[2][2]))
                    t5.insert(END,'\nThe direction vector for the strain you want to solve is (press Enter, split with comma) :\n\n\n')
                    var = StringVar()
                    var.set('0.0')
                    endir2 = Entry(root, width=30, bd=5, textvariable=var)
                    endir2.place(x=600, y=565)

                    def solveForGivenStrain():  #求解给定平面上的正应力与切应力。
                        endir2.delete(0,END)

                        def return_endir2(en):  #将回车与程序的下一部分建立联系，使程序继续进行。
                            t5.delete('0.0',END)
                            t5.insert(0.0, 'The strain metrix is :\n\n')
                            t5.insert(END, '%15f\t\t%15f\t\t%15f\n\n' % (strain[0][0], strain[0][1], strain[0][2]))
                            t5.insert(END, '%15f\t\t%15f\t\t%15f\n\n' % (strain[1][0], strain[1][1], strain[1][2]))
                            t5.insert(END, '%15f\t\t%15f\t\t%15f\n' % (strain[2][0], strain[2][1], strain[2][2]))
                            t5.insert(END,
                                      '\nThe direction vector for the strain you want to solve is (press Enter, split with comma) :\n\n\n')
                            dir_vec_ain = np.zeros((1, 3))
                            dir_vec2 = endir2.get()
                            dir_vec2 = dir_vec2.split(',')
                            for i in range(3):
                                dir_vec2[i] = float(dir_vec2[i])
                                dir_vec_ain[0][i] = dir_vec2[i]
                            length = math.sqrt(
                                dir_vec_ain[0][0] * dir_vec_ain[0][0] + dir_vec_ain[0][1] * dir_vec_ain[0][1] +
                                dir_vec_ain[0][2] * dir_vec_ain[0][2])
                            for i in range(0, 3):
                                dir_vec_ain[0][i] = dir_vec_ain[0][i] / length
                            strain_x = strain[0][0] * dir_vec_ain[0][0] + strain[0][1] * dir_vec_ain[0][1] + strain[0][2]*dir_vec_ain[0][2]
                            strain_y = strain[1][0] * dir_vec_ain[0][0] + strain[1][1] * dir_vec_ain[0][1] + strain[1][2] * dir_vec_ain[0][2]
                            strain_z = strain[2][0] * dir_vec_ain[0][0] + strain[2][1] * dir_vec_ain[0][1] + strain[2][2]*dir_vec_ain[0][2]
                            strain_square = strain_x * strain_x + strain_y * strain_y + strain_z * strain_z
                            normal_strain = strain[0][0] * dir_vec_ain[0][0] * dir_vec_ain[0][0] + strain[1][1] * \
                                                                                                   dir_vec_ain[0][1] * \
                                                                                                   dir_vec_ain[0][1] + \
                                            strain[2][2] * dir_vec_ain[0][2] * dir_vec_ain[0][2] + \
                                            2 * strain[0][1] * dir_vec_ain[0][0] * dir_vec_ain[0][1] + 2 * strain[0][
                                2] * dir_vec_ain[0][0] * dir_vec_ain[0][2] + 2 * strain[2][1] * dir_vec_ain[0][2] * \
                                                                             dir_vec_ain[0][1]

                            try:
                                shear_strain = math.sqrt(strain_square - math.pow(normal_strain, 2))
                            except ValueError:
                                shear_strain = 0
                            shear__strain = 2 * shear_strain
                            t5.insert(END, '\tThe normal strain is {}.\n'.format(normal_strain))
                            t5.insert(END, '\tThe shear strain is {}.\n\n'.format(shear__strain))

                            def principal_strain():  #求解主应变与主方向，最大切应变与所在平面方向。
                                t6.delete('0.0', END)
                                b, c = np.linalg.eig(strain)
                                e = [[0 for i in range(4)] for j in range(3)]
                                for i in range(0, 3):
                                    e[i][0] = b[i]
                                    e[i][1] = c[i][0]
                                    e[i][2] = c[i][1]
                                    e[i][3] = c[i][2]
                                for i in range(0, 3):
                                    c[i] = c[i] / math.sqrt(c[i][0] * c[i][0] + c[i][1] * c[i][1] + c[i][2] * c[i][2])
                                e.sort()
                                temp = [0 for i in range(4)]
                                temp = e[2]
                                e[2] = e[0]
                                e[0] = temp
                                t6.insert(END, '\tThe principal strains and their principal direction vectors are:\n\n')
                                for i in range(0, 3):
                                    t6.insert(END, '\tPrincipal strain[{0}] ={1}.\n'.format(i + 1, e[i][0]))
                                    t6.insert(END, '\tIts principal direction vector:\n')
                                    for j in range(1, 4):
                                        t6.insert(END, '\t{}\n'.format(e[i][j]))
                                    t6.insert(END, '\n')
                                gamma = (e[0][0] - e[2][0])*2
                                t6.insert(END, '\tMaximum shear strain is {0}.\n'.format(gamma/2))
                                t6.insert(END, '\tIts principal direction vector:\n')
                                for i in range(0, 3):
                                    t6.insert(END, '\t{}\n'.format(-(e[0][i + 1] + e[2][i + 1]) / math.sqrt(2)))

                            principal_strain()


                            def visualize():  #实现可视化。
                                b, c = np.linalg.eig(strain)
                                e = [[0 for i in range(4)] for j in range(3)]
                                for i in range(0, 3):
                                    e[i][0] = b[i]
                                    e[i][1] = c[i][0]
                                    e[i][2] = c[i][1]
                                    e[i][3] = c[i][2]
                                for i in range(0, 3):
                                    c[i] = c[i] / math.sqrt(c[i][0] * c[i][0] + c[i][1] * c[i][1] + c[i][2] * c[i][2])
                                e.sort()
                                temp = [0 for i in range(4)]
                                temp = e[2]
                                e[2] = e[0]
                                e[0] = temp
                                fig1=plt.figure('Strain circles')
                                r1 = (e[0][0] - e[2][0]) / 2
                                r2 = (e[0][0] - e[1][0]) / 2
                                r3 = (e[1][0] - e[2][0]) / 2
                                str1 = str(e[0][0])
                                str2 = str(e[1][0])
                                str3 = str(e[2][0])
                                strr1 = str(r1)
                                strr2 = str(r2)
                                strr3 = str(r3)
                                theta = np.linspace(0, 2 * np.pi, 3600)
                                x1 = (e[0][0] + e[2][0]) / 2 + r1 * np.cos(theta)
                                y1 = r1 * np.sin(theta)
                                x2 = (e[0][0] + e[1][0]) / 2 + r2 * np.cos(theta)
                                y2 = r2 * np.sin(theta)
                                x3 = (e[1][0] + e[2][0]) / 2 + r3 * np.cos(theta)
                                y3 = r3 * np.sin(theta)
                                g1 = plt.plot(x1, y1, color='red')
                                g2 = plt.plot(x2, y2, color='blue')
                                g3 = plt.plot(x3, y3, color='green')
                                plt.gca().invert_yaxis()
                                plt.plot(x1, y1, color="red", linewidth=1.5, linestyle="-", label='Circle 1')
                                plt.plot(x2, y2, color="blue", linewidth=1.5, linestyle="-", label='Circle 2')
                                plt.plot(x3, y3, color="green", linewidth=1.5, linestyle="-", label='Circle 3')
                                plt.legend(loc='upper left')
                                plt.axis('equal')
                                plt.xlabel('epsilon_x')
                                plt.ylabel('gamma_xy/2')
                                plt.title('Strain circles')
                                axes=plt.subplot(111)
                                axes.spines['top'].set_color('none')
                                axes.spines['right'].set_color('none')
                                axes = plt.subplot(111)
                                axes.spines['top'].set_color('none')
                                axes.spines['right'].set_color('none')
                                axes.annotate('A', xy=(e[0][0], 0), xytext=(e[0][0] - 0.09 * r1, 0), fontsize=15)
                                axes.plot(e[0][0], 0, 'ok')
                                axes.annotate('B', xy=(e[1][0], 0), xytext=(e[1][0] + 0.03 * r1, 0), fontsize=15)
                                axes.plot(e[1][0], 0, 'ok')
                                axes.annotate('C', xy=(e[2][0], 0), xytext=(e[2][0] - 0.09 * r1, 0), fontsize=15)
                                axes.plot(e[2][0], 0, 'ok')
                                axes.annotate('D', xy=((e[0][0] + e[2][0]) / 2, r1),
                                              xytext=((e[0][0] + e[2][0]) / 2, 1.09 * r1), fontsize=15)
                                axes.plot((e[0][0] + e[2][0]) / 2, r1, 'ok')
                                axes.annotate('E', xy=((e[0][0] + e[1][0]) / 2, r1),
                                              xytext=((e[0][0] + e[1][0]) / 2, r2 - 0.03 * r1), fontsize=15)
                                axes.plot((e[0][0] + e[1][0]) / 2, r2, 'ok')
                                axes.annotate('F', xy=((e[2][0] + e[1][0]) / 2, r3),
                                              xytext=((e[2][0] + e[1][0]) / 2, r3 - 0.03 * r1), fontsize=15)
                                axes.plot((e[1][0] + e[2][0]) / 2, r3, 'ok')
                                ax = fig1.add_subplot(111)
                                ax.text(0.7, 0.99, 'Circle 1\n' + 'Center (' + str((e[0][0] + e[2][
                                    0]) / 2) + ', 0)\n' + 'A (' + str1 + ', 0)\n' + 'C (' + str3 + ', 0)\n' + 'D (' + str(
                                    (e[0][0] + e[2][0]) / 2) + ', ' + strr1 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='red', fontsize=8)
                                ax.text(0.7, 0.2, 'Circle 2\n' + 'Center (' + str((e[0][0] + e[1][
                                    0]) / 2) + ', 0)\n' + 'A (' + str1 + ', 0)\n' + 'B (' + str2 + ', 0)\n' + 'E (' + str(
                                    (e[0][0] + e[1][0]) / 2) + ', ' + strr2 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='blue', fontsize=8)
                                ax.text(0.01, 0.2, 'Circle 3\n' + 'Center (' + str((e[1][0] + e[2][
                                    0]) / 2) + ', 0)\n' + 'B ( ' + str2 + ', 0)\n' + 'C (' + str3 + ', 0)\n' + 'F (' + str(
                                    (e[1][0] + e[2][0]) / 2) + ', ' + strr3 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='green', fontsize=8)
                                figManager = plt.get_current_fig_manager()
                                figManager.resize(*figManager.window.maxsize())

                                b, c = np.linalg.eig(stress)
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
                                fig2 = plt.figure('Stress circles')
                                r1 = (d[0][0] - d[2][0]) / 2
                                r2 = (d[0][0] - d[1][0]) / 2
                                r3 = (d[1][0] - d[2][0]) / 2
                                str1 = str(d[0][0])
                                str2 = str(d[1][0])
                                str3 = str(d[2][0])
                                strr1 = str(r1)
                                strr2 = str(r2)
                                strr3 = str(r3)
                                theta = np.linspace(0, 2 * np.pi, 3600)
                                x1 = (d[0][0] + d[2][0]) / 2 + r1 * np.cos(theta)
                                y1 = r1 * np.sin(theta)
                                x2 = (d[0][0] + d[1][0]) / 2 + r2 * np.cos(theta)
                                y2 = r2 * np.sin(theta)
                                x3 = (d[1][0] + d[2][0]) / 2 + r3 * np.cos(theta)
                                y3 = r3 * np.sin(theta)
                                g1 = plt.plot(x1, y1, color='red')
                                g2 = plt.plot(x2, y2, color='blue')
                                g3 = plt.plot(x3, y3, color='green')
                                plt.gca().invert_yaxis()
                                plt.plot(x1, y1, color="red", linewidth=1.5, linestyle="-", label='Circle 1')
                                plt.plot(x2, y2, color="blue", linewidth=1.5, linestyle="-", label='Circle 2')
                                plt.plot(x3, y3, color="green", linewidth=1.5, linestyle="-", label='Circle 3')
                                plt.legend(loc='upper left')
                                plt.axis('equal')
                                plt.xlabel('sigma_x/Pa')
                                plt.ylabel('tau_xy/Pa')
                                plt.title('Stress circles')
                                axes = plt.subplot(111)
                                axes.spines['top'].set_color('none')
                                axes.spines['right'].set_color('none')
                                axes.annotate('A', xy=(d[0][0], 0), xytext=(d[0][0] - 0.09 * r1, 0), fontsize=15)
                                axes.plot(d[0][0], 0, 'ok')
                                axes.annotate('B', xy=(d[1][0], 0), xytext=(d[1][0] + 0.03 * r1, 0), fontsize=15)
                                axes.plot(d[1][0], 0, 'ok')
                                axes.annotate('C', xy=(d[2][0], 0), xytext=(d[2][0] - 0.09 * r1, 0), fontsize=15)
                                axes.plot(d[2][0], 0, 'ok')
                                axes.annotate('D', xy=((d[0][0] + d[2][0]) / 2, r1),
                                              xytext=((d[0][0] + d[2][0]) / 2, 1.09 * r1), fontsize=15)
                                axes.plot((d[0][0] + d[2][0]) / 2, r1, 'ok')
                                axes.annotate('E', xy=((d[0][0] + d[1][0]) / 2, r1),
                                              xytext=((d[0][0] + d[1][0]) / 2, r2 - 0.03 * r1), fontsize=15)
                                axes.plot((d[0][0] + d[1][0]) / 2, r2, 'ok')
                                axes.annotate('F', xy=((d[2][0] + d[1][0]) / 2, r3),
                                              xytext=((d[2][0] + d[1][0]) / 2, r3 - 0.03 * r1), fontsize=15)
                                axes.plot((d[1][0] + d[2][0]) / 2, r3, 'ok')
                                ax = fig2.add_subplot(111)
                                ax.text(0.7, 0.99, 'Circle 1\n' + 'Center (' + str((d[0][0] + d[2][
                                    0]) / 2) + ', 0)\n' + 'A (' + str1 + ', 0)\n' + 'C (' + str3 + ', 0)\n' + 'D (' + str(
                                    (d[0][0] + d[2][0]) / 2) + ', ' + strr1 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='red', fontsize=8)
                                ax.text(0.7, 0.2, 'Circle 2\n' + 'Center (' + str((d[0][0] + d[1][
                                    0]) / 2) + ', 0)\n' + 'A (' + str1 + ', 0)\n' + 'B (' + str2 + ', 0)\n' + 'E (' + str(
                                    (d[0][0] + d[1][0]) / 2) + ', ' + strr2 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='blue', fontsize=8)
                                ax.text(0.01, 0.2, 'Circle 3\n' + 'Center (' + str((d[1][0] + d[2][
                                    0]) / 2) + ', 0)\n' + 'B ( ' + str2 + ', 0)\n' + 'C (' + str3 + ', 0)\n' + 'F (' + str(
                                    (d[1][0] + d[2][0]) / 2) + ', ' + strr3 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='green', fontsize=8)
                                figManager = plt.get_current_fig_manager()
                                figManager.resize(*figManager.window.maxsize())

                                plt.show()

                            bu3 = Button(root, width=15, text='Visualize', command=visualize)
                            bu3.place(x=1050, y=450)


                        endir2.bind('<Return>',return_endir2)

                    solveForGivenStrain()

                principal_stress()
                tranformation_from_stress_to_strain()

            endir1.bind('<Return>', return_endir1)

        solveForGivenStress()

#以下求解过程为上述过程的对偶过程，不再注释。

    def completeStrainMetrix():
        strain = np.zeros((3, 3))
        epsilon_x = enEss1.get()
        strain[0][0] = float(epsilon_x)
        epsilon_y = enEss2.get()
        strain[1][1] = float(epsilon_y)
        epsilon_z = enEss3.get()
        strain[2][2] = float(epsilon_z)
        gamma__xy = enEss4.get()
        strain[0][1] =0.5* float(gamma__xy)
        strain[1][0] =0.5* float(gamma__xy)
        gamma__xz = enEss5.get()
        strain[0][2] = 0.5*float(gamma__xz)
        strain[2][0] =0.5* float(gamma__xz)
        gamma__yz = enEss6.get()
        strain[2][1] =0.5* float(gamma__yz)
        strain[1][2] = 0.5*float(gamma__yz)

        t3.delete('0.0', END)
        t3.insert(0.0, 'The strain matrix is:\n\n')
        t3.insert(END, '\t%15f\t\t%15f\t\t%15f\n\n' % (strain[0][0], strain[0][1], strain[0][2]))
        t3.insert(END, '\t%15f\t\t%15f\t\t%15f\n\n' % (strain[1][0], strain[1][1], strain[1][2]))
        t3.insert(END, '\t%15f\t\t%15f\t\t%15f\n' % (strain[2][0], strain[2][1], strain[2][2]))
        t3.insert(END, '\n  The direction vector for the strain you want to solve is (press Enter,  split with comma) :\n\n\n')
        var = StringVar()
        var.set('0.0')
        endir1 = Entry(root, width=30, bd=5, textvariable=var)
        endir1.place(x=110, y=565)

        def solveForGivenStrain():
            endir1.delete(0, END)

            def return_endir1(en):
                t3.delete('0.0', END)
                t3.insert(0.0, 'The strain matrix is:\n\n')
                t3.insert(END, '\t%15f\t\t%15f\t\t%15f\n\n' % (strain[0][0], strain[0][1], strain[0][2]))
                t3.insert(END, '\t%15f\t\t%15f\t\t%15f\n\n' % (strain[1][0], strain[1][1], strain[1][2]))
                t3.insert(END, '\t%15f\t\t%15f\t\t%15f\n' % (strain[2][0], strain[2][1], strain[2][2]))
                t3.insert(END,
                          '\n  The direction vector for the strain you want to solve is (press Enter,  split with comma) :\n\n\n')
                dir_vec_ain = np.zeros((1, 3))
                dir_vec2 = endir1.get()
                dir_vec2 = dir_vec2.split(',')
                for i in range(3):
                    dir_vec2[i] = float(dir_vec2[i])
                    dir_vec_ain[0][i] = dir_vec2[i]

                length = math.sqrt(
                    dir_vec_ain[0][0] * dir_vec_ain[0][0] + dir_vec_ain[0][1] * dir_vec_ain[0][1] +
                        dir_vec_ain[0][2] * dir_vec_ain[0][2])
                for i in range(0, 3):
                    dir_vec_ain[0][i] = dir_vec_ain[0][i] / length
                strain_x = strain[0][0] * dir_vec_ain[0][0] + strain[0][1] * dir_vec_ain[0][1] + strain[0][2] * dir_vec_ain[0][2]
                strain_y = strain[1][0] * dir_vec_ain[0][0] + strain[1][1] * dir_vec_ain[0][1] + strain[1][2] * dir_vec_ain[0][2]
                strain_z = strain[2][0] * dir_vec_ain[0][0] + strain[2][1] * dir_vec_ain[0][1] + strain[2][2] * dir_vec_ain[0][2]
                strain_square = strain_x * strain_x + strain_y * strain_y + strain_z * strain_z
                normal_strain =strain[0][0]*dir_vec_ain[0][0]*dir_vec_ain[0][0]+strain[1][1]*dir_vec_ain[0][1]*dir_vec_ain[0][1]+strain[2][2]*dir_vec_ain[0][2]*dir_vec_ain[0][2]+2*strain[0][1]*dir_vec_ain[0][0]*dir_vec_ain[0][1]+2*strain[0][2]*dir_vec_ain[0][0]*dir_vec_ain[0][2]+2*strain[2][1]*dir_vec_ain[0][2]*dir_vec_ain[0][1]
                try:
                    shear_strain = math.sqrt(strain_square - math.pow(normal_strain, 2))
                except ValueError:
                    shear_strain = 0
                shear__strain=2*shear_strain
                t3.insert(END, '\tThe normal strain is {}.\n'.format(normal_strain))
                t3.insert(END, '\tThe shear strain is {}.\n\n'.format(shear__strain))

                def principal_strain():
                    t4.delete('0.0', END)
                    b, c = np.linalg.eig(strain)
                    d = [[0 for i in range(4)] for j in range(3)]
                    for i in range(0, 3):
                        d[i][0] = b[i]
                        d[i][1] = c[i][0]
                        d[i][2] = c[i][1]
                        d[i][3] = c[i][2]
                    t4.insert(1.0, 'The principal strains and their principal direction vectors are:\n\n')
                    for i in range(0, 3):
                        c[i] = c[i] / math.sqrt(c[i][0] * c[i][0] + c[i][1] * c[i][1] + c[i][2] * c[i][2])
                    d.sort()
                    temp = [0 for i in range(4)]
                    temp = d[2]
                    d[2] = d[0]
                    d[0] = temp
                    for i in range(0, 3):
                        t4.insert(END, '\tPrincipal strain[{0}] ={1}.\n'.format(i + 1, d[i][0]))
                        t4.insert(END, '\tIts principal direction vector:\n')
                        for j in range(1, 4):
                            t4.insert(END, '\t{}\n'.format(d[i][j]))
                        t4.insert(END, '\n')
                    gamma = (d[0][0] - d[2][0]) * 2
                    t4.insert(END, '\tMaximum shear strain is {0}.\n'.format(gamma / 2))
                    t4.insert(END, '\tIts principal direction vector:\n')
                    for i in range(0, 3):
                        t4.insert(END, '\t{}\n'.format(-(d[0][i + 1] + d[2][i + 1]) / math.sqrt(2)))

                def tranformation_from_strain_to_stress():
                    t5.delete('0.0', END)
                    G = E / (2 + 2 * v)
                    sig_x = ((1 - v * v) * strain[0][0] + (v + v * v) * (strain[1][1] + strain[2][2])) * E / (
                    1 - 2 * v * v * v - 3 * v * v)
                    sig_y = ((1 - v * v) * strain[1][1] + (v + v * v) * (strain[0][0] + strain[2][2])) * E / (
                    1 - 2 * v * v * v - 3 * v * v)
                    sig_z = ((1 - v * v) * strain[2][2] + (v + v * v) * (strain[1][1] + strain[0][0])) * E / (
                    1 - 2 * v * v * v - 3 * v * v)
                    tau_xy = G * 2 * strain[0][1]
                    tau_xz = G * 2 * strain[0][2]
                    tau_yz = G * 2 * strain[1][2]
                    stress = np.zeros((3, 3))
                    stress[0][0] = sig_x
                    stress[1][1] = sig_y
                    stress[2][2] = sig_z
                    stress[0][1] = tau_xy
                    stress[1][0] = tau_xy
                    stress[0][2] = tau_xz
                    stress[2][0] = tau_xz
                    stress[1][2] = tau_yz
                    stress[2][1] = tau_yz
                    t5.insert(0.0, 'The stress metrix is :\n\n')
                    t5.insert(END, '%15f\t\t%15f\t\t%15f\n\n' % (stress[0][0], stress[0][1], stress[0][2]))
                    t5.insert(END, '%15f\t\t%15f\t\t%15f\n\n' % (stress[1][0], stress[1][1], stress[1][2]))
                    t5.insert(END, '%15f\t\t%15f\t\t%15f\n' % (stress[2][0], stress[2][1], stress[2][2]))
                    t5.insert(END,'\nThe direction vector for the stress you want to solve is (press Enter, split with comma) :\n\n\n')
                    var = StringVar()
                    var.set('0.0')
                    endir2 = Entry(root, width=30, bd=5, textvariable=var)
                    endir2.place(x=600, y=565)

                    def solveForGivenStress():
                        endir2.delete(0, END)

                        def return_endir2(en):
                            t5.delete('0.0',END)
                            t5.insert(0.0, 'The stress metrix is :\n\n')
                            t5.insert(END, '%15f\t\t%15f\t\t%15f\n\n' % (stress[0][0], stress[0][1], stress[0][2]))
                            t5.insert(END, '%15f\t\t%15f\t\t%15f\n\n' % (stress[1][0], stress[1][1], stress[1][2]))
                            t5.insert(END, '%15f\t\t%15f\t\t%15f\n' % (stress[2][0], stress[2][1], stress[2][2]))
                            t5.insert(END,
                                      '\nThe direction vector for the stress you want to solve is (press Enter, split with comma) :\n\n\n')
                            dir_vec_ess = np.zeros((1, 3))
                            dir_vec2 = endir2.get()
                            dir_vec2 = dir_vec2.split(',')
                            for i in range(3):
                                dir_vec2[i] = float(dir_vec2[i])
                                dir_vec_ess[0][i] = dir_vec2[i]
                            length = math.sqrt(
                                dir_vec_ess[0][0] * dir_vec_ess[0][0] + dir_vec_ess[0][1] * dir_vec_ess[0][1] +
                                dir_vec_ess[0][2] * dir_vec_ess[0][2])
                            for i in range(0, 3):
                                dir_vec_ess[0][i] = dir_vec_ess[0][i] / length
                            stress_x = stress[0][0] * dir_vec_ess[0][0] + stress[0][1] * dir_vec_ess[0][1] + \
                                       stress[0][2] * dir_vec_ess[0][2]
                            stress_y = stress[1][0] * dir_vec_ess[0][0] + stress[1][1] * dir_vec_ess[0][1] + \
                                       stress[1][2] * dir_vec_ess[0][2]
                            stress_z = stress[2][0] * dir_vec_ess[0][0] + stress[2][1] * dir_vec_ess[0][1] + \
                                       stress[2][2] * dir_vec_ess[0][2]
                            stress_square = stress_x * stress_x + stress_y * stress_y + stress_z * stress_z
                            normal_stress = stress[0][0] * dir_vec_ess[0][0] * dir_vec_ess[0][0] + stress[1][
                                                                                                       1] * \
                                                                                                   dir_vec_ess[
                                                                                                       0][1] * \
                                                                                                   dir_vec_ess[
                                                                                                       0][1] + \
                                            stress[2][2] * dir_vec_ess[0][2] * dir_vec_ess[0][2] + \
                                            2 * stress[0][1] * dir_vec_ess[0][0] * dir_vec_ess[0][1] + 2 * \
                                                                                                       stress[
                                                                                                           0][
                                                                                                           2] * \
                                                                                                       dir_vec_ess[
                                                                                                           0][
                                                                                                           0] * \
                                                                                                       dir_vec_ess[
                                                                                                           0][
                                                                                                           2] + 2 * \
                                                                                                                stress[
                                                                                                                    2][
                                                                                                                    1] * \
                                                                                                                dir_vec_ess[
                                                                                                                    0][
                                                                                                                    2] * \
                                                                                                                dir_vec_ess[
                                                                                                                    0][
                                                                                                                    1]

                            try:
                                shear_stress = math.sqrt(stress_square - math.pow(normal_stress, 2))
                            except ValueError:
                                shear_stress = 0

                            t5.insert(END, '\tThe normal stress is {}Pa.\n'.format(normal_stress))
                            t5.insert(END, '\tThe shear stress is {}Pa.\n\n'.format(shear_stress))

                            def principal_stress():
                                t6.delete('0.0', END)
                                b, c = np.linalg.eig(stress)
                                e = [[0 for i in range(4)] for j in range(3)]
                                for i in range(0, 3):
                                    e[i][0] = b[i]
                                    e[i][1] = c[i][0]
                                    e[i][2] = c[i][1]
                                    e[i][3] = c[i][2]
                                for i in range(0, 3):
                                    c[i] = c[i] / math.sqrt(
                                        c[i][0] * c[i][0] + c[i][1] * c[i][1] + c[i][2] * c[i][2])
                                e.sort()
                                temp = [0 for i in range(4)]
                                temp = e[2]
                                e[2] = e[0]
                                e[0] = temp
                                t6.insert(END, 'The principal stresses and their principal direction vectors are:\n\n')
                                for i in range(0, 3):
                                    t6.insert(END, '\tPrincipal stress[{0}] ={1}Pa.\n'.format(i + 1, e[i][0]))
                                    t6.insert(END, '\tIts principal direction vector:\n')
                                    for j in range(1, 4):
                                        t6.insert(END, '\t{}Pa\n'.format(e[i][j]))
                                    t6.insert(END, '\n')
                                tau = (e[0][0] - e[2][0]) /1
                                t6.insert(END, '\tMaximum shear stress is {0} Pa.\n'.format(tau / 2))
                                t6.insert(END, '\tIts principal direction vector:\n')
                                for i in range(0, 3):
                                    t6.insert(END, '\t{}\n'.format(-(e[0][i + 1] + e[2][i + 1]) / math.sqrt(2)))

                            principal_stress()


                            def visualize():
                                b, c = np.linalg.eig(stress)
                                e = [[0 for i in range(4)] for j in range(3)]
                                for i in range(0, 3):
                                    e[i][0] = b[i]
                                    e[i][1] = c[i][0]
                                    e[i][2] = c[i][1]
                                    e[i][3] = c[i][2]
                                for i in range(0, 3):
                                    c[i] = c[i] / math.sqrt(
                                        c[i][0] * c[i][0] + c[i][1] * c[i][1] + c[i][2] * c[i][2])
                                e.sort()
                                temp = [0 for i in range(4)]
                                temp = e[2]
                                e[2] = e[0]
                                e[0] = temp
                                fig1 = plt.figure('Stress circles')
                                r1 = (e[0][0] - e[2][0]) / 2
                                r2 = (e[0][0] - e[1][0]) / 2
                                r3 = (e[1][0] - e[2][0]) / 2
                                str1 = str(e[0][0])
                                str2 = str(e[1][0])
                                str3 = str(e[2][0])
                                strr1 = str(r1)
                                strr2 = str(r2)
                                strr3 = str(r3)
                                theta = np.linspace(0, 2 * np.pi, 3600)
                                x1 = (e[0][0] + e[2][0]) / 2 + r1 * np.cos(theta)
                                y1 = r1 * np.sin(theta)
                                x2 = (e[0][0] + e[1][0]) / 2 + r2 * np.cos(theta)
                                y2 = r2 * np.sin(theta)
                                x3 = (e[1][0] + e[2][0]) / 2 + r3 * np.cos(theta)
                                y3 = r3 * np.sin(theta)
                                g1 = plt.plot(x1, y1, color='red')
                                g2 = plt.plot(x2, y2, color='blue')
                                g3 = plt.plot(x3, y3, color='green')
                                plt.gca().invert_yaxis()
                                plt.plot(x1, y1, color="red", linewidth=1.5, linestyle="-", label='Circle 1')
                                plt.plot(x2, y2, color="blue", linewidth=1.5, linestyle="-", label='Circle 2')
                                plt.plot(x3, y3, color="green", linewidth=1.5, linestyle="-", label='Circle 3')
                                plt.legend(loc='upper left')
                                plt.xlabel('sigma_x/Pa')
                                plt.ylabel('tau_xy/Pa')
                                plt.title('Stress circles')
                                plt.axis('equal')
                                axes = plt.subplot(111)
                                axes.spines['top'].set_color('none')
                                axes.spines['right'].set_color('none')
                                axes.annotate('A', xy=(e[0][0], 0), xytext=(e[0][0] - 0.09 * r1, 0), fontsize=15)
                                axes.plot(e[0][0], 0, 'ok')
                                axes.annotate('B', xy=(e[1][0], 0), xytext=(e[1][0] + 0.03 * r1, 0), fontsize=15)
                                axes.plot(e[1][0], 0, 'ok')
                                axes.annotate('C', xy=(e[2][0], 0), xytext=(e[2][0] - 0.09 * r1, 0), fontsize=15)
                                axes.plot(e[2][0], 0, 'ok')
                                axes.annotate('D', xy=((e[0][0] + e[2][0]) / 2, r1),
                                              xytext=((e[0][0] + e[2][0]) / 2, 1.09 * r1), fontsize=15)
                                axes.plot((e[0][0] + e[2][0]) / 2, r1, 'ok')
                                axes.annotate('E', xy=((e[0][0] + e[1][0]) / 2, r1),
                                              xytext=((e[0][0] + e[1][0]) / 2, r2 - 0.03 * r1), fontsize=15)
                                axes.plot((e[0][0] + e[1][0]) / 2, r2, 'ok')
                                axes.annotate('F', xy=((e[2][0] + e[1][0]) / 2, r3),
                                              xytext=((e[2][0] + e[1][0]) / 2, r3 - 0.03 * r1), fontsize=15)
                                axes.plot((e[1][0] + e[2][0]) / 2, r3, 'ok')
                                ax = fig1.add_subplot(111)
                                ax.text(0.7, 0.99, 'Circle 1\n' + 'Center (' + str((e[0][0] + e[2][
                                    0]) / 2) + ', 0)\n' + 'A (' + str1 + ', 0)\n' + 'C (' + str3 + ', 0)\n' + 'D (' + str(
                                    (e[0][0] + e[2][0]) / 2) + ', ' + strr1 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='red', fontsize=8)
                                ax.text(0.7, 0.2, 'Circle 2\n' + 'Center (' + str((e[0][0] + e[1][
                                    0]) / 2) + ', 0)\n' + 'A (' + str1 + ', 0)\n' + 'B (' + str2 + ', 0)\n' + 'E (' + str(
                                    (e[0][0] + e[1][0]) / 2) + ', ' + strr2 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='blue', fontsize=8)
                                ax.text(0.01, 0.2, 'Circle 3\n' + 'Center (' + str((e[1][0] + e[2][
                                    0]) / 2) + ', 0)\n' + 'B ( ' + str2 + ', 0)\n' + 'C (' + str3 + ', 0)\n' + 'F (' + str(
                                    (e[1][0] + e[2][0]) / 2) + ', ' + strr3 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='green', fontsize=8)

                                figManager = plt.get_current_fig_manager()
                                figManager.resize(*figManager.window.maxsize())

                                b, c = np.linalg.eig(strain)
                                d = [[0 for i in range(4)] for j in range(3)]
                                for i in range(0, 3):
                                    d[i][0] = b[i]
                                    d[i][1] = c[i][0]
                                    d[i][2] = c[i][1]
                                    d[i][3] = c[i][2]

                                for i in range(0, 3):
                                    c[i] = c[i] / math.sqrt(
                                        c[i][0] * c[i][0] + c[i][1] * c[i][1] + c[i][2] * c[i][2])
                                d.sort()
                                temp = [0 for i in range(4)]
                                temp = d[2]
                                d[2] = d[0]
                                d[0] = temp
                                fig2 = plt.figure('Strain circles')
                                r1 = (d[0][0] - d[2][0]) / 2
                                r2 = (d[0][0] - d[1][0]) / 2
                                r3 = (d[1][0] - d[2][0]) / 2
                                str1=str(d[0][0])
                                str2=str(d[1][0])
                                str3=str(d[2][0])
                                strr1=str(r1)
                                strr2=str(r2)
                                strr3=str(r3)
                                theta = np.linspace(0, 2 * np.pi, 3600)
                                x1 = (d[0][0] + d[2][0]) / 2 + r1 * np.cos(theta)
                                y1 = r1 * np.sin(theta)
                                x2 = (d[0][0] + d[1][0]) / 2 + r2 * np.cos(theta)
                                y2 = r2 * np.sin(theta)
                                x3 = (d[1][0] + d[2][0]) / 2 + r3 * np.cos(theta)
                                y3 = r3 * np.sin(theta)
                                g1=plt.plot(x1, y1,color='red')
                                g2=plt.plot(x2, y2,color='blue')
                                g3=plt.plot(x3, y3,color='green')
                                plt.gca().invert_yaxis()
                                plt.plot(x1, y1, color="red", linewidth=1.5, linestyle="-",label='Circle 1')
                                plt.plot(x2, y2, color="blue", linewidth=1.5, linestyle="-",label='Circle 2')
                                plt.plot(x3, y3, color="green", linewidth=1.5, linestyle="-",label='Circle 3')
                                plt.legend(loc='upper left')
                                plt.xlabel('epsilon_x/Pa')
                                plt.ylabel('gamma_xy/2')
                                plt.title('Strain circles')
                                plt.axis('equal')
                                axes = plt.subplot(111)
                                axes.spines['top'].set_color('none')
                                axes.spines['right'].set_color('none')
                                axes.annotate('A', xy=(d[0][0], 0), xytext=(d[0][0]-0.09*r1,0 ),fontsize=15)
                                axes.plot(d[0][0],0,'ok')
                                axes.annotate('B', xy=(d[1][0], 0), xytext=(d[1][0]+0.03*r1, 0),fontsize=15)
                                axes.plot(d[1][0], 0, 'ok')
                                axes.annotate('C', xy=(d[2][0], 0), xytext=(d[2][0]-0.09*r1, 0),fontsize=15)
                                axes.plot(d[2][0], 0, 'ok')
                                axes.annotate('D', xy=((d[0][0] + d[2][0]) / 2, r1), xytext=((d[0][0] + d[2][0]) / 2, 1.09*r1),fontsize=15)
                                axes.plot((d[0][0] + d[2][0]) / 2, r1, 'ok')
                                axes.annotate('E', xy=((d[0][0] + d[1][0]) / 2, r1),
                                              xytext=((d[0][0] + d[1][0]) / 2, r2-0.03*r1),fontsize=15)
                                axes.plot((d[0][0] + d[1][0]) / 2, r2,'ok')
                                axes.annotate('F', xy=((d[2][0] + d[1][0]) / 2, r3),
                                              xytext=((d[2][0] + d[1][0]) / 2, r3 - 0.03 * r1),fontsize=15)
                                axes.plot((d[1][0] + d[2][0]) / 2, r3, 'ok')
                                ax=fig2.add_subplot(111)
                                ax.text(0.7, 0.99, 'Circle 1\n'+'Center ('+str((d[0][0] + d[2][0]) / 2)+', 0)\n'+'A ('+str1+', 0)\n'+'C ('+str3+', 0)\n'+'D ('+str((d[0][0] + d[2][0]) / 2)+', '+strr1+')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='red',fontsize=8)
                                ax.text(0.7, 0.2, 'Circle 2\n'+'Center ('+str((d[0][0] + d[1][0]) / 2)+', 0)\n'+'A ('+str1+', 0)\n'+'B ('+str2+', 0)\n'+'E ('+str((d[0][0] + d[1][0]) / 2)+', '+strr2+')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='blue',fontsize=8)
                                ax.text(0.01, 0.2, 'Circle 3\n' + 'Center (' + str((d[1][0] + d[2][
                                    0]) / 2) + ', 0)\n' + 'B ( ' + str2 + ', 0)\n' + 'C (' + str3 + ', 0)\n' + 'F (' + str(
                                    (d[1][0] + d[2][0]) / 2) + ', ' + strr3 + ')',
                                        verticalalignment='top', horizontalalignment='left',
                                        transform=ax.transAxes,
                                        color='green', fontsize=8)
                                figManager = plt.get_current_fig_manager()
                                figManager.resize(*figManager.window.maxsize())
                                plt.show()

                            bu3 = Button(root, width=15, text='Visualize', command=visualize)
                            bu3.place(x=1050, y=450)

                        endir2.bind('<Return>', return_endir2)

                    solveForGivenStress()

                principal_strain()
                tranformation_from_strain_to_stress()

            endir1.bind('<Return>',return_endir1)

        solveForGivenStrain()


    if MType == 'ess':  #针对应力模式的界面初始设置。
        t2.insert(1.0,'\t\tComplete the stress metrix:\n')
        t2.insert(END,'\t\tsigma_x =\n\n')
        t2.insert(END,'\t\tsigma_y =\n\n')
        t2.insert(END,'\t\tsigma_z =\n\n')
        t2.insert(END,'\t\ttau_xy =\n\n')
        t2.insert(END,'\t\ttau_xz =\n\n')
        t2.insert(END,'\t\ttau_yz =\n\n')
        enEss1=Entry(root,width=10,bd=5)
        enEss1.place(x=200,y=140)
        enEss2=Entry(root,width=10,bd=5)
        enEss2.place(x=200,y=178)
        enEss3=Entry(root,width=10,bd=5)
        enEss3.place(x=200,y=216)
        enEss4=Entry(root,width=10,bd=5)
        enEss4.place(x=200,y=254)
        enEss5=Entry(root,width=10,bd=5)
        enEss5.place(x=200,y=292)
        enEss6=Entry(root,width=10,bd=5)
        enEss6.place(x=200,y=330)
        bu2=Button(root,width=12,text='Complete',command=completeStressMetrix)
        bu2.place(x=300,y=330)


    if MType == 'ain':  #针对应变模式的界面初始设置。
        t2.insert(1.0,'\t\tComplete the strain metrix:\n')
        t2.insert(END,'\t           epsilon_x =\n\n')
        t2.insert(END,'\t           epsilon_y =\n\n')
        t2.insert(END,'\t           epsilon_z =\n\n')
        t2.insert(END,'\t         gamma_xy =\n\n')
        t2.insert(END,'\t         gamma_xz =\n\n')
        t2.insert(END,'\t         gamma_yz =\n\n')
        enEss1=Entry(root,width=10,bd=5)
        enEss1.place(x=200,y=140)
        enEss2=Entry(root,width=10,bd=5)
        enEss2.place(x=200,y=178)
        enEss3=Entry(root,width=10,bd=5)
        enEss3.place(x=200,y=216)
        enEss4=Entry(root,width=10,bd=5)
        enEss4.place(x=200,y=254)
        enEss5=Entry(root,width=10,bd=5)
        enEss5.place(x=200,y=292)
        enEss6=Entry(root,width=10,bd=5)
        enEss6.place(x=200,y=330)
        bu2=Button(root,width=12,text='Complete',command=completeStrainMetrix)
        bu2.place(x=300,y=330)


#界面公共初始设置。
bu1=Button(root,width=8,text='SET',command=setParameter)
bu1.place(x=400,y=75)

t2=Text(root,height=13,width=60,font=('Times New Roman',12))
t2.place(x=0,y=120)  #此处为像素。

t3=Text(root,height=15,width=60,font=('Times New Roman',12))
t3.place(x=0,y=370)

t4=Text(root,height=26,width=68,font=('Times New Roman',8))
t4.place(x=485,y=0)

t5=Text(root,height=15,width=55,font=('Times New Roman',12))
t5.place(x=485,y=370)

t6=Text(root,height=26,width=68,font=('Times New Roman',8))
t6.place(x=900,y=0)

buquit=Button(root,width=8,text='Quit',command=quit)
buquit.place(x=1075,y=550)
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
root.geometry("%dx%d+0+0" % (w, h))
root.mainloop()

