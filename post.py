import matplotlib.pyplot as plt
import numpy as np

#drying time plot
def drytime():
    T = np.array([200, 220, 250]) #C
    dry1 = np.array([580, 259.5, 92]) #ms
    dry3 = np.array([1996, 882, 307]) #ms
    dry5 = np.array([3436, 1515, 525.5]) #ms
    dry10 = np.array([7140, 3142, 1090])#ms

    plt.semilogy(T,dry1,T,dry3,T,dry5,T,dry10)
    plt.xlabel('T(C)')
    plt.ylabel('Drying time (ms)')
    plt.text(240,150,r'$1\mu m$')
    plt.text(230,700,r'$3\mu m$')
    plt.text(220,1500,r'$5\mu m$')
    plt.text(210,5000,r'$10\mu m$')
    plt.show()

#deposit comparison
def engineDeposit():
    depositPS = \
    np.array([  7.37430210e-12,   4.09507323e-12,   2.60070965e-13,
             3.21853653e-11,   2.90816554e-11,   2.45223436e-11,
             5.66219974e-11,   5.40940203e-11,   6.45934680e-11,
             1.17455286e-10,   1.16157368e-10,   2.05133000e-10])
    depositold = \
    np.array([[  2.26933280e-08,   1.74121690e-08,   5.72627386e-09],
           [  7.42649876e-08,   6.10429836e-08,   3.32431479e-08],
           [  1.26381366e-07,   1.06042398e-07,   7.61528222e-08],
           [  2.56982493e-07,   2.19898691e-07,   2.48522006e-07]])
    area = 2.199e-7 * 1e6 #mm2
    depositold = depositold/area
    depositPS = depositPS/area
    depositold = depositold.reshape((4,3))
    depositPS = depositPS.reshape((4,3))
    T=np.array([200,220,250])
    plt.semilogy(T,depositold.T,'-s')
    plt.xlabel('T(C)')
    plt.ylabel(r'Deposit Mass ($mg/mm^2$)')
    plt.text(240,5e-8,r'$1\mu m$')
    plt.text(230,2.5e-7,r'$3\mu m$')
    plt.text(220,5e-7,r'$5\mu m$')
    plt.text(210,1e-6,r'$10\mu m$')
    plt.show()
    plt.figure(2)
    T=np.array([200,220,250])
    plt.semilogy(T,depositPS.T,'-o')
    plt.xlabel('T(C)')
    plt.ylabel(r'Deposit Mass ($mg/mm^2$)')
    plt.text(240,5e-12,r'$1\mu m$')
    plt.text(230,1.5e-10,r'$3\mu m$')
    plt.text(220,2.5e-10,r'$5\mu m$')
    plt.text(210,5e-10,r'$10\mu m$')
    plt.show()

# deposit in lab test, 8Hz, 40h,52min (the pure fuel at 8Hz)
def TestDeposit():
    depositPS = \
    np.array([  7.52310311e-09,   7.63377165e-08,   1.92614759e-07,
         5.32859934e-07,   8.77460235e-07])
        
    depositold = \
    np.array([  2.62821762e-06,   2.55577507e-05,   6.46081374e-05,
         1.88039419e-04,   3.22856065e-04])
        
    area = 1.571e-5 * 1e6 #mm2
    # depositold = depositold/area
    # depositPS = depositPS/area
    h=np.array([1,3,5,10,15])
    plt.semilogy(h,depositold,'--o',h,depositPS,'-s')
    plt.xlabel('Film thickness ($\mu m$)')
    plt.ylabel(r'Deposit Mass (mg)')
    plt.show()

    
