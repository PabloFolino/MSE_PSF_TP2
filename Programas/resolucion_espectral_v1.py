# Autor : Pablo D. Folino
# Ejercitación lee información de una archivo "resolucion_espectral.txt"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sc
import os

# Valores del enunciado
Fs  =   200      # frecuencia de muestreo
Ts  =   1/Fs
N   =   100

# Se utiliza la técnica de zero padding, se agregan M1 veces N
# la nueva resolución espectral es Fs/(M1+1).N
M1   = 9


# Se lee el archivo
with open("../Informe/Archivos de enunciados/resolucion_espectral.txt","r") as ins:
    cont = ins.read()  # Esto devuelve el contenido completo, no linea por linea
    arr = eval(cont)
    señal=arr

os.system("clear")

# Se verifica que se leyeron los N elementos
if(len(señal)==N):
    print("\tSe leyeron N={} elementos del archivo".format(N))
    print("\tEl tipo de elemento es={}".format(type(señal)))
    print("\tEl rimer elemento del archivo es={}".format(señal[1]))
    print("\tEl tipo del primer elemento es={}".format(type(señal[1])))
else:
    print("\t ERROR-No se pudieron leer los N={} elementos del archivo".format(N))

#--------------------Señal original-------------------------------------
#----------Discreto------------------
ts = Ts*np.arange(0, N,1)
#----------Tranformada---------------
fft = np.fft.fft(señal)/N                        
M_fft = np.fft.fftshift(abs(fft))                            
F=np.fft.fftfreq(int(Fs/2),Ts)  
F=np.fft.fftshift(F) 
#----------Potencia promedio---------
pot_promedio=np.sum(M_fft**2)

#--------------------Señal con Zero Padding-----------------------------
N_zp=(M1+1)*N
señal_zp=np.zeros(N_zp)
for i in range(0, len(señal)):
    señal_zp[i]=señal[i]

#----------Discreto------------------
ts_zp = Ts*np.arange(0,N_zp,1)
#----------Tranformada---------------
fft_zp = np.fft.fft(señal_zp)/N_zp                        
M_fft_zp = np.fft.fftshift(abs(fft_zp))                            
F_zp=np.fft.fftfreq(int(Fs*(M1+1)/2),Ts)  
F_zp=np.fft.fftshift(F_zp) 
#----------Potencia promedio---------
pot_promedio_zp=np.sum(M_fft_zp**2)





# Se grafica para probar las señales
fig = plt.figure(1)
plt.subplots_adjust(left=None, bottom=0.05, right=None, top=0.95, wspace=0.4, hspace=0.5)

#--------------------Señal original-------------------------------------
s1 = fig.add_subplot(3,2,1)
plt.title("Señal de entrada")
plt.xlabel("Tiempo")
plt.ylabel("Amplitud")
s1.grid(True)
plt.xlim(0,N/Fs)
s1.plot(ts,señal,'ro')

s5 = fig.add_subplot(3,2,3)
plt.title("Transformada de la señal de entrada")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Amplitud FFT")
s5.grid(True)
plt.xlim(-Fs/2,Fs/2-Fs/N)
#s5.plot(F, M_fft)
s5.stem(F, M_fft, markerfmt=" ", basefmt="-b")


s9 = fig.add_subplot(3,2,5)
plt.title("Valores")
plt.xlim(0,10)
plt.ylim(0,10)
plt.axis('off')
s9.spines['right'].set_visible(False)
s9.spines['top'].set_visible(False)
s9.spines['bottom'].set_visible(False)
s9.spines['left'].set_visible(False)
plt.text(0,9,"Resolución temporal="+str(f'{Ts*1000:.{2}f}')+"mseg",fontsize=10)
plt.text(0,7,"Resolución espectral="+str(f'{Fs/N:.{2}f}')+"Hz",fontsize=10)
plt.text(0,5,"Pot prom="+str(f'{pot_promedio:.{2}f}')+"w",fontsize=10)
plt.text(0,3,"FrecSampl.="+str(Fs)+"Hz",fontsize=10)
plt.text(0,1,"Muestras="+str(N),fontsize=10)


#--------------------Señal con Zero Padding-----------------------------
s1 = fig.add_subplot(3,2,2)
plt.title("Señal de entrada Zero Padding")
plt.xlabel("Tiempo")
plt.ylabel("Amplitud")
s1.grid(True)
plt.xlim(0,N_zp/Fs)
s1.plot(ts_zp,señal_zp,'ro')

s5 = fig.add_subplot(3,2,4)
plt.title("Transformada de la señal con Zero Padding")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Amplitud FFT")
s5.grid(True)
plt.xlim(-Fs/2,Fs/2-Fs/N)
#s5.plot(F, M_fft)
s5.stem(F_zp, M_fft_zp, markerfmt=" ", basefmt="-b")


s9 = fig.add_subplot(3,2,6)
plt.title("Valores")
plt.xlim(0,10)
plt.ylim(0,10)
plt.axis('off')
s9.spines['right'].set_visible(False)
s9.spines['top'].set_visible(False)
s9.spines['bottom'].set_visible(False)
s9.spines['left'].set_visible(False)
plt.text(0,9,"Resolución temporal="+str(f'{Ts*1000:.{2}f}')+"mseg",fontsize=10)
plt.text(0,7,"Resolución espectral="+str(f'{Fs/N_zp:.{2}f}')+"Hz",fontsize=10)
plt.text(0,5,"Pot prom="+str(f'{pot_promedio_zp:.{2}f}')+"w",fontsize=10)
plt.text(0,3,"FrecSampl.="+str(Fs)+"Hz",fontsize=10)
plt.text(0,1,"Muestras="+str(N_zp)+"    Se incremtan en="+str(M1+1)+"veces",fontsize=10)


plt.get_current_fig_manager().window.showMaximized() #para QT5
plt.show()
