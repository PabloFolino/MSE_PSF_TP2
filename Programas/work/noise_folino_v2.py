# Autor : Pablo D. Folino
# Ejercitación lee información de una archivo "noise_folino_v1.py"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sc
import os
import simpleaudio as sa
from numpy import genfromtxt

# Valores del enunciado
Fs  =   8000        # frecuencia de muestreo, desconozco con que
Ts  =   1/Fs        # frecuencia se obtuvieron los datos
fs   =  8000        # Frecuencia de muestreo en Hz, ver frecuencias soportadas de
                    # la place de sonido


def reproducir(note):
    audio = note.astype(np.int16)               #tranforma la variable note a entero de 16bits y lo guarda en audio
    play_obj = sa.play_buffer(audio, 1, 2, fs)  # sale el audio
    play_obj.wait_done()                        # espera que termine la linea anterior

def X_ifft_conFiltro(X_fft,filtro):
    #--------------------Señal filtrada en frecuencia---------------------------
    X_fft_filtro=np.copy(X_fft)
    F_filtro=np.arange(-len(X_fft)/2,len(X_fft)/2-Fs/len(X_fft))
    centro=int(len(X_fft)/2)                                    # Acordrse que N=len(X_fft)
    # for i in range(0,centro-filtro):
    #     X_fft_filtro[i]=0
    #     X_fft_filtro[i+centro+filtro]=0
    for i in range(int(-filtro/2),int(filtro/2)):
        X_fft_filtro[i+centro]=0
    M_fft_filtro= abs(X_fft_filtro)**2/len(X_fft)
    #----------Discreto------------------
    ts_filtro = Ts*np.arange(0,len(X_fft),1)           
    #-------Anti Tranformada-------------
    señal_filtro = np.fft.ifft(X_fft_filtro)
    return ts_filtro,señal_filtro,F_filtro,X_fft_filtro,M_fft_filtro


# Se leen los archivo
tData=np.load("../Informe/Archivos de enunciados/chapu_noise.npy")[::1]
hData=(genfromtxt("./xx/filtro_pasabajos.txt",skip_header=1)).astype(np.int16)

os.system("clear")

# Se calcula las longitudes
N=len(tData)
M=len(hData)
size_vector=M+N-1

# Se verifica que se leyeron los N elementos
print("\t Se pudieron leer los N={} elementos del archivo".format(N))
print("\t El formato de los elementos es={}".format(type(tData[0])))
print("\t Se pudieron leer los M={} elementos del filtro".format(M))
print("\t El formato de los elementos es={}".format(type(hData[0])))
print("\t La longitud de los vectores a antitransformar es={}".format(size_vector))


# Datos extendidos
tData_ext=np.zeros(size_vector)
hData_ext=np.zeros(size_vector)
for i in range(0,N,1):
    tData_ext[i]=tData[i]
for i in range(0,M,1):
    hData_ext[i]=hData[i]

# Se reproduce la señal original
consulta=input("\nDesea reproducir la señal S o N[Enter] :")    
if consulta=='S' or consulta =='s':
    reproducir(tData)


# Genero el vector de estados temporales
n_Data      = np.arange(0,N,1) #arranco con numeros enteros para evitar errores de float
t_Data      = n_Data*Ts

# Transformada
f_Data=np.arange(-N/2,N/2+M-1)
fftData  = np.fft.fft(tData_ext)
ffthData = np.fft.fft(hData_ext)
YData=fftData*ffthData/10000                # Multiplico en frecuencia

ifftData = np.fft.ifft(YData)               # Antitransformada
fftData  = np.fft.fftshift(fftData)
YData    = np.fft.fftshift(YData)

salida=np.zeros(N)
for i in range(0,N,1):
    salida[i]=ifftData[i]


consulta=input("\nDesea reproducir la señal filtrada S o N[Enter] :")    
if consulta=='S' or consulta =='s':
    reproducir((np.real(salida)).astype(np.int16) )

# Se grafica para probar las señales
fig = plt.figure(1)
plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.95, wspace=0.4, hspace=0.8)


    

#--------------------Señal original-------------------------------------
s1 = fig.add_subplot(4,2,1)
plt.title("Señal de entrada x_fft")
plt.xlabel("Tiempo (seg)")
plt.ylabel("Amplitud x_fft")
s1.grid(True)
s1.plot(t_Data, tData,"g-")

s9 = fig.add_subplot(4,2,2)
plt.title("Valores")
plt.xlim(0,10)
plt.ylim(0,10)
plt.axis('off')
s9.spines['right'].set_visible(False)
s9.spines['top'].set_visible(False)
s9.spines['bottom'].set_visible(False)
s9.spines['left'].set_visible(False)
plt.text(0,9,"Se supone que FrecSampl.="+str(Fs)+"Hz",fontsize=10)
plt.text(0,7,"Resolución temporal="+str(f'{Ts*1000:.{2}f}')+"mseg",fontsize=10)
plt.text(0,5,"Resolución espectral="+str(f'{Fs/N:.{4}f}')+"Hz",fontsize=10)
plt.text(0,3,"Muestras="+str(N),fontsize=10)


s1 = fig.add_subplot(4,2,3)
plt.title("Transformada de Fouirer (PC)")
plt.xlabel("Frecuencia en Hz")
plt.ylabel("Amplitud |X_fft|")
s1.grid(True)
s1.plot(f_Data,np.abs(fftData)**2/N,'r-')

s1 = fig.add_subplot(4,2,5)
plt.title("Transformada filtrada fc="+"Hz (PC)" )
plt.xlabel("Frecuencia en Hz")
plt.ylabel("Amplitud |X_fft|")
s1.grid(True)
s1.plot(f_Data,np.abs(YData)**2/N,'r-')

s1 = fig.add_subplot(4,2,7)
plt.title("Transformada filtrada con fc="+"Hz (PC)")
plt.xlabel("Frecuencia en Hz")
plt.ylabel("Amplitud |X_fft|")
s1.grid(True)
s1.plot(t_Data,salida,'r-')

plt.get_current_fig_manager().window.showMaximized() #para QT5
plt.show()


