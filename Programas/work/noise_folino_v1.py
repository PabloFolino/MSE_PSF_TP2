# Autor : Pablo D. Folino
# Ejercitación lee información de una archivo "noise_folino_v1.py"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sc
import os
import simpleaudio as sa

# Valores del enunciado
Fs  =   8000        # frecuencia de muestreo, desconozco con que
Ts  =   1/Fs        # frecuencia se obtuvieron los datos
fs   = 44100        # Frecuencia de muestreo en Hz, ver frecuencias soportadas de
                    # la place de sonido


def reproducir(note):
    audio = note.astype(np.int16)               #tranforma la variable note a entero de 16bits y lo guarda en audio
    play_obj = sa.play_buffer(audio, 1, 2, fs)  # sale el audio
    play_obj.wait_done()                        # espera que termine la linea anterior

def X_ifft_conFiltro(X_fft,filtro):
    #--------------------Señal filtrada en frecuencia---------------------------
    X_fft_filtro=np.zeros(len(X_fft),dtype=np.complex128)
    #X_fft_filtro=np.arange(filtro,dtype=np.complex128)
    F_filtro=np.arange(-Fs/2,Fs/2,Fs/len(X_fft))
    centro=int(len(X_fft)/2)                                    # Acordrse que N=len(X_fft)
    for i in range(0,filtro):
        X_fft_filtro[centro-int(filtro/2)+i]=X_fft[centro-int(filtro/2)+i]
    M_fft_filtro= abs(X_fft_filtro)/filtro
    #----------Potencia promedio---------
    pot_promedio_filtro=np.sum(M_fft_filtro**2)
    #----------Discreto------------------
    ts_filtro = Ts*np.arange(0,len(X_fft),1)           
    #-------Anti Tranformada-------------
    señal_filtro = np.fft.ifft(X_fft_filtro)
    return ts_filtro,señal_filtro,F_filtro,M_fft_filtro


# Se lee el archivo
x_fft=np.load("../Informe/Archivos de enunciados/chapu_noise.npy")[::1]

os.system("clear")

N=len(x_fft)

# Se verifica que se leyeron los N elementos
print("\t Se pudieron leer los N={} elementos del archivo".format(N))
print("\t El formato de los elementos es={}".format(type(x_fft[0])))

consulta=input("\nDesea reproducir la señal S o N[Enter] :")    
if consulta=='S' or consulta =='s':
    reproducir(x_fft)


# Genero el vector de estados temporales
nData      = np.arange(0,N,1) #arranco con numeros enteros para evitar errores de float
tData      = nData*Ts

# Potencia romedio calculada en el tiempo
pot_promedio=np.sum(np.abs(x_fft/N)**2)

fftData  = np.fft.fft(x_fft)
ifftData = np.fft.ifft(fftData)
fftData  = np.fft.fftshift(fftData)
fData=np.arange(-N/2,N/2-Fs/N)


delta=8000        # Frecuencias que se restan a cada lado el espectro
filtro=N-2*delta
ts_filtro,señal_filtro,F_filtro,M_fft_filtro=X_ifft_conFiltro(fftData,filtro)

# Se grafica para probar las señales
fig = plt.figure(1)
plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.95, wspace=0.4, hspace=0.8)

# consulta=input("\nDesea reproducir la señal filtrada S o N[Enter] :")    
# if consulta=='S' or consulta =='s':
#     nota1=np.real(señal_filtro)
#     nota1=int(nota1)
#     reproducir(nota1)
    

#--------------------Señal original-------------------------------------
s1 = fig.add_subplot(3,2,1)
plt.title("Señal de entrada x_fft")
plt.xlabel("Tiempo (seg)")
plt.ylabel("Amplitud x_fft")
s1.grid(True)
s1.plot(tData, x_fft,"g-")

s9 = fig.add_subplot(3,2,2)
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
plt.text(0,3,"Pot prom="+str(f'{pot_promedio:.{2}f}')+"w",fontsize=10)
plt.text(0,1,"Muestras="+str(N),fontsize=10)


s1 = fig.add_subplot(3,2,3)
plt.title("Transformada de Fouirer ")
plt.xlabel("Frecuencia en Hz")
plt.ylabel("Amplitud |X_fft|")
s1.grid(True)
s1.plot(fData,np.abs(fftData/N)**2,'r-')

s1 = fig.add_subplot(3,2,4)
plt.title("Transformada de Fouirer filtrada con"+str(filtro)+"muestras")
plt.xlabel("Frecuencia en Hz")
plt.ylabel("Amplitud |X_fft|")
s1.grid(True)
s1.plot(F_filtro,M_fft_filtro,'r-')

s1 = fig.add_subplot(3,2,5)
plt.title("Transformada de Fouirer filtrada con"+str(filtro)+"muestras")
plt.xlabel("Frecuencia en Hz")
plt.ylabel("Amplitud |X_fft|")
s1.grid(True)
s1.plot(ts_filtro,señal_filtro,'r-')



plt.get_current_fig_manager().window.showMaximized() #para QT5
plt.show()


