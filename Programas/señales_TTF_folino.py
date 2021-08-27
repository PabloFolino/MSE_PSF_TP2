# Autor : Pablo D. Folino
# Ejercitación de señales senoidales, cuadradas, triangulares y impulso

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sc

#Valores a probar
f   =   2       # frecuencia de la señal a realizar la transformada
Fs  =   40      # frecuencia de muestreo
Ts  =   1/Fs
fase =  0
N   =   40
amp =   2.5

Pu  =   100     # Cantidad de puntos para generar la continua




# Función senoidal
# Parámetros:
#       fs      --> fecuencia de sampleo
#       f       --> frecuencia de la señal de entrada
#       amp     --> amplitud del señal de 0 a 1.
#       muestras--> cantidad de veces que se repite la señal
#       fase    --> fase de la señal en radianes
# Devuelve:
#       ts      --> vector de tienpos de sampling
#       w1      --> vector de señal de la señal  
#       tsC     --> vector de tienpos de continuos
#       w1C     --> vector de valores contínuos de la señal
#       F       --> vector de frecuancia de samplig
#       M_fft   --> vector de valores de fracuancia de la señal 
def senoidal(fs,f,amp,muestras,fase):
    #----------Discreto------------------
    ts = Ts*np.arange(0, fs,1)
    w1 =amp*np.sin(2*np.pi*f*ts+fase)         
    #----------Contínuo------------------
    tsC = Ts*np.arange(0, fs,1/Pu)
    w1C = amp*np.sin(2*np.pi*f*tsC+fase)         
    #----------Tranformada---------------
    fft = np.fft.fft(w1)/muestras                        
    M_fft = np.fft.fftshift(abs(fft))                            
    #F = Fs*np.arange(0, len(w1))/len(w1) 
    F=np.fft.fftfreq(fs,1/muestras)  
    F=np.fft.fftshift(F) 
    #----------Potencia promedio---------
    pot_promedio=np.sum(M_fft**2)
    return ts,w1,tsC,w1C,F,M_fft,pot_promedio

# Función cuadrada
# Parámetros:
#       fs      --> fecuencia de sampleo
#       f       --> frecuencia de la señal de entrada
#       amp     --> amplitud del señal de 0 a 1.
#       muestras--> cantidad de veces que se repite la señal
# Devuelve:
#       ts      --> vector de tienpos de sampling
#       w1      --> vector de señal de la señal  
#       tsC     --> vector de tienpos de continuos
#       w1C     --> vector de valores contínuos de la señal
#       F       --> vector de frecuancia de samplig
#       M_fft   --> vector de valores de fracuancia de la señal     
def cuadrada(fs,f,amp,muestras,fase):
    #----------Discreto------------------
    ts = Ts*np.arange(0, fs,1)
    w1 =amp*sc.square(2*np.pi*ts*f+fase)        
    #----------Contínuo------------------
    tsC = Ts*np.arange(0, fs,1/Pu)
    w1C = amp*sc.square(2*np.pi*tsC*f+fase)        
    #----------Tranformada---------------
    fft = np.fft.fft(w1)/muestras                        
    M_fft = np.fft.fftshift(abs(fft))                            
    #F = Fs*np.arange(0, len(w1))/len(w1) 
    F=np.fft.fftfreq(fs,1/muestras)  
    F=np.fft.fftshift(F) 
    #----------Potencia promedio---------
    pot_promedio=np.sum(M_fft**2)
    return ts,w1,tsC,w1C,F,M_fft,pot_promedio


# Función triangular
# Parámetros:
#       fs      --> fecuencia de sampleo
#       f       --> frecuencia de la señal de entrada
#       amp     --> amplitud del señal de 0 a 1.
#       muestras--> cantidad de veces que se repite la señal
# Devuelve:
#       ts      --> vector de tienpos de sampling
#       w1      --> vector de señal de la señal  
#       tsC     --> vector de tienpos de continuos
#       w1C     --> vector de valores contínuos de la señal
#       F       --> vector de frecuancia de samplig
#       M_fft   --> vector de valores de fracuancia de la señal      amp*sc.sawtooth(2*np.pi*f*ts+fase,1)
def triangular(fs,f,amp,muestras,fase):
    #----------Discreto------------------
    ts = Ts*np.arange(0, fs,1)
    w1 =amp*sc.sawtooth(2*np.pi*f*ts+fase,1)        
    #----------Contínuo------------------
    tsC = Ts*np.arange(0, fs,1/Pu)
    w1C = amp*sc.sawtooth(2*np.pi*f*tsC+fase,1)        
    #----------Tranformada---------------
    fft = np.fft.fft(w1)/muestras                        
    M_fft = np.fft.fftshift(abs(fft))                            
    #F = Fs*np.arange(0, len(w1))/len(w1) 
    F=np.fft.fftfreq(fs,1/muestras)  
    F=np.fft.fftshift(F) 
    #----------Potencia promedio---------
    pot_promedio=np.sum(M_fft**2)
    return ts,w1,tsC,w1C,F,M_fft,pot_promedio

# Función impulso
# Parámetros:
#       fs      --> fecuencia de sampleo
#       f       --> frecuencia de la señal de entrada  T tiende a infinito
#       amp     --> amplitud del señal de 0 a 1.
# Devuelve:
#       ts      --> vector de tienpos de sampling
#       w1      --> vector de señal de la señal  
#       tsC     --> vector de tienpos de continuos
#       w1C     --> vector de valores contínuos de la señal
#       F       --> vector de frecuancia de samplig
#       M_fft   --> vector de valores de fracuancia de la señal      
def impulso(fs,f,amp,muestras,fase):
    #----------Discreto------------------
    ts = Ts*np.arange(0, fs,1)
    w1=np.zeros(len(ts))                            
    for i in range(1,len(ts)):
        if( ts[i] < 1/amp ):
            w1[i]=amp
        else:
            w1[i]=0
    #----------Contínuo------------------
    tsC = Ts*np.arange(0, fs,1/Pu)
    w1C=np.zeros(len(tsC))                            
    for i in range(1,len(tsC)):
        if( tsC[i] < 1/amp ):
            w1C[i]=amp
        else:
            w1C[i]=0
    #----------Tranformada---------------
    fft = np.fft.fft(w1)/muestras                        
    M_fft = np.fft.fftshift(abs(fft))                            
    #F = Fs*np.arange(0, len(w1))/len(w1) 
    F=np.fft.fftfreq(fs,1/muestras)  
    F=np.fft.fftshift(F) 
    #----------Potencia promedio---------
    pot_promedio=0
    #pot_promedio=np.sum(M_fft**2)
    return ts,w1,tsC,w1C,F,M_fft,pot_promedio




# Se grafica para probar las señales
fig = plt.figure(1)
plt.subplots_adjust(left=None, bottom=0.05, right=None, top=0.95, wspace=0.4, hspace=0.5)

#--------------------Senoidal-------------------------------------
s1 = fig.add_subplot(3,4,1)
plt.title("Senoidal")
plt.xlabel("Tiempo")
plt.ylabel("Amplitud")
s1.grid(True)
ts,w1,tsC,w1C,F,M_fft,pot_promedio=senoidal(Fs,f,amp,N,fase)
plt.xlim(0,N/Fs)
s1.plot(ts,w1,'ro')
s1.plot(tsC,w1C,'b-')


s5 = fig.add_subplot(3,4,5)
plt.title("Transformada")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Amplitud FFT")
s5.grid(True)
plt.xlim(-Fs/2,Fs/2-Fs/N)
#s5.plot(F, M_fft)
s5.stem(F, M_fft, markerfmt=" ", basefmt="-b")


s9 = fig.add_subplot(3,4,9)
plt.xlim(0,10)
plt.ylim(0,10)
plt.axis('off')
s9.spines['right'].set_visible(False)
s9.spines['top'].set_visible(False)
s9.spines['bottom'].set_visible(False)
s9.spines['left'].set_visible(False)
plt.text(0,9,"Frecuencia="+str(f)+"Hz",fontsize=10)
plt.text(0,7,"Amplitud="+str(f'{amp:.{2}f}')+"v",fontsize=10)
plt.text(0,5,"Pot prom="+str(f'{pot_promedio:.{2}f}')+"w",fontsize=10)
plt.text(0,3,"FrecSampl.="+str(Fs)+"Hz",fontsize=10)
plt.text(0,1,"Muestras="+str(N),fontsize=10)

#--------------------Cuadrada-------------------------------------
s2 = fig.add_subplot(3,4,2)
plt.title("Cuadrada")
plt.xlabel("Tiempo")
plt.ylabel("Amplitud")
s2.grid(True)
ts,w1,tsC,w1C,F,M_fft,pot_promedio=cuadrada(Fs,f,amp,N,fase)
plt.xlim(0,N/Fs)
s2.plot(ts,w1,'ro')
s2.plot(tsC,w1C,'b-')


s6 = fig.add_subplot(3,4,6)
plt.title("Transformada")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Amplitud FFT")
s6.grid(True)
plt.xlim(-Fs/2,Fs/2-Fs/N)
#s6.plot(F, M_fft)
s6.stem(F, M_fft, markerfmt=" ", basefmt="-b")

s10 = fig.add_subplot(3,4,10)
plt.xlim(0,10)
plt.ylim(0,10)
plt.axis('off')
s10.spines['right'].set_visible(False)
s10.spines['top'].set_visible(False)
s10.spines['bottom'].set_visible(False)
s10.spines['left'].set_visible(False)
plt.text(0,9,"Frecuencia="+str(f)+"Hz",fontsize=10)
plt.text(0,7,"Amplitud="+str(f'{amp:.{2}f}')+"v",fontsize=10)
plt.text(0,5,"Pot prom="+str(f'{pot_promedio:.{2}f}')+"w",fontsize=10)
plt.text(0,3,"FrecSampl.="+str(Fs)+"Hz",fontsize=10)
plt.text(0,1,"Muestras="+str(N),fontsize=10)


#--------------------Triangular-------------------------------------
s3 = fig.add_subplot(3,4,3)
plt.title("Triangular")
plt.xlabel("Tiempo")
plt.ylabel("Amplitud")
s3.grid(True)
ts,w1,tsC,w1C,F,M_fft,pot_promedio=triangular(Fs,f,amp,N,fase)
plt.xlim(0,N/Fs)
s3.plot(ts,w1,'ro')
s3.plot(tsC,w1C,'b-')


s7 = fig.add_subplot(3,4,7)
plt.title("Transformada")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Amplitud FFT")
s7.grid(True)
plt.xlim(-Fs/2,Fs/2-Fs/N)
#s7.plot(F, M_fft)
s7.stem(F, M_fft, markerfmt=" ", basefmt="-b")

s11 = fig.add_subplot(3,4,11)
plt.xlim(0,10)
plt.ylim(0,10)
plt.axis('off')
s11.spines['right'].set_visible(False)
s11.spines['top'].set_visible(False)
s11.spines['bottom'].set_visible(False)
s11.spines['left'].set_visible(False)
plt.text(0,9,"Frecuencia="+str(f)+"Hz",fontsize=10)
plt.text(0,7,"Amplitud="+str(f'{amp:.{2}f}')+"v",fontsize=10)
plt.text(0,5,"Pot prom="+str(f'{pot_promedio:.{2}f}')+"w",fontsize=10)
plt.text(0,3,"FrecSampl.="+str(Fs)+"Hz",fontsize=10)
plt.text(0,1,"Muestras="+str(N),fontsize=10)


#--------------------Impulso-------------------------------------
amp=5
Fs=40
N=40

s4 = fig.add_subplot(3,4,4)
plt.title("Impulso")
plt.xlabel("Tiempo")
plt.ylabel("Amplitud")
s4.grid(True)
ts,w1,tsC,w1C,F,M_fft,pot_promedio=impulso(Fs,f,amp,N,fase)
plt.xlim(0,N/Fs)
s4.plot(ts,w1,'ro')
s4.plot(tsC,w1C,'b-')


s8 = fig.add_subplot(3,4,8)
plt.title("Transformada")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Amplitud FFT")
s8.grid(True)
plt.xlim(-Fs/2,Fs/2-Fs/N)
plt.ylim(0,1.1)
#s8.plot(F, M_fft)
s8.stem(F, M_fft, 'b', markerfmt=" ", basefmt="-b")


s12 = fig.add_subplot(3,4,12)
plt.xlim(0,10)
plt.ylim(0,10)
plt.axis('off')
s12.spines['right'].set_visible(False)
s12.spines['top'].set_visible(False)
s12.spines['bottom'].set_visible(False)
s12.spines['left'].set_visible(False)
plt.text(0,9,"Frecuencia="+str(f)+"Hz",fontsize=10)
plt.text(0,7,"Amplitud="+str(f'{amp:.{2}f}')+"v",fontsize=10)
plt.text(0,5,"Pot prom="+"NO CORRESPONDE",fontsize=10)
plt.text(0,3,"FrecSampl.="+str(Fs)+"Hz",fontsize=10)
plt.text(0,1,"Muestras="+str(N),fontsize=10)


plt.get_current_fig_manager().window.showMaximized() #para QT5
plt.show()
