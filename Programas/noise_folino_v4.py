# Autor : Pablo D. Folino
# Ejercitación lee información de una archivo "noise_folino_v1.py"

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sc
import os
import simpleaudio as sa
from numpy import genfromtxt
import serial

# Valores iniciales
amp =1.0
header = { "pre": b"*header*", "id": 0, "N": 256, "fs": 10000, "cutFrec":0,"cutFrec2":0,"señal":0,"M":10,"pos":b"end*" }
STREAM_FILE=("/dev/ttyUSB2","serial")
escala=15000

# Valores del enunciado
Fs  =   8000        # Frecuencia de Sampling de la PC
Ts  =   1/Fs        # frecuencia se obtuvieron los datos
fs   =  8000        # Frecuencia de muestreo en Hz, ver frecuencias soportadas de
                    # la place de sonido


def flushStream(f):
    f.flushInput()

def readSamples(adc,synth,N,trigger=False,th=0):
    state="waitLow" if trigger else "sampling"
    i=0
    for t in range(N):
        sample = (readInt4File(streamFile,sign = True)*1.65)/(2**6*512)
        ciaaFFT = (readInt4File(streamFile,sign = True)*1.65)/(2**1*512)
        state,nextI= {
                "waitLow" : lambda sample,i: ("waitHigh",0) if sample<th else ("waitLow" ,0),
                "waitHigh": lambda sample,i: ("sampling",0) if sample>th else ("waitHigh",0),
                "sampling": lambda sample,i: ("sampling",i+1)
                }[state](sample,i)
        adc[i]=sample
        synth[i]=ciaaFFT
        i=nextI

def readInt4File(f,size=2,sign=False):
    raw=f.read(1)
    while( len(raw) < size):
        raw+=f.read(1)
    return (int.from_bytes(raw,"little",signed=sign))

def findHeader(f,h):
    data=bytearray(b'12345678')
    while data!=h["pre"]:
        data+=f.read(1)
        if len(data)>len(h["pre"]):
            del data[0]
    h["id"]      = readInt4File(f,4)
    h["N" ]      = readInt4File(f)
    h["fs"]      = readInt4File(f)
    h["cutFrec"] = readInt4File(f)
    h["cutFrec2"] = readInt4File(f)
    h["señal"]   = readInt4File(f)
    h["M"]       = readInt4File(f)
    data=bytearray(b'1234')
    while data!=h["pos"]:
        data+=f.read(1)
        if len(data)>len(h["pos"]):
            del data[0]
    #print({k:round(v,2) if isinstance(v,float) else v for k,v in h.items()})
    return h["id"],h["N"],h["fs"],h["cutFrec"],h["cutFrec2"],h["señal"],h["M"]


def rx_edu_ciaa():
    global señal_ciaa
    global fData_rx

    flushStream(streamFile)
    reproducir(amp*1.5*tData)
    id,N_rx,fs_rx,cutFrec,cutFrec2,señal,M_rx=findHeader ( streamFile,header )
    nData_rx = np.arange(0,N_rx+M_rx-1,1) #arranco con numeros enteros para evitar errores de float
    adc_rx   = np.zeros(N_rx+M_rx-1)
    ciaaFFT_rx = np.zeros(N_rx+M_rx-1).astype(complex)
    tData_rx=nData_rx/fs_rx
    fData_rx=nData_rx*fs_rx/(N_rx+M_rx-1)-fs_rx/2

    readSamples(adc_rx,ciaaFFT_rx,N_rx+M_rx-1,False,0)

    for i in range(0,N_rx+M_rx-1,1):
        señal_ciaa[i]=0

    conta_rx=0
    while(conta_rx<50):
        conta_rx=conta_rx+1
        señal_ciaa=ciaaFFT_rx+señal_ciaa
        id,N_rx,fs_rx,cutFrec,cutFrec2,señal,M_rx=findHeader ( streamFile,header ) 
        readSamples(adc_rx,ciaaFFT_rx,N_rx+M_rx-1,False,0)  
    
    print(" Se detectaron={} muestras".format(conta_rx))

def reproducir(note):
    audio = note.astype(np.int16)               #tranforma la variable note a entero de 16bits y lo guarda en audio
    play_obj = sa.play_buffer(audio, 1, 2, fs)  # sale el audio
    #play_obj.wait_done()                        # espera que termine la linea anterior

def X_ifft_conFiltro(X_fft,filtro):
    global señal_ciaa
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

def graficar():
    global señal_ciaa,fData_rx

    plt.clf()
#--------------------Señal original-------------------------------------
    s1 = fig.add_subplot(3,2,1)
    plt.title("Señal de entrada x_fft")
    plt.xlabel("Tiempo (seg)")
    plt.ylabel("Amplitud x_fft")
    s1.grid(True)
    s1.plot(t_Data, tData,"g-")

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
    plt.text(0,3,"Muestras="+str(N),fontsize=10)


    s1 = fig.add_subplot(3,2,3)
    plt.title("Transformada de Fouirer (PC)")
    plt.xlabel("Frecuencia en Hz")
    plt.ylabel("Amplitud |X_fft|")
    s1.grid(True)
    s1.plot(f_Data,np.abs(fftData)**2/N,'r-')
#--------------------Señal filtrada-------------------------------------

    s1 = fig.add_subplot(3,2,4)
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
    plt.ion()
    #plt.show()


def valores():
    global N,fs,amp,fase,cantidad,fsec,f,amp2
    os.system("clear")
    print("Los valores actuales son:")
    print("----------------------------------------------------------------------") 
    print("La frecuencia de sampling \t\tfs={}Hz --> \tts={:.4f}mseg".format(fs,1/fs*1000)) 
    print("La amplitud de la señal de red \t\tamp={}%".format(amp*100))
    consulta=input("Desea cambiar los valores S o N[Enter] :")
    if consulta=='S' or consulta =='s':
        valor=input("Ingrese la frecuencia Hz de sampling \t\t\tfs=")
        if valor.isdigit():
            fs=int(valor)
            print("\t\tEl tiempo de sampleo es ts={:.4f}mseg".format(1/fs*1000))
        valor=input("Ingrese la amplitud de la señal de red(0-100)% \t\tamp=")
        if valor.isdigit():
            if(int (valor)>100):
                valor=100
            amp=float(valor)/100  
        input("Presiona cualquier tecla para continuar")
        os.system("clear")


#================================================================
# Inicio del programa principal
#================================================================
global fData_rx,señal_ciaa
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


# Genero el vector de estados temporales
n_Data      = np.arange(0,N,1) #arranco con numeros enteros para evitar errores de float
t_Data      = n_Data*Ts

# Transformada
f_Data=np.arange(-N/2,N/2+M-1)
fftData  = np.fft.fft(tData_ext)
ffthData = np.fft.fft(hData_ext)
YData=fftData*ffthData/escala              # Multiplico en frecuencia

ifftData = np.fft.ifft(YData)               # Antitransformada
fftData  = np.fft.fftshift(fftData)
YData    = np.fft.fftshift(YData)

salida=np.zeros(N)
for i in range(0,N,1):
    salida[i]=ifftData[i]

streamFile = serial.Serial(port=STREAM_FILE[0],baudrate=460800,timeout=None)


# Se grafica para probar las señales
fig = plt.figure(1)
plt.subplots_adjust(left=None, bottom=0.1, right=None, top=0.95, wspace=0.4, hspace=0.8)

# Se guarda espacios para los datos
señal_ciaa=np.zeros(256).astype(complex)
#señal_ciaa=np.zeros(len(t_Data))

menu="""
Programa para Tx señales a la EDU-CIAA por placa de sonido:

[1] Escuchar señal original(TP2)
[2] Escuchar señal filtrada(PC)(TP2)
[3] Visualizar gráfico(PC) 
[4] Rx EDU-CIAA
[5] Escuchar Señal dela EDU-CIAA
[6] Borrar buffer de Rx

[7] TBD
[8] Seteo de frecuencia de la señal de entrada, número de muestras, 
    frecuencia de sampling.
[9] Salir
"""

while(True):
    os.system("clear")
    print(menu)

    opcion=input("Elija una opción: ")

    if opcion== '1':
        # Se reproduce la señal original
        reproducir(amp*2*tData)
    elif opcion== '2':
        reproducir((np.real(amp*2*salida)).astype(np.int16) )
    elif opcion== '3':
        graficar()  
    elif opcion== '4':
        rx_edu_ciaa()
    elif opcion== '5':
        pass
    elif opcion== '6':
        flushStream(streamFile)
    elif opcion== '7':
        pass
    elif opcion== '8':
        valores()
    elif opcion== '9':
        os.system("clear")
        print("Gracias por usar el programa !!!")
        exit (0)
    else:
        print("No selecionó una opción válida\n\r")

