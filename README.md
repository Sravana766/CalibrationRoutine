# CalibrationRoutine

## Introduction
This Readme provides a breif explanation and an insight into this routine. Please refer Accoustic Calibration Manual in the repository for detailed steps. 

## Calibration 
It is essential to calibrate a sound device to ensure that the device is reliable. This process also ensures that the measurements of the device are accurate and can be trusted. If not, the team can take steps to compensate for these errors. 

The calibration process described in this paper measures the frequency response of an input stimulus played by the speakers. The frequency response of an audio system represents the range of frequencies that the system can produce. It measures if the audio system reproduces the sound in all frequencies or if the system caused changes to the input stimulus. It also measures the sound level variations caused at each frequency. Ideally, one would expect a flat frequency response. Meaning that the audio system does not affect the frequency of the sound it is playing. However, this is not the case in many audio systems such as speakers, amplifiers, headphones, and microphones. There are some variations in the sound level caused by the audio system at different frequencies.  

After measuring the frequency response, this process then uses inverse filtering techniques to recompensate for the variations that occurred at various frequencies. This way, the sound played by the audio system is consistent and reliable. 

The hardware utilized in this process is KEMAR and NI-DAQ. KEMAR is a “head and torso simulator for acoustic measurements” (KEMAR Website). It is a mannequin that contains a microphone and a preamplifier inside the ears. The audio system plays the input stimulus from the computer to this mannequin. This mannequin sends the altered frequency back to the computer. NI-DAQ is a data acquisition card that records and digitizes audio before sending them to the computer. The hardware used, are compatible with the software. 
  
The software utilized in this process is Python and MATLAB. There are many libraries available in Python that allows for convenient interface between the hardware and software. The Python script used in this process utilizes the EasyDAQmx from the Python library to read data from the NI-DAQ. The intermediary interface between the hardware and the Python script is NI-DAQmx. This software reads the information from the NI-DAQ hardware. The Python script is responsible for playing the input stimulus and storing the altered frequency output collected from the KEMAR mannequin. The process uses MATLAB to analyze the frequency response of the KEMAR output and performing the inverse filtering later.
It is important to note that the calibration process described in this paper is intended for research and scientists. It is not available for mainstream purposes such as calibrating a musical device. This process is a possible alternative to the calibration techniques performed today. Due to the low cost of the hardware used in this process, it is a cheaper alternative. As mentioned earlier, although the technique is not new, the software is not always available for researchers or scientists. This software provides a cheaper and available alternative of sound calibration to researchers.

## iPad Accoustical Calibration 

The current system involves an iPad connected to a pair of headphones. The iPad will play pure tones in an automated fashion. Depending on the user input, the algorithm will either increase or decrease the pure tone's sound intensity on the next iteration. The data from each trial will then be analyzed. Hence it is essential to ensure that the user hears the sound at the same intensity that the algorithm sets. 

![system](https://user-images.githubusercontent.com/62814852/153280509-6c3df0de-dc41-4b33-9d21-2170b718e428.PNG)


The figure above describes system. The unity code is the input to the system that is responsible for playing the sounds. The output is the voltage recordings that the computer recieves. 

The calibration process first involved scaling the sound intensity measurement to a dB SPL scale by dividing the signal by its own RMS vale. This property allows for scaling and manipulating the sound played from the iPad.

However, the exact value of the dB SPL scale was unknown. This value was measured using a sound level meter. The sound level meter measured the sound level of output sound from the headphones. The system was then set up in a soundproof booth to ensure that outside noise does not affect this measurement. The equation below was then used to manipulate the sounds and play them at a desired dB SPL level.

The next step is playing the calibrated sounds from the iPad to KEMAR. KEMAR is a head torso stimulator that is used for acoustic research. The calibrated pure tone with a frequency of 1kHz was played to KEMAR via headphones. The computer recorded the output response from KEMAR. 

The output signal was analyzed using a transfer function, and maximum length sequence. MLS is a useful technique used in signal processing for measuring the impulse response of a linear system using the MLS sequence as an input to the system. The MLS sequence itself is a pseudorandom signal that consists of values from -1 to 1. The pseudorandom samples of the maximum length sequence were collected from the output response from the KEMAR. 
![image](https://user-images.githubusercontent.com/62814852/153279346-3c384d41-dcf9-4a2c-9a48-51cee72db1e0.png)

Figure 1: MLS sequence for 1kHz pure tone

The above picture is the MLS sequence of the 1kHz pure tone. The system's input is the MLS sequence, and the impulse response is the output of the system. The impulse response of the system is in the frequency domain. The system uses Fourier transform to transform the function into the frequency domain. 

![image](https://user-images.githubusercontent.com/62814852/153279398-9f618280-daef-4791-aa91-05440a0339bc.png)

Figure 2: Impulse response of the system to an input of 1kHz pure tone

The above picture is the impulse response of playing a 1kHz sound for a certain period. The impulse response is expected to have a flat frequency response.  However, in the picture it can be observed that the frequency response is not flat. For instance, the amplitude at 2kHz frequency and 5 kHz frequency are not the same. It is expected that both frequencies have the same magnitude. The transfer function applied in the system helps identify which frequencies need to be compensated to have the desired magnitude. This compensation is done in the original iPad algorithm. 




